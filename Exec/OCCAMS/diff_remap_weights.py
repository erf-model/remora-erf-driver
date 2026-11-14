#!/usr/bin/env python3
"""Quantify differences between two remap_weights/ dumps (driver.remap_dump_weights).

check_remap_weights.sh compares the dumped arrays byte for byte. That is the
right test for one binary under different settings, but two builds that round
the host-side clipping differently (the Perlmutter Cray CPU and CUDA builds,
#34) never match exactly. This script says by how much they differ.

Usage:
  diff_remap_weights.py <dir-a> <dir-b> [--tol T]

<dir-a> and <dir-b> are remap_weights/ directories, or case directories that
contain one, or check_remap_weights.sh result directories (every case found in
both is compared).

For each weight/index pair (a2o_weight with a2o_index, and o2a_weight with
o2a_index) it reports:
  slot      max |w_a - w_b| slot by slot, as the exchanges read them
  same-src  max |w_a - w_b| per destination cell, matching entries by source
            cell rather than by slot. Edge-touching overlaps whose area is pure
            rounding (~1e-15 of the cell) can appear in one build and not the
            other, which shifts later entries into different slots; "slot" then
            reports differences of order one that are not weight differences.
  only-a/b  entries present in one build only, with the largest such weight
  reordered destination cells whose source list differs in order or content
All other arrays (valid flags, uncoupled mask, donor flags) are compared exactly.

Exit status is 1 when a same-src difference exceeds --tol (default 1e-12), an
only-a/b weight exceeds --tol, a non-weight array differs, or arrays are
missing or have different shapes; 0 otherwise.
"""
import argparse
import os
import sys

import numpy as np

DTYPES = {"float64": np.float64, "float32": np.float32, "int32": np.int32, "int64": np.int64}


def read_array(d, name):
    meta = {}
    with open(os.path.join(d, name + ".txt")) as f:
        for line in f:
            key, _, val = line.partition(" ")
            meta[key] = val.strip()
    ncomp = int(meta["ncomp"])
    data = np.fromfile(os.path.join(d, name + ".bin"), dtype=DTYPES[meta["type"]])
    return data.reshape(ncomp, -1), meta["box"]


def stencils(W, I):
    # Slot s holds weight W[s] and source index I[3s], I[3s+1] (I[3s+2] is k = 0).
    nslot, ncell = W.shape
    out = []
    for c in range(ncell):
        out.append([((int(I[3 * s, c]), int(I[3 * s + 1, c])), W[s, c])
                    for s in range(nslot) if W[s, c] != 0.0])
    return out


def compare_dir(a, b, tol):
    names_a = {f[:-4] for f in os.listdir(a) if f.endswith(".txt")}
    names_b = {f[:-4] for f in os.listdir(b) if f.endswith(".txt")}
    bad = False
    for n in sorted(names_a ^ names_b):
        print(f"  {n}: only in {'a' if n in names_a else 'b'}")
        bad = True
    names = names_a & names_b

    done = set()
    for w in sorted(n for n in names if "_weight" in n):
        i = w.replace("_weight", "_index")
        if i not in names:
            continue
        done |= {w, i}
        (Wa, box_a), (Wb, box_b) = read_array(a, w), read_array(b, w)
        Ia, _ = read_array(a, i)
        Ib, _ = read_array(b, i)
        if Wa.shape != Wb.shape or Ia.shape != Ib.shape:
            print(f"  {w}: SHAPE MISMATCH {box_a} vs {box_b}")
            bad = True
            continue
        dW = np.abs(Wa - Wb)
        slot = float(dW.max())
        # Cells whose index lists match slot for slot need no matching by source.
        nslot = Wa.shape[0]
        Ia2 = Ia.reshape(nslot, 3, -1)[:, :2]
        Ib2 = Ib.reshape(nslot, 3, -1)[:, :2]
        live_a, live_b = Wa != 0.0, Wb != 0.0
        idx_same = (live_a == live_b).all(axis=0) & \
                   ((Ia2 == Ib2).all(axis=1) | ~live_a).all(axis=0)
        same = float(dW[:, idx_same].max(initial=0.0))
        cells = np.flatnonzero(~idx_same)
        sa, sb = stencils(Wa[:, cells], Ia[:, cells]), stencils(Wb[:, cells], Ib[:, cells])
        only_a, only_b = [], []
        reordered = 0
        for ca, cb in zip(sa, sb):
            if [k for k, _ in ca] != [k for k, _ in cb]:
                reordered += 1
            da, db = dict(ca), dict(cb)
            for k in da.keys() & db.keys():
                same = max(same, abs(da[k] - db[k]))
            only_a += [da[k] for k in da.keys() - db.keys()]
            only_b += [db[k] for k in db.keys() - da.keys()]
        max_oa = max(only_a, default=0.0)
        max_ob = max(only_b, default=0.0)
        flag = "" if max(same, max_oa, max_ob) <= tol else "   <-- exceeds tol"
        bad |= bool(flag)
        print(f"  {w:20s} cells={Wa.shape[1]:8d} slot={slot:.3e} same-src={same:.3e} "
              f"only-a={len(only_a)} (max {max_oa:.1e}) only-b={len(only_b)} (max {max_ob:.1e}) "
              f"reordered={reordered}{flag}")

    for n in sorted(names - done):
        (Xa, box_a), (Xb, box_b) = read_array(a, n), read_array(b, n)
        if Xa.shape != Xb.shape:
            print(f"  {n}: SHAPE MISMATCH {box_a} vs {box_b}")
            bad = True
            continue
        ndiff = int((Xa != Xb).sum())
        print(f"  {n:20s} {'identical' if ndiff == 0 else f'{ndiff} entries DIFFER'}")
        bad |= ndiff > 0
    return bad


def weight_dirs(top):
    """Map case name -> remap_weights dir under top."""
    if os.path.isfile(os.path.join(top, "a2o_weight.txt")):
        return {".": top}
    if os.path.isdir(os.path.join(top, "remap_weights")):
        return {".": os.path.join(top, "remap_weights")}
    return {c: os.path.join(top, c, "remap_weights") for c in sorted(os.listdir(top))
            if os.path.isdir(os.path.join(top, c, "remap_weights"))}


def main():
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("a")
    p.add_argument("b")
    p.add_argument("--tol", type=float, default=1e-12)
    args = p.parse_args()

    da, db = weight_dirs(args.a), weight_dirs(args.b)
    if not da:
        sys.exit(f"error: no remap_weights found under {args.a}")
    bad = False
    for case in sorted(da.keys() | db.keys()):
        if case not in da or case not in db:
            print(f"{case}: only under {'a' if case in da else 'b'}")
            bad = True
            continue
        print(f"{case}:")
        bad |= compare_dir(da[case], db[case], args.tol)
    print("RESULT: " + ("DIFFERENCES ABOVE TOLERANCE" if bad else f"all within {args.tol:g}"))
    sys.exit(1 if bad else 0)


if __name__ == "__main__":
    main()

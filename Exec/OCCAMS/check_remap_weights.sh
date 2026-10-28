#!/usr/bin/env bash
#
# Before/after checks for the spherical remap-weight build (issue #31: make the
# host-side weight construction safe when the default arena is device memory).
#
# The Upwelling decks only reach the Cartesian remap. The code that issue #31
# changes in the spherical path -- the ERF corner mesh and its ring extension,
# both staging plans, donor fill, the O2A owner match -- only runs on a case
# with lon/lat on both grids, which is OCCAMS. This script runs a fixed set of
# short OCCAMS cases and compares two sets of results.
#
# Usage:
#   check_remap_weights.sh run <out-dir> [case ...]
#   check_remap_weights.sh compare <baseline-dir> <candidate-dir>
#
# Cases (default: smoke smoke_ring smoke_coawst twoway):
#   smoke          inputs_OCCAMS_smoke: one-way, spherical. ERF corner mesh, ERF
#                  staging plan, A2O weights and normalization.
#   smoke_ring     smoke + driver.remap_extend_corner_ring=1 (off by default).
#   smoke_coawst   smoke + driver.remap_donor_fill=coawst: the donor search.
#   twoway         inputs_OCCAMS_state at the smoke cadence (2 ERF steps per
#                  ocean step): REMORA staging plan, O2A weights, owner match.
#   twoway_timing  inputs_OCCAMS_state at its own cadence (30 ERF steps), no
#                  plotfiles. For timing managed vs device memory; not compared.
#
# Environment:
#   DRIVER      erf-remora-driver to run (required)
#   FCOMPARE    amrex_fcompare (default: next to DRIVER, which needs the build
#               configured with -DERF_REMORA_ENABLE_FCOMPARE=ON)
#   LAUNCH      launcher prefix (default: "srun -N 1 -n 16")
#   WRAP        placed between LAUNCH and DRIVER, e.g. "compute-sanitizer --tool memcheck"
#   EXTRA_ARGS  key=value overrides appended to every case. Pass
#               amrex.the_arena_is_managed explicitly: the decks set it to 1
#               before the fix and not at all after it.
#   FCOMPARE_TOL  extra fcompare arguments for compare (default: none, i.e.
#               exact). GPU runs are not bitwise reproducible past the first
#               exchange, so GPU comparisons need e.g.
#               "--rel_tol 1e-9 --abs_tol 1e-9". The dumped remap arrays
#               (remap_weights/*.bin) are always compared exactly.
#
# The script reads the decks (inputs_OCCAMS_*, inputs_ERF, inputs_REMORA) from
# its own directory, so it can live in a copy of Exec/OCCAMS anywhere, e.g. on
# scratch. The erfbdy cache is kept in that directory too.
#
# Results: <out-dir>/<case>/{run.log,erfplt*,remoraplt*,exchange_*}, plus
# <out-dir>/MANIFEST recording the driver, launcher and overrides. Each case
# directory holds symlinks to the decks it ran.

set -uo pipefail

here=$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)

# name|deck|case overrides separated by ';'
case_table=(
  "smoke|inputs_OCCAMS_smoke|"
  "smoke_ring|inputs_OCCAMS_smoke|driver.remap_extend_corner_ring=1"
  "smoke_coawst|inputs_OCCAMS_smoke|driver.remap_donor_fill=coawst"
  "twoway|inputs_OCCAMS_state|driver.remap_frame=spherical;erf.fixed_dt=10.0;remora.fixed_dt=20.0;remora.ndtfast=20"
  "twoway_timing|inputs_OCCAMS_state|driver.remap_frame=spherical;erf.plot_int_1=-1;remora.plot_int=-1;driver.plot_int_exchange=-1"
)
default_cases=(smoke smoke_ring smoke_coawst twoway)

# Every compared case writes its exchange plotfiles at each exchange and the
# solver states at start and end, and prints the full remap diagnostics. The
# erfbdy cache is shared with run_occams.sh so it is not rebuilt per case.
common_args=(
  driver.remap_v=2
  driver.remap_dump_weights=1
  driver.plot_int_exchange=1
  erf.plot_int_1=2 remora.plot_int=1
  erf.check_int=-1 remora.check_int=-1
  erf.plot_file_1=erfplt remora.plot_file=remoraplt
  driver.plot_file_exchange_a2o=exchange_a2o
  driver.plot_file_exchange_o2a=exchange_o2a
  "erf.erfbdy_file=${here}/erfbdy"
)

# The lines that describe the weights: frame, staging, mesh, coverage, owner
# match, stencil counts, SCRIP comparison and warnings. Timings are dropped.
pick_remap_lines () {
  grep -E '^\[driver\]|coverage|conservation|WARNING' "$1" | grep -viE 'time|elapsed|seconds|weights written'
}

run_cases () {
  local out=$1; shift
  local cases=("$@")
  [[ ${#cases[@]} -eq 0 ]] && cases=("${default_cases[@]}")

  : "${DRIVER:?set DRIVER to the erf-remora-driver to run}"
  local driver; driver=$(realpath "$DRIVER")
  [[ -x $driver ]] || { echo "error: no executable at $driver" >&2; exit 1; }
  local launch=${LAUNCH:-srun -N 1 -n 16}
  local wrap=${WRAP:-}
  local extra=${EXTRA_ARGS:-}

  if [[ -f $here/erfbdy/Header && $here/erfbdy/Header -ot $driver ]]; then
    echo "WARNING: $here/erfbdy predates the driver; move it aside if ERF changed (see run_occams.sh)." >&2
  fi

  mkdir -p "$out"
  out=$(realpath "$out")
  {
    echo "date      : $(date '+%Y-%m-%d %H:%M:%S')"
    echo "host      : $(hostname)"
    echo "driver    : $driver"
    echo "built     : $(date -r "$driver" '+%Y-%m-%d %H:%M:%S')"
    echo "launch    : $launch"
    echo "wrap      : ${wrap:-none}"
    echo "extra     : ${extra:-none}"
    echo "CUDA_LAUNCH_BLOCKING=${CUDA_LAUNCH_BLOCKING:-unset}"
  } >> "$out/MANIFEST"

  local failed=()
  for name in "${cases[@]}"; do
    local entry="" e
    for e in "${case_table[@]}"; do [[ ${e%%|*} == "$name" ]] && entry=$e; done
    [[ -n $entry ]] || { echo "error: unknown case $name" >&2; exit 1; }

    local deck args_str
    IFS='|' read -r _ deck args_str <<< "$entry"
    local case_args=()
    [[ -n $args_str ]] && IFS=';' read -r -a case_args <<< "$args_str"

    local dir=$out/$name
    if [[ -e $dir ]]; then echo "error: $dir exists; refusing to overwrite" >&2; exit 1; fi
    mkdir -p "$dir"
    # The decks pull in inputs_ERF and inputs_REMORA with FILE = <name>, which
    # ParmParse opens relative to the working directory, not to the deck. Link
    # all three into the case directory and run from there.
    local f
    for f in "$deck" inputs_ERF inputs_REMORA; do
      [[ -f $here/$f ]] || { echo "error: $here/$f not found; keep it next to this script" >&2; exit 1; }
      ln -s "$here/$f" "$dir/$f"
    done

    echo "=== $name: $deck ${case_args[*]:-} ${extra}"
    local start=$SECONDS
    # ParmParse keeps the last value given, so case overrides follow the common
    # arguments and EXTRA_ARGS follows both.
    # shellcheck disable=SC2086  # LAUNCH, WRAP and EXTRA_ARGS are word lists
    if (cd "$dir" && $launch $wrap "$driver" "$deck" \
          "${common_args[@]}" ${case_args[@]+"${case_args[@]}"} $extra > run.log 2>&1); then
      echo "    ok ($((SECONDS - start)) s)"
    else
      echo "    FAILED ($((SECONDS - start)) s); see $dir/run.log"
      grep -m3 -E 'Abort|Error|SIGSEGV|SIGABRT|error' "$dir/run.log" | sed 's/^/    /'
      failed+=("$name")
    fi
  done
  [[ ${#failed[@]} -eq 0 ]] || { echo "failed: ${failed[*]}"; return 1; }
}

compare_cases () {
  local base=$1 cand=$2
  local fcompare=${FCOMPARE:-$(dirname "$(realpath "${DRIVER:-$base}")")/amrex-build/Tools/Plotfile/amrex_fcompare}
  [[ -x $fcompare ]] || { echo "error: no fcompare at $fcompare; set FCOMPARE" >&2; exit 1; }
  local tol=${FCOMPARE_TOL:-}
  local fail=0

  # fcompare exits 1 both on a real difference and when it cannot do its job, so
  # first make sure it runs at all: a plotfile compared with itself must pass.
  # A build made with the Cray wrappers and cudatoolkit loaded links libcudart
  # even for CPU, and then needs that module environment to start.
  local probe
  probe=$(ls -d "$base"/*/*/Header 2>/dev/null | head -1)
  if [[ -n $probe ]]; then
    probe=$(dirname "$probe")
    local probe_out
    if ! probe_out=$("$fcompare" "$probe" "$probe" 2>&1); then
      echo "error: fcompare cannot run ($fcompare):" >&2
      echo "$probe_out" | head -3 | sed 's/^/    /' >&2
      echo "       load the module environment the build used (e.g. the Perlmutter profile) and re-run compare" >&2
      exit 1
    fi
  fi

  for bdir in "$base"/*/; do
    local name; name=$(basename "$bdir")
    [[ $name == twoway_timing ]] && continue
    local cdir=$cand/$name
    [[ -d $cdir ]] || { echo "$name: missing in candidate"; fail=1; continue; }

    local n=0 bplt plt
    for bplt in "$bdir"/*/Header; do
      plt=$(basename "$(dirname "$bplt")")
      n=$((n + 1))
      if [[ ! -f $cdir/$plt/Header ]]; then echo "$name/$plt: missing"; fail=1; continue; fi
      # shellcheck disable=SC2086
      if ! "$fcompare" --abort_if_not_all_found $tol "$bdir/$plt" "$cdir/$plt" > "$cdir/$plt.fcompare" 2>&1; then
        # A real difference comes with fcompare's error table; anything else is
        # fcompare failing (missing variable, unreadable plotfile, ...).
        if grep -q "variable name" "$cdir/$plt.fcompare"; then
          echo "$name/$plt: DIFFERS (see $cdir/$plt.fcompare)"
        else
          echo "$name/$plt: fcompare FAILED, not a comparison (see $cdir/$plt.fcompare)"
        fi
        fail=1
      fi
    done
    # The remap arrays themselves (driver.remap_dump_weights): always exact, since
    # they are built on the host whatever the GPU does downstream. Binaries that
    # predate the option write none, and the check is skipped.
    local nw=0 bw
    if [[ -d $bdir/remap_weights && -d $cdir/remap_weights ]]; then
      for bw in "$bdir"/remap_weights/*.bin; do
        nw=$((nw + 1))
        if ! cmp -s "$bw" "$cdir/remap_weights/$(basename "$bw")"; then
          echo "$name/remap_weights/$(basename "$bw"): DIFFERS"; fail=1
        fi
      done
    fi
    if ! diff <(pick_remap_lines "$bdir/run.log") <(pick_remap_lines "$cdir/run.log") > "$cdir/remap_lines.diff"; then
      echo "$name: remap diagnostics differ (see $cdir/remap_lines.diff)"; fail=1
    fi
    echo "$name: compared $n plotfiles, $nw weight arrays and $(pick_remap_lines "$bdir/run.log" | wc -l) remap lines"
  done

  if (( fail )); then echo "RESULT: DIFFERENCES FOUND"; return 1; fi
  echo "RESULT: identical${tol:+ (within $tol)}"
}

case ${1:-} in
  run)     shift; [[ $# -ge 1 ]] || { echo "usage: $0 run <out-dir> [case ...]" >&2; exit 2; }
           run_cases "$@" ;;
  compare) shift; [[ $# -eq 2 ]] || { echo "usage: $0 compare <baseline-dir> <candidate-dir>" >&2; exit 2; }
           compare_cases "$@" ;;
  *)       sed -n '2,40p' "${BASH_SOURCE[0]}"; exit 2 ;;
esac

# Exec Reproducer Index

This directory contains the input decks for the ERF-REMORA driver, plus a guide
to reading the driver's own diagnostic output.

Each case is split into three files: a coupled deck that sets the driver
controls, and the two solver decks it pulls in with `FILE =`. Run the coupled
deck, not the solver decks.

## Inputs

### The Upwelling case

Both Upwelling decks are the coastal upwelling problem of Perlin et al. (2007),
<https://doi.org/10.1175/JPO3104.1>, used for validation as well as for MPI
regression. They descend from `inputs_perlin_ERF_REMORA`, which was split into
the one-way and two-way decks below.

### `Upwelling/inputs_oneway`

Purpose:

- `ERFtoREMORA` idealized upwelling. ERF forces REMORA; no SST returns to the
  atmosphere.

Run:

```bash
cd Exec/Upwelling
../../build/erf-remora-driver inputs_oneway
```

Pulls in `inputs_ERF` and `inputs_REMORA` from the same directory. Needs no
external data, so this is the case to run first on a new build.

---

### `Upwelling/inputs_twoway`

Purpose:

- `two_way` idealized upwelling. REMORA SST feeds ERF's surface layer through
  the coupled-SST path (`erf.use_coupled_sst`, set by the driver).

Run:

```bash
cd Exec/Upwelling
../../build/erf-remora-driver inputs_twoway
```

What to verify:

- SST arrives in ERF and is physical, not zero and not a placeholder,
- coupling clock and ocean subcycling advance together.

Coverage caveat: the Upwelling grids are matched 1:1 and fully wet, so this deck
reports complete coverage by construction. It is a no-regression check on the
remap, **not** a test of coverage, masking, or the coverage threshold — none of
those code paths are exercised here. Use OCCAMS for that.

---

### `OCCAMS/inputs_OCCAMS_oneway`

Purpose:

- `ERFtoREMORA` on the real OCCAMS domain. Atmosphere forces ocean only.

Run:

```bash
cd Exec/OCCAMS
srun -N 1 -n 16 ../../build/erf-remora-driver inputs_OCCAMS_oneway
```

---

### `OCCAMS/inputs_OCCAMS_state`

Purpose:

- `two_way` OCCAMS with the **state-passing** atmosphere-to-ocean contract.
  Exercises the ocean-to-atmosphere SST return path in the spherical frame.

Run:

```bash
cd Exec/OCCAMS
srun -N 1 -n 16 ../../build/erf-remora-driver inputs_OCCAMS_state
```

---

### `OCCAMS/inputs_OCCAMS_flux`

Purpose:

- `two_way` OCCAMS with the **flux-passing** contract. Same SST return path as
  `inputs_OCCAMS_state`, with the surface-flux lanes on top.

Run:

```bash
cd Exec/OCCAMS
srun -N 1 -n 16 ../../build/erf-remora-driver inputs_OCCAMS_flux
```

### `OCCAMS/inputs_OCCAMS_smoke`

Purpose:

- Fast guard on everything expensive to reach: ERF WRF-init, REMORA NetCDF init
  (grid, masks, boundaries, forcings), spherical remap weight construction, the
  ERF-to-REMORA staging copy, and one full pack/apply exchange. One ocean step,
  two ERF steps — minutes rather than a queue slot, so it can be rerun after
  every pull.

Run:

```bash
cd Exec/OCCAMS
srun -N 1 -n 16 ../../build_ncutil/erf-remora-driver inputs_OCCAMS_smoke
```

Unlike the production decks it sets `amrex.fpe_trap_invalid = 1`, so floating
point problems abort instead of propagating. It also pins
`driver.remap_frame = spherical` and `driver.remap_v = 1`, and writes no
plotfiles or checkpoints — the coupling diagnostics are the pass/fail signal.

What it does not do: two ERF steps are not enough to expose dynamical
instability, so a blowup that appears at step 5 of a production run will not
show up here. Use `inputs_OCCAMS_oneway` for that.

Status: committed but not yet run. Treat the first execution as a test of the
deck as much as of the code.

---

#### OCCAMS notes, common to all four

Deck naming: the two-way decks are `inputs_OCCAMS_state` and
`inputs_OCCAMS_flux`, distinguished by contract. There is no
`inputs_OCCAMS_twoway`. `inputs_OCCAMS_smoke` and `inputs_OCCAMS_oneway` are
both `ERFtoREMORA`.

External data: these decks reference absolute NERSC paths for `wrfinput_d01`,
`wrfbdy_d01`, the REMORA initial/grid/boundary files, and the forcing file.
None of it is in the repository, so the decks are not runnable elsewhere without
editing those paths.

Geometry: ERF runs `init_type = "WRFInput"` and the REMORA grid file carries
`lon_psi`/`lat_psi`, so the remap resolves to the **spherical** frame. Do not
force `driver.remap_frame = cartesian` here. The `prob_lo`/`prob_hi` in these
decks are REMORA's internal Cartesian geometry only; they are not a lon/lat
extent, and remapping against them silently produces a wrong overlap.

Grids are non-conformal by design — ERF is `637 x 754`, REMORA is `563 x 751`
over the same extent — so partial coverage is expected and the coverage
diagnostics below carry real information.

---

## Reading driver output

Driver diagnostics are tagged `[driver]` or `[coupling-layout]`. Solver output
from ERF and REMORA is untagged.

Most of these lines are not pass/fail. Several report degradation that does not
stop the run and will not appear anywhere else, so the useful question for each
is not "what does it say" but "what value is normal here".

### Verbosity

- `driver.v = 1` — driver lifecycle and layout queries.
- `driver.remap_v = 1` — remap geometry detail: corner ranges, extents, stencil
  occupancy, and per-direction coverage summaries.

Warnings and the coverage summary print regardless. Raise `remap_v` when a
coverage number looks wrong and you need the geometry behind it.

### What contract is running

```
Atmos->ocean contract mode: <name> (<N> active layers)
Atmos->ocean active lanes: <names>
```

Under the flux-passing contract the lane set changes after the first ocean step,
which is reported separately:

```
Atmos->ocean bootstrap lanes (ocean_step=0): <names>
Atmos->ocean steady-state lanes (ocean_step>=1): <names>
```

A lane can also be deliberately withheld, reported per lane as
`<lane>: withheld (not supplied by ERF)`. Withheld is not a failure: it means
the driver declined to pass a value it does not have, so REMORA keeps whatever
its own configuration provides. The two common cases announce themselves:

- `[driver] ERF reports no active radiation model, so the SWrad/LWrad
  atmos->ocean lanes are withheld rather than passed as zeros.`
- `[driver] ERF has no cloud condensate ... so the Cloud atmos->ocean lane is
  withheld ...`

Both are correct behaviour when ERF has no radiation or moisture model. They
matter because the alternative — passing zeros — would look like real forcing.

### Which geometry the remap chose

```
[driver] BuildRemapWeights: remap frame = spherical (requested auto)
```

On OCCAMS this must say `spherical`. If it says `cartesian`, detection of
`lon_psi`/`lat_psi` failed; the run will complete and the answer will be wrong.
This is the single most important line in the startup output.

Under `remap_v = 1` the frame is followed by the geometry it derived — corner
ranges in the spherical case, `remap Cartesian extents` otherwise, and:

```
[driver] remap ERF corner mesh: nodes=... max corner shift vs raw edge arrays=... deg
```

`max corner shift` is how far the reconstructed staggered corners moved relative
to the raw `XLONG_U`/`XLAT_V` arrays. Small is expected. Large means the corner
reconstruction disagrees with the file's own staggered coordinates.

Longitude convention is reconciled automatically:

```
[driver] remap: shifting ERF longitudes by <n> deg to match the REMORA convention
```

That is routine — one grid is on `[-180,180]` and the other on `[0,360]`. The
failure form is not:

```
[driver] WARNING remap: ERF longitudes [...] and REMORA longitudes [...] do not
overlap under any 360-degree shift.
```

That means the two grids do not describe a common region at all, and every
coverage count below will be zero.

### Coverage, ocean to atmosphere

```
[driver] O2A coverage: uncovered <A> outside ocean footprint, <B> source-masked,
                       <C> below coverage threshold (0.001)
```

Three distinct reasons an ERF cell gets no ocean SST, and they mean different
things:

| Bucket | Meaning | Normal on OCCAMS | If it is large |
|---|---|---|---|
| `outside ocean footprint` | No REMORA cell overlaps this ERF cell at all. | **Large.** Roughly the atmosphere area not over the ocean nest. | Expected; not a signal. |
| `source-masked` | REMORA cells overlap, but all are dry in REMORA's mask. | **Thin, coastline-shaped.** | Grids misaligned — the ocean land mask is landing in the wrong place. This is the number that catches a bad spherical corner reconstruction. |
| `below coverage threshold` | Real overlap, too little of it. | **Near zero.** Slivers at the nest edge. | The two grids are offset by a fraction of a cell. |

Uncovered cells are not a hole: ERF keeps its own `wrflowinp` SST there. That is
the withholding contract, and it is why the count being large is informative
rather than fatal.

The threshold is `driver.remap_min_coverage_frac`, default `0.001`. A sliver
overlap is worse than no overlap: normalization would turn 1e-6 of coverage into
a confident-looking answer. Discarding negligible destination fractions is
ordinary practice in conservative remapping — see Jones (1999) for the scheme
this remap follows, and the CESM and OASIS3-MCT couplers for production
implementations.

The specific value is inherited from the driver's source comment, which cites
CESM at a destination fraction of 0.001 and OASIS/CERFACS at 1/1000. Those two
numbers have not been checked against either code's source; they are plausible
and conventional, but treat them as unverified until someone reads the couplers
directly.

### Coverage, atmosphere to ocean

```
[driver] WARNING remap (<label>): <N> of <M> REMORA wet cells have no ERF source
overlap.
```

Those cells fall back to standard-atmosphere placeholders rather than coupled
values — forced, but not coupled. Expect a small count along coastlines where
the two land masks disagree, or where the ocean grid reaches past the
atmosphere. A large count means the atmosphere does not cover the ocean nest.

Under `remap_v = 1`, the per-direction summary:

```
[driver] remap coverage (<label>): cells=... uncovered=... no_stencil=...
                                   frac min/mean/max=.../.../... max|frac-1|=...
```

`max|frac-1|` is the conservation check. The remapped weights should sum to one
over each destination cell; a large deviation means the weights are not
conservative and the coupled fluxes are not either.

### Warnings that degrade silently

These do not stop the run and nothing downstream reports them:

```
[driver] remap: <N> wet destination cells had no valid source overlap and were
filled from their nearest non-land source cell (COAWST check_weights behaviour).
```

By design, and it matches COAWST (Warner et al. 2010), the coupled system this
driver's remap was ported from. Those cells carry a neighbour's atmospheric
state rather than their own. Watch the trend: a count that *rises* between runs
means the two land masks are drifting apart.

```
[driver] WARNING BuildOceanToAtmosWeights: <N> ocean overlaps dropped because a
destination cell exceeded MAX_STENCIL_SIZE=<K> (max overlaps seen = <M>);
weights normalized over the largest-area overlaps kept.
```

The stencil overflowed. The run continues with weights renormalized over the
largest-area overlaps retained, so conservation is degraded — silently. It
implies one destination cell covers more source cells than the fixed stencil
holds, i.e. a resolution ratio larger than the remap was built for.

The non-warning form under `remap_v = 1` is the one to check proactively:

```
[driver] remap stencils: max overlaps per destination cell = <M> of <K>
```

If `<M>` is near `<K>`, the next grid change will overflow.

### Parallel layout

```
[driver] remap staging plan: remora_boxes=... erf_boxes=... staged_cells=...
         (erf domain ...) max_box_cells=... remora_boxes_without_overlap=...
[driver] O2A remap staging plan: erf_boxes=... remora_boxes=... staged_cells=...
```

The two solvers own different `DistributionMapping`s, so remap inputs are staged
across them. `staged_cells` far exceeding the domain cell count means boxes are
being replicated more than the decomposition requires. `[coupling-layout]` lines
under `driver.v = 1` dump the BoxArrays and process maps behind this.

The staging path asserts three assumptions about solver layouts: that REMORA psi
boxes are nodal conversions of the rho mask boxes, that `lon_psi` and `lat_psi`
share box indexing and processor mapping, and that ERF corners are nodal
conversions of the O2A weight boxes. If a real layout violates one, the run
aborts on the assertion rather than mis-mapping quietly.

## Tests

The CTest suite drives the Upwelling decks with runtime overrides — 74 tests,
all from `Exec/Upwelling`. OCCAMS is not in the suite; its external data is not
in the repository.

### Run the suite at `-j1`

```bash
ctest --test-dir build --output-on-failure
```

**The suite is only valid at `-j1`.** The harnesses `file(REMOVE_RECURSE)` their
working directories, so concurrent tests delete each other's output and fail in
ways that look like solver bugs — a stale plotfile, a missing checkpoint, a
comparison against a directory that vanished mid-run. Two suites run at once
produced 49 phantom failures; `-j2` alone produced 5. Do not add `-j` here, and
do not start a second suite in another worktree while one is running.

### Test name structure

Names are generated from the coupling axes, so you can select any slice with
`-R`:

```
reproducer_restart_<direction>_<contract>_<cadence>_nopbl
parallel_upwelling_<direction>_<contract>_<cadence>_nopbl_np<N>_<kind>
```

- `<direction>` — `oneway` or `twoway`
- `<contract>` — `state` or `flux`
- `<cadence>` — `matched` (one ocean step per atmosphere step) or `substep`
- `<N>` — `2`, `4`, `6`, or `8`
- `<kind>` — `smoke` (does it run) or `compare_to_np1` (does it match serial)

### Running one test

A single restart test:

```bash
ctest --test-dir build -R '^reproducer_restart_twoway_state_matched_nopbl$' \
  --output-on-failure
```

A single np4-to-np1 comparison:

```bash
ctest --test-dir build \
  -R '^parallel_upwelling_twoway_state_matched_nopbl_np4_compare_to_np1$' \
  --output-on-failure
```

Anchor the regex with `^...$`. Without it, `-R parallel_upwelling_twoway_state`
matches all cadences and all four rank counts — 16 tests, not one — and `-R
..._np4` also matches the paired `_smoke` test.

Use `-N` to see what a selection would run before running it:

```bash
ctest --test-dir build -N -R np4_compare_to_np1
```

### Running one slice

By label:

```bash
ctest --test-dir build -L restart --output-on-failure     # the 8 restart tests
ctest --test-dir build -L parallel --output-on-failure    # the 64 MPI tests
```

Available labels: `compile`, `coupling`, `parallel`, `reproducer`, `restart`.

By axis, when you want one contract or one rank count across the matrix:

```bash
ctest --test-dir build -R '_flux_.*_np4_compare_to_np1' --output-on-failure
ctest --test-dir build -R 'reproducer_restart_twoway' --output-on-failure
```

### The two tests that are not Upwelling matrix entries

- `compile_erf_remora_driver` — build gate; other tests depend on it.
- `reproducer_coupled_sst_coverage` — the coverage-contract test. Drives
  `inputs_twoway` four times with different `driver.remap_min_coverage_frac`
  and domain overrides, and asserts on **covered-cell counts**, not on
  temperature. It is the only test in the suite that exercises the coverage
  threshold, because the Upwelling grids are otherwise fully covered.

```bash
ctest --test-dir build -R reproducer_coupled_sst_coverage --output-on-failure
```

### What the suite does not cover

The Upwelling decks are matched 1:1, all sea, with uniform SST and uniform
surface stress. That makes them structurally blind to the spherical remap, the
wet-mask exclusion, non-conformal coverage, and donor fill. A green suite means
"no regression against the idealized case" — it is not evidence that the
coverage or masking paths work. Those first execute on OCCAMS.

## References

All DOIs below were verified against CrossRef.

- Perlin, N., E. D. Skyllingstad, R. M. Samelson, and P. L. Barbour (2007).
  "Numerical Simulation of Air–Sea Coupling during Coastal Upwelling."
  *Journal of Physical Oceanography*. <https://doi.org/10.1175/JPO3104.1>
  — the Upwelling case.
- Warner, J. C., B. Armstrong, R. He, and J. B. Zambon (2010). "Development of a
  Coupled Ocean–Atmosphere–Wave–Sediment Transport (COAWST) Modeling System."
  *Ocean Modelling*. <https://doi.org/10.1016/j.ocemod.2010.07.010>
  — the coupled system this driver's remap and donor fill were ported from.
- Jones, P. W. (1999). "First- and Second-Order Conservative Remapping Schemes
  for Grids in Spherical Coordinates." *Monthly Weather Review*.
  <https://doi.org/10.1175/1520-0493(1999)127%3C2204:FASOCR%3E2.0.CO;2>
  — SCRIP; the conservative remapping scheme and the `fracarea`/`dstarea`
  normalization distinction.
- Craig, A., S. Valcke, and L. Coquart (2017). "Development and performance of a
  new version of the OASIS coupler, OASIS3-MCT_3.0." *Geoscientific Model
  Development*. <https://doi.org/10.5194/gmd-10-3297-2017>
- Danabasoglu, G., et al. (2020). "The Community Earth System Model Version 2
  (CESM2)." *Journal of Advances in Modeling Earth Systems*.
  <https://doi.org/10.1029/2019MS001916>

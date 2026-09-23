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
../../build_cpu_noahmp/erf-remora-driver inputs_oneway
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
../../build_cpu_noahmp/erf-remora-driver inputs_twoway
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
srun -N 1 -n 16 ../../build_cpu_noahmp/erf-remora-driver inputs_OCCAMS_oneway
```

---

### `OCCAMS/inputs_OCCAMS_state`

Purpose:

- `two_way` OCCAMS with the **state-passing** atmosphere-to-ocean contract.
  Exercises the ocean-to-atmosphere SST return path in the spherical frame.

Run:

```bash
cd Exec/OCCAMS
srun -N 1 -n 16 ../../build_cpu_noahmp/erf-remora-driver inputs_OCCAMS_state
```

---

### `OCCAMS/inputs_OCCAMS_flux`

Purpose:

- `two_way` OCCAMS with the **flux-passing** contract. Same SST return path as
  `inputs_OCCAMS_state`, with the surface-flux lanes on top.

Run:

```bash
cd Exec/OCCAMS
srun -N 1 -n 16 ../../build_cpu_noahmp/erf-remora-driver inputs_OCCAMS_flux
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
srun -N 1 -n 16 ../../build_cpu_noahmp/erf-remora-driver inputs_OCCAMS_smoke
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

Boundary cache: on its first run ERF converts `wrfbdy_d01` to its own levels and
writes the result to `Exec/OCCAMS/erfbdy/`, then every later run loads that
directory instead of the NetCDF file. Nothing records which ERF wrote it, and
the only consistency check is the variable count. A cache written by an older
ERF is loaded silently and can be wrong — one written before ERF `f0ac0b388`
carries a height-interpolation bug and blasts all four lateral edges within the
first step. After any ERF submodule update, move
`erfbdy` aside; the next run regenerates it in a few minutes. The log line that
tells you which happened is `Loading boundary data from erfbdy file` (cache) versus
`Initialized erfbdy file` (fresh conversion). `run_occams.sh` warns when the cache
predates the driver binary.

Cadence and stability: `max_step` counts ocean steps. With the deck's
`remora.fixed_dt = 300` and `erf.fixed_dt = 10`, one `max_step` is thirty ERF
steps and one ocean step, about ten minutes on sixteen ranks. Three settings that
differ from ERF's and REMORA's defaults are load-bearing and are explained in
`translation_map.tsv`:

- `erf.vert_implicit_fac = 1 1 1` in `inputs_ERF` makes the PBL vertical
  diffusion implicit in every RK stage, as WRF's YSU is. With the default the
  third stage is explicit and a two-cell vertical mode grows over the steepest
  terrain and kills ERF step 10.
- `erf.avg_grid_faces_to_nodes = true` in `inputs_ERF` keeps WRF's vertical
  levels. ERF's default rebuilds a grid from the 20 m surface layer that is about
  three times finer than WRF's between 500 m and 2 km, and resolved afternoon
  thermals over the Nevada desert then reach the vertical CFL limit at ERF step
  64. The boundary cache depends on the levels, so move `erfbdy` aside whenever
  this flag changes.
- `remora.ndtfast = 300` in each coupled deck keeps the barotropic step at the
  1 s the ROMS `ocean.in` used. Inheriting the base deck's 15 with a 300 s step
  makes it 20 s and the ocean's first step fails. Any deck that changes
  `remora.fixed_dt` must scale `remora.ndtfast` with it.

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
  occupancy, the resolved donor-fill mode, and per-direction coverage summaries.
- `driver.remap_v = 2` — additionally the per-REMORA-box staged lon/lat bounding
  box, one line per box. That box is the geographic clip the `coawst` donor search
  is restricted to, so this is the view to raise when a donor choice looks wrong.

Warnings and the coverage summary print regardless. So do the donor-fill mode
warnings: opting out of the default is never silent, whatever `remap_v` says.
Raise `remap_v` when a coverage number looks wrong and you need the geometry
behind it.

### What contract is running

```
Atmos->ocean contract mode: <name> (<N> active layers)
Atmos->ocean active lanes: <names>
```

The lane set is fixed for the whole run. Flux decks used to substitute the
state contract at `ocean_step = 0`, because the atmos->ocean pack happened before
ERF advanced and the stress tensor was still zero at that point. The pack now
follows the atmosphere substeps, so real stress exists at the first ocean step and
the substitution has been removed.

One consequence for a deck with no ERF radiation model: the old step-0 state pass
withheld SWrad/LWrad, so REMORA used its own `remora.swrad_type`/`lwrad_type`
forcing for that one step. Flux mode passes zeros deliberately, so those lanes are
now zero from the start rather than for every step but the first.

### Where each forcing lane's data comes from

Startup prints one row per atmosphere-to-ocean lane:

```
Atmos->ocean lane provenance (configured mode=state)
  lane    source  reason
  Uwind   ERF     packed by ERF
  Qair    REMORA  ERF does not fill it; REMORA default stands
  rain    REMORA  remora.rain_type or its value set in deck
  LWrad   ERF     packed by ERF; REMORA longwave set to external downward
```

Each lane resolves in order: one the deck configured explicitly stays with
REMORA; otherwise ERF supplies it if ERF's pack fills it; otherwise it is
withheld and REMORA's own default stands. Today ERF never fills `Qair` or
`rain`, fills `SWrad`/`LWrad` only with a radiation model, and `cloud` only with
cloud liquid or ice.

Withholding is how the driver declines a lane, not a failure: REMORA treats any
non-null lane as authoritative, so a lane of zeros would read as real forcing and
displace REMORA's own. Withheld lanes also appear during validation as
`<lane>: withheld (not supplied by ERF)`. The flux contract is exempt — REMORA
drops the whole apply if any flux lane is null.

`LWrad` gets one extra step: `remora.lwrad_type` defaults to `computed`, under
which REMORA's bulk fluxes ignore the field entirely, so when that lane resolves
to ERF the driver points REMORA at it. A deck setting `lwrad_type` owns the lane.

Override with `driver.atm2ocn_lane_source.<lane> = erf | remora | auto`, using
the printed lane name. Asking for `erf` on a lane ERF cannot fill aborts rather
than passing an unwritten lane.

### Atmosphere-to-ocean time averaging

The driver uses the atmosphere state at the end of each ocean interval by
default. Two independent, default-off controls replace selected lanes with the
trapezoidal mean over that same interval:

- `driver.time_average_atm2ocn_state = 1` averages `Uwind`, `Vwind`, `Patm`,
  `Qair`, `Tair`, and `cloud`.
- `driver.time_average_atm2ocn_flux = 1` averages `rain`, `SWrad`, `LWrad`,
  `tau_x`, `tau_y`, `SHflux`, `LHflux`, and `evap`.

Lane character, not `driver.atm2ocn_mode`, selects the control. In particular,
`rain`, `SWrad`, and `LWrad` follow the flux-lane control even though they also
appear in the state contract. The controls can be enabled separately or
together. Non-averaged lanes retain the interval-end sample, and averaged
exchange plotfiles are stamped at the interval midpoint. The accumulator is
reset and consumed within one ocean macro-step; it is not checkpoint state.

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
[driver] remap ERF corner mesh: nodes=... max corner-mean offset from cell centers=... deg
```

The corners are averaged from ERF's cell-centered `XLONG`/`XLAT`. The
`corner-mean offset` is how far each cell's four corners, averaged, land from
that cell's reported center. Small is expected (a curvature term). A sizable
fraction of a cell means the mesh is misplaced, e.g. a staggering mismatch
between ERF and the driver.

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

Reported once per weight family — `A2O scalar/rho`, `A2O tau_x/u-face`,
`A2O tau_y/v-face`, `A2O wind->rho` — because each is built against its own
destination stagger and each can be uncovered differently.

```
[driver] WARNING remap (<label>): <N> of <M> REMORA wet cells have no ERF source
overlap. ApplyConservativeRemap in Submodules/ERF/Source/Coupling/ERF_to_REMORA.cpp
branches on the -1 index sentinel and assigns a per-lane fallback: ...
```

Do not read that as "these cells fall back to a standard atmosphere". Two lanes
get a standard-atmosphere value. Every other lane gets zero, or nothing:

| Lane | What an uncovered wet cell receives |
|---|---|
| Patm | 101325 Pa — sea-level standard |
| Tair | 288.15 K — 15 C, standard |
| Uwind, Vwind | 0 m/s — calm, not a standard-atmosphere value |
| Cloud, SWrad, LWrad | literal 0.0 |
| all eight flux lanes | literal 0.0 |
| Qair, Rain | not remapped at all; they keep whatever `PrefillAtmosToOceanFallbacks` wrote from `remora.air_humidity` / `remora.rain`, both defaulting to 0.0 |

So "forced but not coupled" is true of Patm and Tair. For the wind and flux lanes
an uncovered cell is neither forced nor coupled. It is also an assignment, not
retention: the remap zeroes the destination first and then writes the fallback,
so a prefill does not survive on the lanes that are remapped.

Expect a small count along coastlines where the two land masks disagree, or where
the ocean grid reaches past the atmosphere. A large count means the atmosphere
does not cover the ocean nest.

Under `remap_v = 1`, the per-family summary:

```
[driver] remap coverage (<label>): cells=... uncovered=... no_stencil=...
                                   donor_filled=... frac min/mean/max=.../.../...
                                   max|frac-1|=...
```

- `uncovered` — destination cells whose weight sum is zero.
- `no_stencil` — destination cells whose stored source index is still the `-1`
  sentinel. The ERF-side fallback keys on this one, not on the weight sum. The two
  must agree and the driver asserts it, because a cell with an index but no weight
  would receive 0 Pa / 0 K rather than its lane fallback.
- `donor_filled` — cells whose stencil is a donor graft rather than a real
  overlap. Always `0` under the default `driver.remap_donor_fill = withhold` and
  in the cartesian frame; see the next section. The field is absent from the
  `O2A SST` line, which has no donor mechanism at all.
- `max|frac-1|` is the coverage error, not a conservation error. Before
  normalization the weight sum *is* the covered fraction, so a destination cell
  the source grid only half covers reports `0.5` here and is working correctly.
  What it does catch is a `frac_max` above 1, which is polygon-clipping round-off.

The all-clear, also under `remap_v = 1`:

```
[driver] remap coverage (<label>): all <M> wet cells covered by real source
overlap (donor_filled=0), min wet frac=<x>
```

It is suppressed whenever any cell was donor-filled, so a clean sheet is a clean
sheet.

### Donor fill: `driver.remap_donor_fill`

What a wet REMORA cell with no overlapping ERF source cell actually gets. The flag
only has an effect in the spherical frame — the donor search never runs in the
cartesian frame — and under `remap_v = 1` the driver says which mode resolved:

```
[driver] BuildRemapWeights: donor fill = withhold (inert in the cartesian frame)
```

**`withhold`, the default.** No stencil is written. The cell is counted uncovered,
the warning above fires, and `ApplyConservativeRemap` assigns the per-lane
fallback in the table above. Nothing is stubbed out: withholding *is* the
behaviour, and it is what the production couplers do with an unmapped point. ESMF
defaults `unmappedaction = ERROR`; CMEPS pre-fills and maps with
`zeroregion = REGION_SELECT`; YAC does not interpolate unmapped points
(`partial_coverage` defaults false); CDO pre-fills with a missing value; CESM
withholds by merge weight. SCRIP and OASIS are the only comparable codes that fill
at all, and both fill with zero — OASIS documents that as a deficiency.

**`coawst`.** A full-weight stencil is grafted pointing at the nearest non-land
ERF cell, reproducing COAWST's `check_weights` (`Lib/SCRIP_COAWST/scrip.f:343`).
Selecting it prints two lines unconditionally, not gated on `remap_v`:

```
[driver] WARNING remap: driver.remap_donor_fill=coawst grafts a full-weight
nearest-neighbour stencil into the fracarea-normalized conservative matrix ...
The resulting atmosphere->ocean operator is NOT conservative.
```

```
[driver] remap: coawst donors are restricted to source cells whose centroid lies
inside the destination box's own grown lon/lat bounding box, and to within
driver.remap_donor_max_diagonals=<r> donor diagonals ...
```

On a cartesian deck a third line follows, so a request that cannot take effect
does not look like one that did:

```
[driver] remap: driver.remap_donor_fill=coawst has no effect in the cartesian
frame; the donor search only runs when the frame is spherical.
```

The non-conservation claim is not a hedge. ESMF refuses to combine extrapolation
with conservative regridding at all, because "doing so would result in
non-conservative weights" — and that is precisely what this graft is: a
full-weight nearest-neighbour row spliced into a `fracarea`-normalized
conservative matrix. Bounding the donor, as below, makes *which* cell is grafted
more reproducible. It does not make the graft conservative, and no bound will.

Donor-filled cells are reported as their own category and are never counted as
covered:

```
[driver] WARNING remap (<label>): <D> of <M> REMORA wet cells have no ERF source
cell of their own and were donor-filled from the nearest non-land ERF cell
(driver.remap_donor_fill=coawst). ...
```

Their weight sum is nonzero, so they are not in `uncovered` or `no_stencil` —
those two keep exactly the meaning they have under `withhold`. They are in
`donor_filled`, they are excluded from `min wet frac` (a graft normalizes to a
flawless 1.0 and would flatter the statistic), and their presence suppresses the
all-clear line.

#### What bounds the donor search, and what that does and does not buy

Three restrictions apply, in this order. Only the first two are principled.

1. **A destination box with no geographic ERF overlap is not donor-filled at
   all.** For such a box the staged ERF region is not a real region: it is a
   one-cell placeholder parked at the ERF domain corner, because a `BoxArray`
   cannot hold an empty box. Without this test the search selects that corner, and
   every wet cell in the box takes the ERF domain corner's atmospheric state. It
   is rejected on a recorded flag rather than a distance threshold, because no
   threshold rejects it reliably — two nearly-touching grids could put the corner
   within any radius.
2. **Candidates are clipped to the destination box's own grown lon/lat bounding
   box.** This one is a proof, not a heuristic. Every ERF cell inside that box
   belongs to some ERF box whose own bounding box intersects it, and the staging
   plan unions every such ERF box into the staged region, so the surviving
   candidate set is identical under **every ERF decomposition** — whatever
   `NProcs()` and `amr.max_grid_size*` do to ERF.

   It removes dependence on **REMORA's** decomposition **not at all.** The
   bounding box is itself a per-REMORA-box quantity, REMORA's boxes are chopped by
   the same `NProcs()`, and a destination cell near the edge of one can have the
   same donor accepted under one REMORA box layout and rejected under another.
   `coawst` mode is bounded and empirically checkable; it is not proven
   rank-invariant. Only a global donor search would close that, and it is not
   implemented.
3. **An acceptance radius**, `driver.remap_donor_max_diagonals` (default `3.0`),
   measured against the donor cell's own control-volume diagonal in the same
   `cos(latitude)`-corrected metric the search uses.

   **This bound is a deliberate departure from precedent, not borrowed practice.**
   No production coupler bounds a nearest-source fill by distance. ESMF's
   `NEAREST_STOD` is explicitly unbounded — when no local candidate exists it
   queries its process directory with an infinite box, so it always finds the
   globally nearest source however far away — and creep fill bounds by iteration
   count rather than by distance. `3.0` is a judgement call and has not been
   calibrated against any production grid pair. If the number matters to your
   configuration, sweep it and watch `donor_filled` and the uncovered counts move;
   a plateau means it is not critical, a steep curve means it is.

   Set under `withhold` it does nothing, and the driver says so rather than
   letting a set value look effective:

   ```
   [driver] remap: driver.remap_donor_max_diagonals is inert under
   driver.remap_donor_fill=withhold ...
   ```

A rejected donor is not a third outcome: the cell falls through to the same `-1`
index sentinel and takes the same per-lane fallback it would under `withhold`.

Ties in the donor search are broken by smallest global source ID, which is what
ESMF and YAC do, rather than by scan order. Exact ties are common rather than
exotic on a structured grid — a destination cell centred between two source cells
gives bitwise-equal distances — and scan order starts at the staged region's
corner, which moves with the decomposition.

### Warnings that degrade silently

These do not stop the run and nothing downstream reports them:

```
[driver] remap: <N> wet destination cells had no valid source overlap and were
filled from their nearest non-land source cell (COAWST check_weights behaviour).
```

Only under `driver.remap_donor_fill = coawst`; the default never prints it. It
matches COAWST (Warner et al. 2010), the coupled system this driver's remap was
ported from, and those cells carry a neighbour's atmospheric state rather than
their own.

`<N>` here is a raw running total, summed over all four weight families and
incremented once per box visit, so a donor-filled face on a box seam can be
counted once per adjacent box and the number moves with the decomposition. Use
the per-family `donor_filled=` fields on the coverage lines when comparing runs:
those apply the same seam-ownership rule as every other coverage count. Watch the
trend either way — a count that *rises* between runs means the two land masks are
drifting apart.

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

`remora_boxes_without_overlap` counts REMORA boxes whose grown geographic
footprint intersected no ERF box at all. Those boxes hold a placeholder staged
region rather than a real one, and they are the boxes the `coawst` donor search
refuses to fill. If it is zero on your configuration, that failure mode never
arises there.

Under `remap_v = 2` each box's own footprint is printed, one line per box:

```
[driver] remap staging box <d> bbox=[(<lon_lo>, <lat_lo>), (<lon_hi>, <lat_hi>)]
```

That is the same box the `coawst` donor search clips its candidates to.

The staging path asserts three assumptions about solver layouts: that REMORA psi
boxes are nodal conversions of the rho mask boxes, that `lon_psi` and `lat_psi`
share box indexing and processor mapping, and that ERF corners are nodal
conversions of the O2A weight boxes. If a real layout violates one, the run
aborts on the assertion rather than mis-mapping quietly.

## Tests

The CTest suite drives the Upwelling decks with runtime overrides — 111 tests,
all from `Exec/Upwelling`. OCCAMS is not in the suite; its external data is not
in the repository.

### Run the suite by label, not in one `ctest` call

```bash
ctest --test-dir build -L compile  --output-on-failure --no-tests=error
ctest --test-dir build -L coupling --output-on-failure --no-tests=error
ctest --test-dir build -L restart  -j4 --output-on-failure --no-tests=error
ctest --test-dir build -L parallel -j4 --output-on-failure --no-tests=error
```

**`compile_erf_remora_driver` must not run alongside anything else.** It invokes
the build, and the driver relinks on every invocation, so the binary is replaced
underneath any test currently executing it. Those tests die on `SIGILL` (exit
132) with an `mpiexec ... could not access or execute` message that points at the
executable rather than at the race. A plain `ctest -jN` over the whole suite
schedules the compile test concurrently with the rest and produces a burst of
exactly these. `.github/workflows/ci.yml` avoids it by invoking `ctest` once per
label, which is why the commands above are split.

Within a label, concurrency is fine. All 108 matrix tests take a distinct
`ERF_REMORA_WORK_SUBDIR`, so no two can `file(REMOVE_RECURSE)` each other's
output. Measured: `-L restart` passes 12/12 at both `-j3` and `-j6` (222 s and
62 s), `-L parallel` passes 96/96 at `-j3`. `PROCESSORS` is set from the rank
count, so `ctest` will not oversubscribe the `np8` tests.

Two suites sharing one build directory is a different matter — they do collide,
because the work subdirs are keyed by test name, not by run. Use separate build
directories.

### Test name structure

Names are generated from the coupling axes, so you can select any slice with
`-R`:

```
reproducer_restart_<direction>_<contract>_<cadence>
parallel_upwelling_<direction>_<contract>_<cadence>_np<N>_<kind>
```

- `<direction>` — `oneway` or `twoway`
- `<contract>` — `state` or `flux`
- `<cadence>` — `matched` (one ocean step per atmosphere step) or `substep`
- `<N>` — `2`, `4`, `6`, or `8`
- `<kind>` — `smoke` (does it run) or `compare_to_np1` (does it match serial)

### Running one test

A single restart test:

```bash
ctest --test-dir build -R '^reproducer_restart_twoway_state_matched$' \
  --output-on-failure
```

A single np4-to-np1 comparison:

```bash
ctest --test-dir build \
  -R '^parallel_upwelling_twoway_state_matched_np4_compare_to_np1$' \
  --output-on-failure
```

Anchor the regex with `^...$`. Without it, `-R parallel_upwelling_twoway_state`
matches all cadences and all four rank counts — 24 tests, not one — and `-R
..._np4` also matches the paired `_smoke` test.

Use `-N` to see what a selection would run before running it:

```bash
ctest --test-dir build -N -R np4_compare_to_np1
```

### Running one slice

By label:

```bash
ctest --test-dir build -L restart --output-on-failure     # the 12 restart tests
ctest --test-dir build -L parallel --output-on-failure    # the 96 MPI tests
```

Available labels: `compile`, `coupling`, `parallel`, `reproducer`, `restart`.

By axis, when you want one contract or one rank count across the matrix:

```bash
ctest --test-dir build -R '_flux_.*_np4_compare_to_np1' --output-on-failure
ctest --test-dir build -R 'reproducer_restart_twoway' --output-on-failure
```

### The three tests that are not Upwelling matrix entries

- `compile_erf_remora_driver` — build gate; other tests depend on it.
- `reproducer_coupled_sst_coverage` — the coverage-contract test. Drives
  `inputs_twoway` four times with different `driver.remap_min_coverage_frac`
  and domain overrides, and asserts on **covered-cell counts**, not on
  temperature. It is the only test in the suite that exercises the coverage
  threshold, because the Upwelling grids are otherwise fully covered.

```bash
ctest --test-dir build -R reproducer_coupled_sst_coverage --output-on-failure
```

- `reproducer_lane_provenance` — the lane-provenance test. Drives
  `inputs_oneway` three times, once per contract and averaging combination, and
  asserts on the **exchange_info sidecar** beside each exchange plotfile: that
  every lane belongs to its contract exactly once in contract order, and that
  each lane's character and averaged columns follow the contract and the
  `driver.time_average_atm2ocn_*` flags rather than a default.

```bash
ctest --test-dir build -R reproducer_lane_provenance --output-on-failure
```

Neither of these carries the `restart` or `parallel` label, so neither is run by
`.github/workflows/ci.yml`, which selects only `compile`, `restart` and
`parallel`. Run them locally or add a `-L coupling` step.

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

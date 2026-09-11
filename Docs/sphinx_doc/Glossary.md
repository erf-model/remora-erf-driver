# Glossary

The driver's own vocabulary, in one place. This exists because the same
concepts currently carry more than one name across the repository: `lane`
appears 75 times in `Source/` but never in `README.md`, while `cache` appears
39 times in `Source/` and 3 times in `README.md` but nowhere in the newer
documents, and `slab` appears 27 times in `Source/` and in no document at all.
These are case-insensitive substring counts, reproducible with `grep -roi`;
matching on word boundaries instead gives 57, 21 and 17, because most uses are
inside identifiers.

**Preferred term is given first; deprecated synonyms are marked.** New prose
should use the preferred term.

## Exchange

**Lane** *(preferred; also seen as* layer*)*
: One exchanged variable's slot in an exchange cache — `Uwind`, `SST`,
  `SHflux` and so on. Lanes are addressed by the index enums in
  `Source/ERFRemoraCouplingContract.H`. A lane is *active* when the running
  contract selects it, and allocated-but-inactive otherwise.

**Cache**
: The driver-owned `Vector<MultiFab*>` holding all lanes for one direction —
  `m_atmos_to_ocean` and `m_ocean_to_atmos`. Allocated as a superset so that
  switching contract mode does not change the allocation.

**Slab** *(source-only; avoid in prose)*
: A single level's `MultiFab` within a lane. The term appears only in the
  source; prefer "lane" in documentation, or name the `MultiFab` directly.

**Contract**
: Which set of lanes is active for atmosphere-to-ocean exchange, selected by
  `driver.atm2ocn_mode`:
  - `state` — the atmosphere exports near-surface states and REMORA computes
    surface fluxes internally through its bulk-flux path.
  - `flux` — the atmosphere exports computed stress, heat, moisture and
    radiation lanes directly.

**A2O / O2A**
: Atmosphere-to-ocean and ocean-to-atmosphere. Used in parameter names
  (`driver.plot_file_exchange_a2o`) and diagnostic output.

**Pack / apply**
: The driver-owned API boundary. The source solver *packs* its state into a
  cache; the destination solver *applies* it. ERF and REMORA never call each
  other directly.

## Remapping

**Frame**
: The coordinate system the remap weights are built in, set by
  `driver.remap_frame` — cartesian or spherical. Several controls are inert in
  the cartesian frame.

**Coverage**
: The fraction of a destination cell's area reached by valid source cells.
  Thresholded on area rather than on any overlap; `driver.remap_min_coverage_frac`
  sets the threshold.

**Withhold** *(the default donor-fill mode)*
: When a wet destination cell has no valid source overlap, write no stencil.
  The cell is counted uncovered and receives its per-lane fallback. This is
  what the production couplers do with an unmapped point.

**Donor fill**
: Grafting a stencil for an otherwise-uncovered cell, pointing at a nearby
  source cell. Selected by `driver.remap_donor_fill`; `coawst` mode reproduces
  COAWST's `check_weights` behaviour. Donor-filled cells are counted separately
  from genuinely covered ones.

**`fracarea` / `dstarea`**
: The two conservative-remap normalizations distinguished by SCRIP
  (Jones, 1999). `fracarea` divides by the covered fraction; `dstarea` divides
  by the full destination area.

**Staged box**
: A REMORA-side region — psi corners and wet mask — copied onto the ERF box
  layout so weight construction is rank-local.

## Wet mask

**Wet / dry**
: Ocean versus land cells in REMORA's mask. Coverage counts only wet
  destination cells; a missing wet mask is reported rather than treated as an
  all-wet sheet.

# Symbol Isolation Policy (`prefixed`)

## Why this exists

When ERF and REMORA are linked into one executable, both projects expose
collision-prone global symbols:

- `amrex_probinit` (weak problem initializer hook)
- `Problem`, `ProblemBase`, `SolverChoice`

In static-library superbuilds, weak `amrex_probinit` can also fail archive
extraction and resolve as a null call target.

## Policy

This repository enforces a fixed full-prefix policy in superbuild mode.

This superbuild uses full explicit symbol prefixing on both sides:

- ERF: `ERFProblem`, `ERFProblemBase`, `ERFSolverChoice`, `erf_probinit`
- REMORA: `REMORAProblem`, `REMORAProblemBase`, `REMORASolverChoice`,
  `remora_probinit`

## Weak probinit archive mitigation

In ERF library-only mode, `amrex_probinit` is weak and can be skipped during
archive extraction. To force extraction of `Exec/ERF_Prob.cpp`, ERF now adds a
small link anchor symbol referenced from `ERF.cpp` in superbuild mode.

## Build examples

```bash
cmake -S . -B build
cmake --build build -j4
```

## Upstream PR slicing

Recommended split for reviewable PRs:

1. Build/plumbing: top-level fixed prefix policy + superbuild CMake wiring.
2. ERF: probinit link-anchor mitigation + full prefixed compile defs.
3. REMORA: policy-aware compile definitions for combined executable safety.
4. Driver/docs: adapter-boundary notes and policy documentation.

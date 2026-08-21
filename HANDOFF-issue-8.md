# HANDOFF issue 8

## Changed

- Source/main.cpp:77-87 changes `RequireCheckpointSelfConsistency` to accept a time offset and reconstruct `offset + step * dt`, so REMORA headers validate on the shifted ocean clock while ERF keeps the zero-offset default.
- Source/main.cpp:91-111 passes `driver.ocean_time_offset` into `BuildRestartState`, subtracts it from REMORA checkpoint time, and uses the elapsed REMORA time for cross-solver time comparison and ocean tick reconstruction.
- Source/main.cpp:214-215 parses `driver.ocean_time_offset` from `pp_driver` with the same key and default behavior used by the coupled container.
- Source/main.cpp:259-264 rejects `amr.restart` before `InitializeSolvers()` can let ERF consume it without coupled restart validation.
- Source/main.cpp:282-283 passes the parsed offset into `BuildRestartState` when both explicit restart paths are present.

## Deliberately not done

- Did not touch `Submodules/`, `CMakeLists.txt`, `Exec/`, or `cmake/`.
- Did not alias `amr.restart` into `erf.restart`; the driver now aborts and requires `erf.restart` plus `remora.restart`.
- Did not move or broaden the existing AMR parameter scrub; this patch only blocks `amr.restart` before solver initialization.
- Did not build or run tests.

## Needs human decision

- Maintainers may choose whether `amr.restart` should remain rejected or become an explicit alias to `erf.restart` in a later change.
- No separate problem outside the issue's declared locations was fixed here.

## UNVERIFIED

- All runtime behavior remains unverified.
- No compile, CMake configure, MPI run, restart run, checkpoint read, or `amr.restart` deck run was executed.
- The offset restart path and `amr.restart` abort path still need the separate verification pass.

# Validated fracture configuration on geomech

Prepared on local branch `feat/geomech-validated-fracture-configuration`, based
on completion-identity PR #40 at `bbed49d8fb2126bf8e1b3401d8d15a298bc979cc`.
This is uncommitted review work. It does not include the old experimental grid
history, mechanics-constraint PR #41 or physics-harness PR #39.

## Purpose and behavior

Validate fracture settings immediately after loading a file, built-in alias or
C++ defaults, before the coupling/solver code reads them. Reject unknown,
misplaced and duplicate keys, invalid scalar conversions, nonfinite numbers,
selected invalid ranges, unsupported modes and combinations. Diagnostics name
the source and full parameter path. Preserve the underlying JSON/file error.

The registry has **267 entries**, reconciled with the clean parameter readers.
Outer `solver.*` settings remain distinct from `fractureparam.solver.*` settings.
Numeric and boolean strings remain compatible. Validation does not add defaults,
change input or convert unit-suffixed values. Scalar `__*` documentation fields
are allowed. FlexibleSolver's `flow_solver` subtree is delegated to that solver.

Compared with the preserved old work, the registry adds
`fractureparam.vem_face_quadrature` (`kite` or `triangulated`) and removes the
unread `fractureparam.solver.damping_p` and `damping_w` options. Newly enabled or
moved linear-solver reads were already covered. Separate sector-grid settings
remain outside the fracture-file schema. Ten compatibility fields are explicitly
labelled as having no active reader, including `hasfractures`, `control.WI`,
`full_intersections` and the legacy method fields; they are retained for built-in
and existing profile compatibility, not activated by this change.

## Input and effective seed reports

Rank zero logs validated input, source provenance, principal resolved modes and
unit annotations. Each locally created WSEED seed logs its exact initialization
tree and separate parent/completion/grid identities from PR #40. These reports
are initialization records, not a global MPI catalog or an exhaustive expansion
of solver-local defaults.

Each seed receives its own settings copy. WSEED vertical extent overrides
`config.axis_scale`; WSEED WIDTH overrides hydraulic `config.min_width`. The
numerical `solver.min_width` remains separate. Reject nonfinite/nonpositive
vertical extent and nonfinite/negative WIDTH. Valid values and equations are not
retuned. JSON pressure values are generally Pa, lengths m and rates m^3/s;
existing bar/day/K exceptions are annotated in the report and reference guide.

Five unused outer `solver.ctf_change_threshold` entries and one unused inner
`well_local_solve` entry are removed from existing profiles. They are removed,
not moved into active namespaces, preserving their prior numerical behavior.

## Focused validation on 23 October 2026

- Production-source configuration suite compiled and linked, then passed without
  simulator execution or MPI initialization.
- All **13 root fracture profiles**, **14 curated profiles**, C++ defaults used
  by the simulator and both aliases pass.
- Negative cases cover misspellings/namespaces, types/shapes/duplicates,
  nonfinite/range errors, unsupported modes/combinations, removed damping keys
  and misplaced sector settings. Both clean VEM quadrature modes are accepted
  and reported faithfully.
- Report parsing, input immutability, independent WSEED overrides, SI values and
  parent/completion/grid identity checks pass.
- The legacy `makeDefaultFractureParam(true)` factory emits unsupported
  `perf_rate`; the suite checks its rejection. The simulator uses the false
  variant. This change does not introduce or silently enable perf_rate.
- The instantiated `flow_energy_geomech` translation unit passes a syntax check
  with the new constructor integration enabled.

Compilation was sequential with an 8 GiB process ceiling. The focused test uses
actual FractureModel defaults, built-in profiles and PropertyTree sources;
function-section linking avoids building/linking the full simulator. Evidence
and the literal-reader audit are under the workspace's
`geomech/validation/configuration/`. The original fifteen uncommitted files are
preserved unchanged in the old checkout and preservation package.

After a configured full build, the normal test is:

```bash
cmake --build build --target test_fracture_configuration -j1
ctest --test-dir build/opm-flowgeomechanics \
  -R '^Unit.FractureConfiguration$' --output-on-failure
```

## Dependencies and remaining merge gates

**Depends on #40** for trustworthy grid-qualified completion reports. Prefer
integrating #40 first, then apply the configuration-only patch onto updated
hnil/geomech and create a separate PR. Opening this branch against geomech before
#40 merges would also show #40's commit; that is not a configuration-only diff.
Build validation also needs the common HYPRE/build prerequisites recorded in
workspace GEOMECH_PR_MIGRATION_REPORT.md. #39 is useful for physics comparisons
and #41 for mechanics correctness, but neither is a code dependency here.

Full clean-line simulator linking, nine regressions, malformed-configuration
simulations on one/four ranks, fresh SIMPLE/WELLREF physics/report comparisons
and TPSA runtime coverage remain pending. Historical old-line successes are not
claimed for this port. Keep any eventual PR as a draft until integration checks
pass. Use the workspace safe-build wrapper and one job; do not repeat the old
six-job build commands.

## Limits and maintenance

Update registry entries and tests alongside new readers. Type conversions follow
PropertyTree, which erases JSON scalar types and cannot distinguish all empty
container/scalar cases. Bounds and required-field checks are selective. The
FlexibleSolver subtree, reserved documentation, compatibility-only fields and
separate sector file are explicit scope boundaries. No public JSON Schema,
unit conversion, physical-validity proof, schedule-update report or new
filter-cake/convergence/restart behavior is introduced.

# Agent notes for gaato/stats

## Development

- Use `moon` on PATH and the toolchain pinned by `moonbit-version`.
- Packages are layered: `num` → `special` → `data` → `dist`, `multi` →
  `describe` → `fit` → `infer` → `regress`, `survival`, `timeseries`,
  `multivariate`, `bayes`. A package imports only layers before it in this
  list, so the numerical base never depends on the statistics built on it.
- These files are generated; edit the generator and rerun it with
  `uv run --no-project --with scipy --with mpmath python -I scripts/<name>.py`,
  because the next run overwrites hand edits:
  - `src/special/constants.mbt` from `gen_special_constants.py`
  - `src/dist/catalog.mbt` and `src/dist/catalog_test.mbt` from `gen_catalog.py`
  - `src/special/special_test.mbt` from `gen_special_tests.py`
  - `src/infer/infer_oracle_test.mbt` from `gen_infer_tests.py`
- A catalogue family's method written by hand (in `closed_forms.mbt` or
  `extras.mbt`) must also be listed in `HAND_METHODS` in `gen_catalog.py`;
  otherwise the generator emits a stub for it too and the duplicate does not
  compile. A family gets oracle tests by having its distribution function in
  `mp_cdfs` there.
- Check new numerics against mpmath at high precision, not against SciPy:
  SciPy was off by up to 2e-10 in tails where the closed forms here were
  exact. Before loosening a tolerance, compute the case's conditioning with
  mpmath; tolerances are about twice the measured worst error
  (`gen_special_tests.py --report` prints it), and an error far above the
  conditioning is a bug.
- Never compute an upper tail as `1 - cdf`; each tail has its own formula.
  R's nmath and GSL are GPL: compare against them, never copy or follow their
  code.
- A path not implemented yet raises `Unsupported` (or is a typed `...`
  body); it does not approximate silently. When a body lands, update its row
  in `DESIGN.md`, which maps every concept to its type and status. Its
  "Language notes" list the MoonBit behaviour that shaped the API (trait
  methods through objects, `pub` versus `pub(all)`, reserved words).
- Do not publish to mooncakes.io or push without the user's request.

## Checks

Run `moon fmt`, `moon check --target all`, and `moon test` on each of
`wasm-gc`, `js`, `native` and `wasm`: the backends use different
elementary-function libraries, so a tolerance that holds on one can fail on
another. Then run `moon info` and review the `pkg.generated.mbti` changes;
do not edit those files by hand.

## Documentation

- `README.md` is for visitors: what the module is, its status, the design
  in brief. Design rationale, accuracy tables and the language notes live in
  `DESIGN.md`.
- Before the first publish, follow the other gaato MoonBit repositories:
  move the README to `src/overview.mbt.md` with `README.md` as a symbolic
  link to it (its `mbt check` examples then compile as tests), add a
  `.moonignore` for repository-only files, and copy
  `scripts/docs_audit.mbtx` from gaato/discord.mbt rather than writing a new
  one.

## API changes

- When a published function's labelled parameter is renamed or removed, an
  optional one becomes required, or a positional one becomes labelled, keep
  existing calls compiling for one minor release with `#label_migration`
  (`alias=old`, `fill=false`, `fill=true` or `allow_positional=true`) and a
  `msg` naming the new form, and remove the attribute in the next minor
  release. Callers then get a warning instead of a build failure; the
  CHANGELOG still lists the change with its migration note. A renamed or
  removed function keeps a `#deprecated` alias for the same period.

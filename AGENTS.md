# Agent notes for gaato/stats

## Development

- Use `moon` on PATH and the toolchain pinned by `moonbit-version`.
- Keep imports pointing down the layer order in `DESIGN.md`, so that the
  numerical base never depends on the statistics built on it.
- A file whose first line names a generator in `scripts/` is generated:
  change the generator and rerun it with
  `uv run --no-project --with scipy --with mpmath python -I scripts/<name>.py`,
  because the next run overwrites hand edits. Code for a generated
  distribution family goes in another file of `src/dist`; the catalogue
  generator leaves out any method it finds implemented there.
- Test new numerics against mpmath at high precision. SciPy checks
  conventions, not accuracy, because its tails are not exact. Set a
  tolerance at about twice the worst measured error; before loosening one,
  compute the case's conditioning with mpmath, since an error far above it
  is a bug.
- Never compute an upper tail as `1 - cdf`; each tail has its own formula.
  R's nmath and GSL are GPL: compare against them, never copy or follow
  their code.
- A path not implemented yet raises `Unsupported` or is a typed `...` body,
  never a silent approximation. When one lands, update its row in
  `DESIGN.md`.

## Checks

Run `moon fmt`, `moon check --target all`, and `moon test` on each of
`wasm-gc`, `js`, `native` and `wasm`: each backend has its own
elementary-function library, so a tolerance can hold on one and fail on
another. Then run `moon info` and review the `pkg.generated.mbti` changes.

## Documentation

- `README.md` is for visitors deciding whether and how to use the module.
  Design rationale, accuracy tables and language notes go in `DESIGN.md`.

## API changes

- When a published function's labelled parameter is renamed or removed, an
  optional one becomes required, or a positional one becomes labelled, keep
  existing calls compiling for one minor release with `#label_migration`
  (`alias=old`, `fill=false`, `fill=true` or `allow_positional=true`) and a
  `msg` naming the new form, and remove the attribute in the next minor
  release. Callers then get a warning instead of a build failure; the
  CHANGELOG still lists the change with its migration note. A renamed or
  removed function keeps a `#deprecated` alias for the same period.

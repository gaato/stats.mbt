# Agent notes for gaato/stats

- Test new numerics against mpmath at high precision; SciPy checks
  conventions, not accuracy, because its tails are not exact. Before
  loosening a tolerance, compute the case's conditioning with mpmath: an
  error far above it is a bug.
- Never compute an upper tail as `1 - cdf`; each tail has its own formula.
- R's nmath and GSL are GPL: compare against them, never copy or follow
  their code.
- A path not implemented yet raises `Unsupported`; it never approximates
  silently.
- Run `scripts/gates.sh` before handing work back.

## API changes

- When a published function's labelled parameter is renamed or removed, an
  optional one becomes required, or a positional one becomes labelled, keep
  existing calls compiling for one minor release with `#label_migration`
  (`alias=old`, `fill=false`, `fill=true` or `allow_positional=true`) and a
  `msg` naming the new form, and remove the attribute in the next minor
  release. Callers then get a warning instead of a build failure; the
  CHANGELOG still lists the change with its migration note. A renamed or
  removed function keeps a `#deprecated` alias for the same period.

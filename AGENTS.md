# Agent notes for gaato/stats

- Test new numerics against mpmath at high precision; SciPy checks
  conventions, not accuracy, because its tails are not exact. Before
  loosening a tolerance, compute the case's conditioning with mpmath: an
  error far above it is a bug.
- Never compute an upper tail as `1 - cdf`; each tail has its own formula.
- R's nmath and GSL are GPL: compare against them, never copy or follow
  their code.
- A path not implemented yet raises an error that says so
  (`NotImplemented` in infer) or is a typed `...` body; it never
  approximates silently. `Unsupported` is for requests that have no
  answer by the method asked for, so callers can tell the two apart.
- Catalogue formulas are private `T::std_*` (continuous, standard form)
  and `T::lattice_*` (discrete, integer-valued Doubles) methods, which
  `scripts/gen_catalog.py` finds by name and wraps; its docstring lists
  them. A hand-written `pub impl` of a trait the generator also emits
  would be a duplicate that does not compile, and wrapping in one place
  is what keeps location, scale and every tail forwarded.
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
- The README is `src/overview.mbt.md`, with `README.md` linked to it, and
  every `pub` declaration has a `///` comment: `scripts/docs_audit.mbtx`
  checks these rules, which every gaato MoonBit repository shares so each
  mooncakes.io page reads the same. The script is a byte-for-byte copy of
  gaato/discord.mbt's; put this repository's settings in
  `scripts/docs_audit.json` and `scripts/docs_audit.allow`, never in the
  script. Generated files say `DO NOT EDIT` on their first line, and the
  audit leaves their docs to the generator.

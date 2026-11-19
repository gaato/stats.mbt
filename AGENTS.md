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

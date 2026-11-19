# gaato/stats

Statistics foundations for MoonBit: probability distributions, special
functions, descriptive statistics and hypothesis tests.

**Status: design stage.** The types are in place and type-check on every
backend; most bodies are still unimplemented, and nothing is published yet.

## Design

- **Distributions are the centre.** Each family is an immutable value; what a
  distribution can compute is expressed by capability traits (`Density`,
  `Cdf`, `Quantile`, `Moments`, sampling), with separate traits for
  continuous and discrete support.
- **Both tails are first class.** `sf`, `logcdf`, `logsf` and `isf` exist
  for every distribution, and families override them with formulas that do
  not lose digits in the tails.
- **Conventions are explicit.** Scale versus rate, the sample quantile
  definition and the variance divisor are chosen by the caller, never by a
  default.
- **Summaries compose.** Descriptive statistics are accumulators with an
  identity and an associative `merge`, so data can be summarised in pieces,
  in parallel or as a stream.
- **Assumptions are visible.** A test that relies on an assumption asks for
  a certificate value, which the caller obtains explicitly.
- **No global state.** Sampling takes the random source as an argument.

## Layout

| Package | Contents |
|---|---|
| `special` | Incomplete beta and gamma functions, normal tail functions |
| `dist` | Capability traits, families, empirical and finite discrete distributions |
| `describe` | Mergeable accumulators |
| `infer` | Test results, intervals, estimates |

## License

Apache-2.0

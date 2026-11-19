# Design map

Which concepts have types, which have working code, and how the ones not
yet typed are expected to fit. The point of typing concepts before
implementing them is to find out early whether the capability traits hold.

Status: **works** = implemented and tested; **typed** = signatures in
place, body is `...`; **planned** = shape described here, no code.

## Univariate distributions (`dist`)

Every distribution on the real line has a distribution function, whatever
its Lebesgue decomposition, so `Cdf`, `Quantile` and `Sample` on `Double`
are common; `Density` (absolutely continuous part) and `Mass` (atoms on
the integers) are extra capabilities.

| Concept | Type | Status |
|---|---|---|
| Distribution function, both tails, logs | `Cdf` | typed (defaults work) |
| Generalised inverse, upper-tail inverse | `Quantile` | typed |
| Sampling with an explicit random source | `Sample`, `RandomSource` | typed |
| Density / mass | `Density`, `Mass` | typed |
| Integer distribution function and inverse | `IntCdf`, `IntQuantile`, `SampleInt` | typed |
| Mean, variance; skewness, kurtosis; entropy | `Moments`, `Shape`, `Entropy` | works for some families |
| Characteristic function | `Characteristic` | works (normal, binomial) |
| Moment and cumulant generating functions | `Mgf` | works (normal, gamma) |
| Fitting to a sample | `Fit`, `FitError` | typed |
| Composite bundle for trait objects | `Univariate` | works |

Families: `Normal`, `StudentT`, `Gamma`, `Binomial`.
Built from data or tables: `Empirical` (quantile definition fixed at
construction), `FiniteDiscrete` (both cumulative sums kept).

## Composition (`dist`)

| Concept | Type | Status |
|---|---|---|
| Atom | `PointMass` | works |
| Mixture of any components | `Mixture` (trait objects) | works except quantile, sampling |
| Location-scale transform | `Affine[D]` | works |
| Truncation | `Truncated[D]` | works for `Cdf` |
| Censoring (atoms at detection limits) | `Censored[D]` | works for `Cdf` |
| Singular distribution | `Cantor` | moments work; rest typed |

## Multivariate (`multi`)

| Concept | Type | Status |
|---|---|---|
| Dimension, density, mass, sampling into buffers | `Dimension`, `MvDensity`, `MvMass`, `MvSample` | typed |
| Mean vector, covariance (row-major d * d) | `MvMoments` | typed |
| Marginal as a univariate trait object | `Marginal` | typed |
| Families | `MvNormal` (stored as Cholesky factor), `Dirichlet`, `Multinomial` | typed |
| Dependence separate from marginals | `Copula`, `GaussianCopula`, `Joint[C]` | sampling works once the copula does |

Needs dense linear algebra (Cholesky, triangular solves) before the
bodies can be written; the covariance stays a flat array until a matrix
type exists.

## Descriptive statistics (`describe`)

| Concept | Type | Status |
|---|---|---|
| Mergeable summaries (monoids) | `Accumulator` | works |
| Count, mean, variance with explicit divisor | `MeanVar` | works |

## Inference (`infer`)

| Concept | Type | Status |
|---|---|---|
| Alternative hypothesis | `Alternative` | works |
| Test result with statistic, p-value, log p-value | `TestResult` | typed |
| Interval with its level and procedure | `Interval` | typed |
| Point estimate with optional uncertainty | `Estimate` | typed |
| Assumptions as certificates | `Assumed[A]`, `Normality`, `EqualVariance` | works |
| Welch t test | `welch_t_test` | typed |

## Special functions (`special`)

Incomplete beta and gamma with their complements and inverses, the normal
tail functions, log-gamma: typed. The incomplete beta is the first body to
write; it unlocks the t, F, beta and binomial distribution functions.

## Planned, not yet typed

| Concept | Expected shape | Depends on |
|---|---|---|
| Directional (von Mises) | `Density`, `Cdf` on [-π, π) with explicit wrap semantics | — |
| Extreme value (GEV, Gumbel, Weibull) | ordinary families | — |
| Kernel density estimate | a `Density` + `Cdf` built from data, bandwidth chosen at construction | — |
| Distribution from a characteristic function | invert φ numerically with an FFT; result as `FiniteDiscrete` or a gridded `Density` | an FFT module |
| Exact rank-test null distributions | polynomial products / FFT convolution into `FiniteDiscrete` | an FFT module |
| Matrix-variate (Wishart, inverse Wishart) | a matrix-valued `MvSample` | a matrix type |
| Bayesian conjugate updates | functions `posterior(prior, data) -> prior`, not a trait: traits have one `Self` | — |
| Likelihood for models | `logpdf` summed over data, then the optimisation layer | an optimisation module |
| Regression and multivariate analysis | `Fit`-like results with covariance and diagnostics | dense linear algebra (QR, SVD) |
| Stochastic processes (Poisson, random walks, Markov chains) | out of scope for the base library | — |

## Language notes (MoonBit 0.10.21)

- Methods of a supertrait are not reachable with dot syntax through a
  trait object; call them as `Cdf::cdf(d, x)`.
- Trait methods implemented for a concrete type are not reachable with dot
  syntax either, even in the same package, unless exposed with
  `pub extend T with Trait::{...}`; otherwise call them qualified.
- Traits have one `Self` and no associated types, so continuous and
  discrete support use separate traits, and a two-type relation (prior and
  likelihood) is a function rather than a trait.
- An empty struct literal is written `T::{ }`; a bare `{ }` is a map.

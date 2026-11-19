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
| Upper tail including an atom, P(X >= x), for p-values from discrete nulls | `Cdf::sf_at_least` | works for atoms, mixtures, binomial |
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

## Families and parameters (`dist`)

A family is a distribution with its parameters left open. Traits cannot
name the type a family produces (no associated types), so a family is a
record of functions, generic in the distribution type; fixed
hyperparameters and alternative parameterisations are separate family
values.

| Concept | Type | Status |
|---|---|---|
| Parameter constraints and the unconstrained scale | `Constraint`, `ParamSpec` | works |
| Family: parameter vector to distribution and back | `Family[D]` | works (normal, t, gamma by rate and by scale, binomial with fixed n) |
| Exponential family: η, T, A, h, mean-matching | `ExpFamily[D]` | works (normal); binomial typed |
| Conjugate prior and likelihood | `Conjugate[P, D]` | works (normal-normal, known variance) |
| Log-likelihood with compensated summation | `loglik`, `loglik_int` | works |
| Score (gradient in the parameters) | `Score` | typed |

## Fitting (`fit`)

| Concept | Type | Status |
|---|---|---|
| Fit result with log-likelihood, k, n, convergence, covariance | `MleResult[D]` | works |
| AIC, AICc, BIC | `MleResult::aic`, `aicc`, `bic` | works |
| Exponential-family MLE by mean matching | `fit_exp_family` | works (densities) |
| General MLE on the unconstrained scale | `fit_mle` | typed; needs an optimiser |
| Model comparison | `Candidate`, `akaike_weights` | works |

Gap: `fit_exp_family` asks for a `Density`, so discrete exponential
families (binomial, Poisson) need an integer-data variant.

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
| Question: estimand and data shape, checked together | `Estimand`, `Sample`, `Question`, `QuestionError` | works |
| Method and null hypothesis | `Method`, `Hypothesis` | typed |
| Estimate for a question | `estimate` | typed |
| Observed statistic with its null distribution | `NullModel`, `null_model` | p-values work |

After R's infer package (specify, hypothesize, generate, calculate): a
null distribution is an ordinary distribution, parametric for a
theory-based test and `Empirical` for a randomisation test, so one
`NullModel::p_value` serves both. Each tail includes the observed value.

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
- A `pub struct` cannot be constructed outside its package; `pub(all)`
  can. Results the library produces (`TestResult`, `MleResult`) stay
  `pub` so that callers cannot fabricate them; plain data callers supply
  (`Candidate`, `NullModel`) is `pub(all)`.
- `method` and `var` are reserved words.
- A function stored in a record field is called as `(f.make)(x)`.

# Design map

Every concept the library is meant to cover, the type that carries it,
and its status. Concepts are typed before they are implemented, so that
the capability traits and data shapes are tested against the whole scope
before a hundred distributions or tests are written against them.

Status: **works** = implemented and tested against an oracle;
**typed** = signatures in place, body `...`; **planned** = described
here, no code.

## Where the scope comes from

Three surveys, kept in `~/works/numeric-proto/stats-survey/` (not part of
this repository):

- **Demand on mooncakes.** 48 statistics modules, each reimplementing what
  it needs (only one depends on another statistics module). Most common:
  mean (31), quantile (31), variance (27), median (23), histogram (22),
  confidence intervals (18), weighted statistics (18), Pearson
  correlation, covariance, robust location and scale (16 each), OLS (14),
  bootstrap (14), normal distribution (13). Formal tests are rare. Bugs
  found include a wrong incomplete beta (t-test p-value 0.372 instead of
  0.628), a chi-square CDF that is not one, a log-rank chi-square tail of
  exp(-x/2) for every df, 1 - cdf upper tails in seven modules, z where t
  is needed, and no module letting the caller choose a quantile definition.
- **Scope of mature libraries.** R `stats`, `scipy.stats`, statsmodels,
  HypothesisTests.jl and StatsBase.jl: 18 categories, 30 input shapes,
  the fields of test results, the resampling API, and the conventions the
  libraries resolve differently.
- **What practitioners are expected to know.** The scope of the Japan
  Statistical Society certificate (grades 2, pre-1 and 1), and NIST StRD
  for accuracy (univariate, ANOVA, linear and nonlinear regression, MCMC,
  with certified values to 15 digits).

## Layers

```
num        log-space, compensated sums, solvers with status, random sources
special    incomplete beta/gamma, normal tails, gamma family, ...
data       observations, shapes, weights, tables, conventions on data
dist       distributions, families and parameters, composition
multi      multivariate distributions, copulas
describe   accumulators (monoids), quantiles, ranks, robust summaries
fit        maximum likelihood, information criteria
infer      questions, null models, tests, multiple testing, power, resampling
regress    linear models, GLM, penalised and nonlinear least squares
survival   Kaplan–Meier, log-rank, Cox, parametric survival
timeseries ACF/PACF, spectra, ARIMA, state space
multivariate PCA, factor analysis, discriminants, clustering, ROC
bayes      MCMC, diagnostics, EM
```

Each layer imports only layers above it in this list.

## num

| Concept | Type | Status |
|---|---|---|
| log(e^a ± e^b), log1p(e^x), log(1 − e^x), log Σ e^x | `log_add_exp`, `log_sub_exp`, `log1p_exp`, `log1m_exp`, `log_sum_exp` | works |
| Compensated summation (Neumaier), mergeable | `CompensatedSum` | works |
| Status of an iteration | `Status`, `Solution`, `Tolerance` | typed |
| Bracketed root finding (Brent; Newton inside a bracket) | `find_root`, `find_root_newton` | typed |
| Continued fractions (modified Lentz), series, quadrature | `continued_fraction`, `sum_series`, `integrate` | typed |
| Domain errors for checked variants | `DomainError` | typed |
| Random source, generator, streams | `RandomSource`, `Xoshiro256`, `Jump` | generator works (reference outputs); `jump` typed |
| Uniforms | `uniform`, `uniform_open`, `uniform_below` | first two work |
| Standard variates with frozen, caller-chosen methods | `NormalMethod`, `ExponentialMethod`, `GammaMethod`, `DiscreteMethod`, `shuffle` | typed |

Reproducibility: for a given generator, seed and method the stream never
changes; a better algorithm is a new enum case. NumPy's NEP 19 chose
differently (bit streams stable, variate algorithms may change); see
open decisions.

## special

Incomplete beta and gamma with complements and inverses, normal tails and
inverse, log-gamma, erf family, digamma and trigamma, log-beta, log
binomial coefficient, Owen's T, the Kolmogorov distribution, the
studentized range: typed. Error policy after Boost.Math: plain functions
return NaN outside their domain; checked variants (planned) raise
`DomainError`. The incomplete beta comes first: it unlocks t, F, beta and
binomial.

## data

| Concept | Type | Status |
|---|---|---|
| Observation with censoring | `Observation` | works |
| Missing values | `MissingPolicy`, `clean` | works |
| Weights and their meaning | `WeightKind`, `Weighted` | works |
| Contingency table | `Table` | works |
| Quantile definition, rank ties | `QuantileMethod`, `RankTies` | works |
| (x, y) pairs, long format with labels, blocks | `Bivariate`, `Grouped`, `Blocks` | `Bivariate` works |
| Published summaries, proportions, category counts | `Summary`, `Proportion`, `CategoryCounts` | typed |
| Survival records (time, event, entry) | `SurvivalRecord` | works |
| Data matrix, time series, angles | `DataMatrix`, `TimeSeries`, `Circular` | containers work |

## dist

Every distribution on the real line has a distribution function, so `Cdf`,
`Quantile` and `Sample` are common; `Density` and `Mass` are extra.

| Concept | Type | Status |
|---|---|---|
| Distribution function, both tails, logs, P(X ≥ x) | `Cdf` (incl. `sf_at_least`) | works |
| Quantile, upper-tail inverse | `Quantile` | typed |
| Density, mass, integer convenience | `Density`, `Mass`, `IntCdf`, `IntQuantile`, `SampleInt` | typed |
| Moments, shape, entropy | `Moments`, `Shape`, `Entropy` | works for some |
| Characteristic and moment generating functions | `Characteristic`, `Mgf` | works (normal, binomial; normal, gamma) |
| Fitting, score | `Fit`, `Score` | typed |
| Composite bundle (trait objects) | `Univariate` | works |
| Families | `Normal`, `StudentT`, `Gamma`, `Binomial` | partly |
| Reference distributions for tests and power | `ChiSquared`, `FisherF`, `Beta`, `Poisson`, `Hypergeometric`, `NoncentralT`, `NoncentralChiSquared`, `NoncentralF` | constructors work |
| From data, from tables | `Empirical`, `FiniteDiscrete` | partly |
| Point mass, mixture, affine, truncated, censored | `PointMass`, `Mixture`, `Affine[D]`, `Truncated[D]`, `Censored[D]` | distribution functions work |
| Singular | `Cantor` | moments work |
| Families as records, constraints, exponential families, conjugacy, log-likelihood | `Family[D]`, `Constraint`, `ExpFamily[D]`, `Conjugate[P, D]`, `loglik` | works for normal; others partly |

The full catalogue (SciPy's ~100 continuous and ~20 discrete families,
checked against Distributions.jl's hierarchy and the Leemis relationship
chart) is generated after this design settles.

## multi

`Dimension`, `MvDensity`, `MvMass`, `MvSample`, `MvMoments`, `Marginal`;
`MvNormal` (Cholesky factor), `Dirichlet`, `Multinomial`; `Copula`,
`GaussianCopula`, `Joint[C]`. Typed; needs dense linear algebra.

## describe

| Concept | Type | Status |
|---|---|---|
| Mergeable accumulators | `Accumulator`, `MeanVar` | works |
| Covariance and correlation, shifted for large offsets | `Comoment` | works (1e-14 at offset 1e9) |
| Weighted mean and variance by weight kind | `WeightedMeanVar` | mean works |
| Skewness and kurtosis with explicit bias correction | `Moments4` | typed |
| Nine quantile definitions | `quantile`, `quantile_sorted` | works (matches NumPy for all nine) |
| Ranks with five tie rules | `ranks` | works (matches SciPy) |
| MAD with explicit scale, trimmed mean, winsorizing | `mad`, `MadScale`, `trimmed_mean`, `winsorize` | typed |
| Histogram with bin rules, out-of-range counts | `histogram`, `Bins`, `Histogram` | typed |
| Spearman, Kendall tau-b | `spearman`, `kendall_tau_b` | typed |

## fit

`MleResult[D]` with AIC, AICc, BIC; `fit_exp_family` (works for the
normal); `fit_mle` (typed, needs an optimiser); `Candidate`,
`akaike_weights` (works).

## infer

| Concept | Type | Status |
|---|---|---|
| Question: estimand × data shape, checked | `Estimand`, `Sample`, `Question` | works |
| Null model: observed statistic + null distribution | `NullModel` | p-values work (`Doubling`) |
| Test result (union of htest, SciPy, HypothesisTests.jl fields) | `TestResult` | typed |
| Interval with kind, estimate | `Interval`, `IntervalKind`, `Estimate` | typed |
| Assumptions as certificates | `Assumed[A]`, `Normality`, `EqualVariance` | works |
| Conventions | `TwoSidedRule`, `PValueKind`, `Exactness`, `Variances`, `Zeros`, `TieTolerance`, `LeveneCenter` | types work |
| t and z tests, from data or summaries | `t_test_one`, `t_test_two`, `t_test_paired`, `t_test_summary` | typed |
| Proportions and counts | `binomial_test`, `proportion_test`, `chi_square_goodness_of_fit`, `chi_square_independence`, `fisher_exact`, `mcnemar` | typed |
| Several groups and post hoc | `anova_one_way`, `levene`, `bartlett`, `post_hoc`, `PostHoc` | typed |
| Rank tests | `mann_whitney`, `wilcoxon_signed_rank`, `kruskal_wallis`, `friedman` | typed |
| Goodness of fit, normality | `ks_one`, `ks_two`, `shapiro_wilk`, `anderson_darling` | typed |
| Correlation, likelihood ratio | `correlation_test`, `likelihood_ratio` | typed |
| Multiple testing | `adjust`, `reject`, `Adjust` | works except Hommel (matches statsmodels) |
| Effect sizes | `Effect`, `EffectKind` | typed |
| Power, solving for one unknown | `PowerQuery`, `PowerDesign`, `solve_power` | query validation works |
| Resampling with user statistics | `Statistic`, `bootstrap`, `BootstrapInterval`, `permutation_test`, `PermutationScheme`, `monte_carlo_test`, `jackknife` | typed |

## Models

- **regress**: `Design` (matrix, intercept, offset, weights), `ols` by QR,
  `gls`, `glm` with `GlmFamily` and `Link` by IRLS, `penalized` (ridge,
  lasso, elastic net), `cross_validate`, `nonlinear_least_squares`;
  results carry diagnostics (leverage, Cook's distance, studentized
  residuals, VIF, Durbin–Watson, condition number) and an ANOVA table.
- **survival**: `kaplan_meier`, `nelson_aalen`, `log_rank`, `cox`
  (Breslow, Efron), `fit_censored` for parametric models.
- **timeseries**: `acf`, `pacf`, `ljung_box`, `periodogram` (needs FFT),
  `arima` with forecasts, `StateSpace` and `kalman_filter`, `adf_test`.
- **multivariate**: `pca` by SVD, `factor_analysis` with rotations,
  `discriminant` (LDA, QDA), `roc`, `kmeans`, `hierarchical`, `mds`,
  `cca`, `correspondence`.
- **bayes**: `metropolis`, `gibbs`, `Chain`, `diagnose` (split R-hat,
  effective size), `summarize` (credible intervals), `em_normal_mixture`.

All typed; most need dense linear algebra.

## Planned, not yet typed

| Concept | Expected shape | Depends on |
|---|---|---|
| Directional statistics (von Mises, Rayleigh test) | families and tests over `Circular` | — |
| Kernel density estimate | a `Density` + `Cdf` from data, bandwidth rule as an enum | — |
| Distribution from a characteristic function | FFT inversion into `FiniteDiscrete` or a gridded density | FFT module |
| Exact rank-test null distributions | polynomial products / FFT into `FiniteDiscrete` | FFT module |
| Matrix-variate (Wishart) | matrix-valued `MvSample` | matrix type |
| Inter-rater agreement (kappa) | over a rater matrix | — |
| Formula language and data frames | out of scope; MoonFrame on mooncakes | — |

## Open decisions

1. **Variate methods: explicit always, or a default?** Explicit methods
   make streams reproducible but every call longer. NEP 19 instead
   promises only the bit stream.
2. **Two-sided p-values by minimum likelihood** need the null's mass
   function, which `&Univariate` does not carry: add a mass closure to
   `NullModel`, or a discrete counterpart of `Univariate`.
3. **Discrete exponential-family fitting** needs an integer-data variant
   of `fit_exp_family`.
4. **Shifted accumulation** should also go into `MeanVar` and
   `WeightedMeanVar`, as in `Comoment`.
5. **Matrix type**: `DataMatrix` is a placeholder until a dense linear
   algebra module exists; the model layer's types will switch to it.
6. **Checked special functions**: names and which functions get a
   checked variant.

## Accuracy targets

Every implemented body is tested against an oracle: NumPy, SciPy,
statsmodels or mpmath for functions, published reference outputs for
generators, and NIST StRD certified values for summaries and regression.

## Language notes (MoonBit 0.10.21)

- Supertrait methods are not reachable with dot syntax through a trait
  object; call `Cdf::cdf(d, x)`.
- Trait methods implemented for a concrete type are not reachable with
  dot syntax, even in the same package, unless exposed with
  `pub extend T with Trait::{...}`.
- Traits have one `Self` and no associated types: families, exponential
  families and conjugate pairs are records of functions.
- `pub struct` cannot be constructed outside its package; `pub(all)` can.
  Library-produced results stay `pub`; caller-supplied data is `pub(all)`.
- `method` and `var` are reserved; an empty struct literal is `T::{ }`;
  a function in a record field is called `(f.make)(x)`.
- Adding a case to a `pub(all)` enum breaks callers' exhaustive matches, so
  shapes such as `Sample` are fixed up front.

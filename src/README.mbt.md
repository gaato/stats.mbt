# gaato/stats

Statistics foundations for MoonBit: probability distributions, special
functions, descriptive statistics and hypothesis tests, on Native, JS, Wasm
and Wasm GC.

> **Status**: pre-release; nothing is published yet. The numerical base,
> the special functions, the data shapes and a large part of the
> distributions, summaries and tests are implemented and checked against
> oracles. The model layers are signatures only. The [roadmap](#roadmap)
> lists both.

## Example

A distribution is a value. Each tail has its own formula, so the upper tail
keeps its digits where `1 - cdf` has none left:

```mbt check
///|
/// The backends differ in the last bits of an elementary function, so the
/// examples compare to a relative tolerance.
fn assert_close(got : Double, want : Double) -> Unit raise {
  assert_true((got - want).abs() <= 1.0e-10 * want.abs())
}

///|
test "both tails of a distribution" {
  let z = @dist.Normal::new()
  assert_eq(1.0 - z.cdf(10.0), 0.0)
  assert_close(z.sf(10.0), 7.619853024160526e-24)
  assert_close(z.logsf(40.0), -804.6084420137538)
  assert_close(z.isf(1.0e-100), 21.273453560965324)
  // Parameters as in SciPy: the shape, then loc and scale.
  let t = @dist.StudentT::new(3.0, loc=1.0, scale=2.0)
  assert_close(t.quantile(0.975), 7.3648926105674155)
  // Scale or rate is said in the call, not left to a default.
  assert_eq(
    @dist.Gamma::with_rate(2.0, rate=4.0),
    @dist.Gamma::new(2.0, scale=0.25),
  )
}
```

A summary is an accumulator that merges, and the conventions (the variance
divisor, the quantile definition) are arguments:

```mbt check
///|
test "summaries in pieces" {
  let xs = [2.0, 3.0, 5.0, 7.0, 11.0]
  let left : @describe.MeanVar = @describe.accumulate(xs[:2])
  let right : @describe.MeanVar = @describe.accumulate(xs[2:])
  let whole = left.merge(right)
  guard (whole.mean(), whole.variance(ddof=1)) is (Some(mean), Some(variance)) else {
    fail("five observations have a mean and a variance")
  }
  assert_close(mean, 5.6)
  assert_close(variance, 12.8)
  // Type7 is R's and NumPy's default, Type6 is Minitab's and SPSS's.
  assert_close(@describe.quantile(xs, 0.9, Type7), 9.4)
  assert_close(@describe.quantile(xs, 0.9, Type6), 11.0)
}
```

A test that relies on an assumption takes it as an argument, so it shows in
the caller's code:

```mbt check
///|
test "a Welch t test" {
  let a = [5.1, 4.9, 5.6, 5.8, 6.0]
  let b = [6.2, 6.4, 5.9, 6.8, 7.1]
  let normality : @infer.Assumed[@infer.Normality] = @infer.Assumed::postulate()
  let r = @infer.t_test_two(
    a,
    b,
    variances=Welch,
    alternative=TwoSided,
    level=0.95,
    normality~,
  )
  assert_close(r.statistic, -3.355780276070123)
  assert_close(r.p_value, 0.010001509249941531)
  guard (r.df, r.interval) is (Some(df), Some(interval)) else {
    fail("a t test has degrees of freedom and an interval")
  }
  assert_close(df, 7.995943945324385)
  assert_close(interval.lo, -1.6872344350767405)
  assert_close(interval.hi, -0.31276556492325946)
}
```

## Design

- **Distributions are the centre.** Each family is an immutable value; what a
  distribution can compute is expressed by capability traits (`Density`,
  `Cdf`, `Quantile`, `Moments`, sampling), with separate traits for
  continuous and discrete support.
- **Both tails are first class.** `sf`, `logcdf`, `logsf` and `isf` exist
  for every distribution. The upper tail has no `1 - cdf` default: each
  distribution states its own, and wrappers (location and scale,
  truncation, mixtures) forward every tail of what they wrap.
- **Parameters as in SciPy.** The catalogue has SciPy's 129 families and
  parameter names; every continuous family takes `loc` and `scale`.
  Discrete distributions answer on Doubles, so a quantile no integer type
  holds is still exact.
- **Conventions are explicit.** Scale versus rate, the sample quantile
  definition and the variance divisor are chosen by the caller, never by a
  default.
- **Summaries compose.** Descriptive statistics are accumulators with an
  identity and an associative `merge`, so data can be summarised in pieces,
  in parallel or as a stream.
- **Assumptions are visible.** A test that relies on an assumption asks for
  a certificate value, which the caller obtains explicitly; a choice that
  brings an assumption (pooled variances) carries its certificate, so it
  cannot be made without one.
- **No global state.** Sampling takes the random source as an argument.
- **Unfinished is never approximate.** A path not written yet fails and
  says so; it does not return a rougher answer.
- **Accuracy is measured.** Every implemented body is tested against an
  oracle on all four backends: mpmath at high precision for numerics,
  SciPy and NumPy for conventions, published reference outputs for the
  generator. [DESIGN.md](DESIGN.md) lists the worst error measured for
  each special function.

## Packages

Each package imports only packages above it in this table.

| Package | Contents |
|---|---|
| `num` | Log-space arithmetic, compensated sums, solvers, random sources |
| `special` | Incomplete beta and gamma, normal tails, gamma family |
| `data` | Observations, data shapes, weights, tables, conventions |
| `describe` | Mergeable accumulators, quantiles, ranks, robust summaries |
| `dist` | Distributions, families and parameters, composition |
| `multi` | Multivariate distributions and copulas |
| `fit` | Maximum likelihood, information criteria |
| `infer` | Questions, null models, tests, multiple testing, power, resampling |
| `regress`, `survival`, `timeseries`, `multivariate`, `bayes` | Models |

## Roadmap

A ticked item is implemented and tested against an oracle. An unticked one
is not written yet; most have their signatures in place, and calling one
fails rather than approximates. [DESIGN.md](DESIGN.md) has the type and
status of every concept.

- [ ] `num`
    - [x] Log-space arithmetic, compensated summation
    - [x] Bracketed root finding, continued fractions, series, quadrature
    - [x] Xoshiro256; uniforms, indices, shuffling
    - [ ] Normal, exponential and gamma variates by named methods
    - [ ] Generator jumps for parallel streams
    - [ ] Checked variants raising `DomainError`
- [ ] `special`
    - [x] Error functions, normal tails and their inverses
    - [x] ln Γ, digamma, trigamma, ln B, ln of binomial coefficients
    - [x] Incomplete gamma and beta, both tails, with inverses
    - [x] Owen's T, Kolmogorov tail
    - [ ] Incomplete gamma at very large parameters (NaN until then)
    - [ ] Studentized range
- [x] `data`: observations with censoring, weights by kind, tables, pairs,
      groups, blocks, published summaries, survival records
- [ ] `describe`
    - [x] Mergeable accumulators: mean, variance, extremes, covariance
    - [x] Sum, mean, variance, standard deviation, median of an array
    - [x] Quantiles by nine definitions; ranks by five tie rules
    - [ ] Weighted variance (the weighted mean works)
    - [ ] Skewness and kurtosis
    - [ ] MAD, trimmed mean, winsorizing
    - [ ] Histograms
    - [ ] Spearman, Kendall tau-b
- [ ] `dist`
    - [x] Capability traits with both tails and their logarithms
    - [x] All 129 families as types, with checked constructors and
          sampling by inversion
    - [x] 44 families with closed forms
    - [x] Normal, t, χ², gamma, beta, F, binomial, Poisson, hypergeometric
    - [ ] The other 76 families
    - [ ] Moments, mode and entropy for every family
    - [ ] Characteristic and moment generating functions beyond a few
    - [x] Empirical, finite discrete and atomic distributions
    - [x] Point mass, mixture, affine, truncated, censored
    - [ ] Exponential families and conjugacy beyond the normal
- [ ] `fit`
    - [x] Exponential-family fit for the normal
    - [x] AIC, AICc, BIC, Akaike weights
    - [ ] General maximum likelihood (needs an optimiser)
    - [ ] Exponential-family fit for discrete data
- [ ] `infer`
    - [x] Questions, estimates, null models; every alternative and either
          two-sided rule
    - [x] t tests: one sample, Welch, pooled, paired, from summaries
    - [x] Exact binomial, χ² goodness of fit and independence, Fisher 2 × 2
    - [x] One-way ANOVA (classic, Welch), Levene, Bartlett, Kruskal–Wallis
    - [x] Pearson, Spearman, likelihood ratio
    - [x] Mann–Whitney (exact without ties, or normal)
    - [x] Bootstrap (percentile, basic), permutation and Monte Carlo
          tests, one-sample jackknife
    - [x] Multiple-testing adjustments other than Hommel
    - [ ] Hommel
    - [ ] BCa intervals, block bootstrap, exact pairings, jackknife beyond
          one sample
    - [ ] Fisher beyond 2 × 2; exact Mann–Whitney with ties
    - [ ] Proportion z, McNemar
    - [ ] Wilcoxon signed-rank, Friedman, Kendall
    - [ ] Kolmogorov–Smirnov, Shapiro–Wilk, Anderson–Darling
    - [ ] Post hoc comparisons (need the studentized range)
    - [ ] Power and sample size
- [ ] On dense linear algebra
    - [ ] A matrix type to replace the `DataMatrix` placeholder
    - [ ] `multi`: multivariate normal, Dirichlet, multinomial, copulas
    - [ ] `regress`: least squares, GLS, GLM, penalised and nonlinear
          fits, diagnostics
    - [ ] `survival`: Kaplan–Meier, Nelson–Aalen, log-rank, Cox, parametric
    - [ ] `multivariate`: PCA, factor analysis, discriminants, clustering,
          MDS, CCA, correspondence analysis, ROC
    - [ ] `bayes`: Metropolis, Gibbs, chain diagnostics, EM
- [ ] On an FFT (`gaato/fft`, from transforms.mbt)
    - [ ] `timeseries`: ACF and PACF, Ljung–Box, ADF, periodogram, ARIMA,
          Kalman filter
    - [ ] Exact null distributions of rank statistics
    - [ ] Distributions from a characteristic function
- [ ] Not yet typed
    - [ ] Directional statistics
    - [ ] Kernel density estimates
    - [ ] Matrix-variate distributions
    - [ ] Inter-rater agreement

Not planned: a formula language and data frames.

## Development

The toolchain is MoonBit 0.10.21, pinned in `moonbit-version`.

```sh
scripts/gates.sh
```

runs the checks every change must pass: formatting, type-checking and the
tests on all four backends, and that the public interface did not move
unreviewed. This file is one of the tests: its examples are compiled and
run. The oracle tests and the catalogue are generated by the other scripts
in `scripts/`.

## License

Apache-2.0

///|
/// The tests, by family. Signatures only for now; each takes the data
/// shape from `data`, every convention from `conventions.mbt`, and the
/// assumptions it relies on as certificates.

// ---------------------------------------------------------------------------
// Location: t and z tests.
// ---------------------------------------------------------------------------

///|
pub fn t_test_one(
  x : ArrayView[Double],
  mu~ : Double,
  alternative~ : Alternative,
  level~ : Double,
  normality~ : Assumed[Normality],
) -> TestResult raise QuestionError {
  ignore((x, mu, alternative, level, normality))
  ...
}

///|
/// Two independent samples; `Pooled` additionally asks for the
/// equal-variance certificate.
pub fn t_test_two(
  x : ArrayView[Double],
  y : ArrayView[Double],
  variances~ : Variances,
  alternative~ : Alternative,
  level~ : Double,
  normality~ : Assumed[Normality],
  equal_variance? : Assumed[EqualVariance],
) -> TestResult raise QuestionError {
  ignore((x, y, variances, alternative, level, normality, equal_variance))
  ...
}

///|
pub fn t_test_paired(
  x : ArrayView[Double],
  y : ArrayView[Double],
  alternative~ : Alternative,
  level~ : Double,
  normality~ : Assumed[Normality],
) -> TestResult raise QuestionError {
  ignore((x, y, alternative, level, normality))
  ...
}

///|
/// From published summaries only.
pub fn t_test_summary(
  x : @data.Summary,
  y : @data.Summary,
  variances~ : Variances,
  alternative~ : Alternative,
  level~ : Double,
  normality~ : Assumed[Normality],
) -> TestResult raise QuestionError {
  ignore((x, y, variances, alternative, level, normality))
  ...
}

// ---------------------------------------------------------------------------
// Proportions and counts.
// ---------------------------------------------------------------------------

///|
/// Exact binomial test.
pub fn binomial_test(
  x : @data.Proportion,
  p0~ : Double,
  alternative~ : Alternative,
  two_sided~ : TwoSidedRule,
  level~ : Double,
) -> TestResult raise QuestionError {
  ignore((x, p0, alternative, two_sided, level))
  ...
}

///|
/// Score (z) test for one or several proportions.
pub fn proportion_test(
  xs : ArrayView[@data.Proportion],
  p0? : Double,
  continuity~ : Bool,
  alternative~ : Alternative,
  level~ : Double,
) -> TestResult raise QuestionError {
  ignore((xs, p0, continuity, alternative, level))
  ...
}

///|
pub fn chi_square_goodness_of_fit(
  counts : @data.CategoryCounts,
  fitted_params~ : Int,
) -> TestResult raise QuestionError {
  ignore((counts, fitted_params))
  ...
}

///|
pub fn chi_square_independence(
  table : @data.Table,
  continuity~ : Bool,
) -> TestResult raise QuestionError {
  ignore((table, continuity))
  ...
}

///|
/// Fisher's exact test for r × c tables (network algorithm beyond 2 × 2).
pub fn fisher_exact(
  table : @data.Table,
  alternative~ : Alternative,
  two_sided~ : TwoSidedRule,
) -> TestResult raise QuestionError {
  ignore((table, alternative, two_sided))
  ...
}

///|
pub fn mcnemar(table : @data.Table, exactness~ : Exactness) -> TestResult raise QuestionError {
  ignore((table, exactness))
  ...
}

// ---------------------------------------------------------------------------
// Several groups.
// ---------------------------------------------------------------------------

///|
/// One-way ANOVA; `Welch` gives Welch's heteroscedastic version.
pub fn anova_one_way(
  groups : ArrayView[ArrayView[Double]],
  variances~ : Variances,
  normality~ : Assumed[Normality],
) -> TestResult raise QuestionError {
  ignore((groups, variances, normality))
  ...
}

///|
pub fn levene(
  groups : ArrayView[ArrayView[Double]],
  center~ : LeveneCenter,
) -> TestResult raise QuestionError {
  ignore((groups, center))
  ...
}

///|
/// Levene's original uses the mean; Brown–Forsythe the median.
#derive(Eq, Debug)
pub(all) enum LeveneCenter {
  Mean
  Median
  Trimmed(Double)
}

///|
pub fn bartlett(
  groups : ArrayView[ArrayView[Double]],
  normality~ : Assumed[Normality],
) -> TestResult raise QuestionError {
  ignore((groups, normality))
  ...
}

///|
/// Pairwise comparisons after ANOVA.
#derive(Eq, Debug)
pub(all) enum PostHoc {
  TukeyHsd
  Dunnett(control~ : Int)
  /// Pairwise t tests with a p-value adjustment.
  Pairwise(Adjust)
}

///|
pub fn post_hoc(
  groups : ArrayView[ArrayView[Double]],
  how : PostHoc,
  level~ : Double,
) -> Array[TestResult] raise QuestionError {
  ignore((groups, how, level))
  ...
}

// ---------------------------------------------------------------------------
// Rank tests.
// ---------------------------------------------------------------------------

///|
pub fn mann_whitney(
  x : ArrayView[Double],
  y : ArrayView[Double],
  alternative~ : Alternative,
  exactness~ : Exactness,
  continuity~ : Bool,
) -> TestResult raise QuestionError {
  ignore((x, y, alternative, exactness, continuity))
  ...
}

///|
pub fn wilcoxon_signed_rank(
  x : ArrayView[Double],
  y? : ArrayView[Double],
  alternative~ : Alternative,
  zeros~ : Zeros,
  exactness~ : Exactness,
  continuity~ : Bool,
) -> TestResult raise QuestionError {
  ignore((x, y, alternative, zeros, exactness, continuity))
  ...
}

///|
pub fn kruskal_wallis(
  groups : ArrayView[ArrayView[Double]],
) -> TestResult raise QuestionError {
  ignore(groups)
  ...
}

///|
pub fn friedman(blocks : @data.Blocks) -> TestResult raise QuestionError {
  ignore(blocks)
  ...
}

// ---------------------------------------------------------------------------
// Goodness of fit and normality.
// ---------------------------------------------------------------------------

///|
/// One-sample Kolmogorov–Smirnov against any distribution with a CDF.
/// The null must not have been fitted to the same data (use Lilliefors
/// then); the certificate makes that visible.
pub fn ks_one(
  x : ArrayView[Double],
  null_dist : &@dist.Univariate,
  alternative~ : Alternative,
  exactness~ : Exactness,
) -> TestResult raise QuestionError {
  ignore((x, null_dist, alternative, exactness))
  ...
}

///|
pub fn ks_two(
  x : ArrayView[Double],
  y : ArrayView[Double],
  alternative~ : Alternative,
  exactness~ : Exactness,
) -> TestResult raise QuestionError {
  ignore((x, y, alternative, exactness))
  ...
}

///|
pub fn shapiro_wilk(x : ArrayView[Double]) -> TestResult raise QuestionError {
  ignore(x)
  ...
}

///|
pub fn anderson_darling(
  x : ArrayView[Double],
  null_dist : &@dist.Univariate,
) -> TestResult raise QuestionError {
  ignore((x, null_dist))
  ...
}

// ---------------------------------------------------------------------------
// Correlation and likelihood.
// ---------------------------------------------------------------------------

///|
#derive(Eq, Debug)
pub(all) enum CorrelationKind {
  Pearson
  Spearman
  KendallTauB
}

///|
/// Test of zero correlation with an interval (Fisher z for Pearson).
pub fn correlation_test(
  xy : @data.Bivariate,
  kind : CorrelationKind,
  alternative~ : Alternative,
  level~ : Double,
) -> TestResult raise QuestionError {
  ignore((xy, kind, alternative, level))
  ...
}

///|
/// Likelihood ratio test of nested models: 2(ℓ₁ − ℓ₀) against χ² with the
/// difference in free parameters.
pub fn likelihood_ratio(
  loglik_null : Double,
  loglik_alt : Double,
  df : Int,
) -> TestResult raise QuestionError {
  ignore((loglik_null, loglik_alt, df))
  ...
}

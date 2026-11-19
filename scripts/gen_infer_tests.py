"""Generate src/infer/infer_oracle_test.mbt: the tests against SciPy.

Fixed random data sets; every statistic, p-value, degrees of freedom and
interval SciPy reports is compared at 1e-10 relative (SciPy's own
p-values come from the same kind of special functions, so the two agree
far below that, but not to the last bit).

    uv run --no-project --with scipy python -I scripts/gen_infer_tests.py
"""

import math
import pathlib

import numpy as np
import scipy
import scipy.stats as st

ROOT = pathlib.Path(__file__).resolve().parent.parent
rng = np.random.default_rng(20261119)


def lit(x):
    r = repr(float(x))
    if r == "inf":
        return "1.0 / 0.0"
    if r == "-inf":
        return "-1.0 / 0.0"
    if "e" in r and "." not in r.split("e")[0]:
        m, e = r.split("e")
        r = m + ".0e" + e
    elif "." not in r and "e" not in r:
        r += ".0"
    return r


def arr(xs):
    return "[" + ", ".join(lit(v) for v in xs) + "]"


ALT = {"two-sided": "TwoSided", "greater": "Greater", "less": "Less"}
TESTS = []


def block(name, body):
    TESTS.append(f"///|\ntest \"infer: {name}\" {{\n" + "\n".join("  " + l for l in body) + "\n}")


def checks(result_var, stat=None, p=None, df=None, df2=None, ci=None, est=None):
    out = []
    if stat is not None:
        out.append(f"assert_close({result_var}.statistic, {lit(stat)}, \"statistic\")")
    if p is not None:
        out.append(f"assert_close({result_var}.p_value, {lit(p)}, \"p\")")
    if df is not None:
        out.append(f"assert_close({result_var}.df.unwrap(), {lit(df)}, \"df\")")
    if df2 is not None:
        out.append(f"assert_close({result_var}.df2.unwrap(), {lit(df2)}, \"df2\")")
    if est is not None:
        out.append(f"assert_close({result_var}.estimate.unwrap(), {lit(est)}, \"estimate\")")
    if ci is not None:
        out.append(f"let iv = {result_var}.interval.unwrap()")
        out.append(f"assert_close(iv.lo, {lit(ci[0])}, \"ci lo\")")
        out.append(f"assert_close(iv.hi, {lit(ci[1])}, \"ci hi\")")
    return out


HELPERS = '''///|
fn assert_close(got : Double, want : Double, what : String) -> Unit raise {
  if got == want {
    return
  }
  let err = (got - want).abs() / (want.abs() + 1.0e-300)
  assert_true(err <= 1.0e-10, msg="\\{what}: \\{got} vs \\{want} (relative \\{err})")
}

///|
fn normality() -> @infer.Assumed[@infer.Normality] {
  @infer.Assumed::postulate()
}
'''

x = rng.normal(5, 2, 15).round(3)
y = rng.normal(6, 3, 12).round(3)
z = (x + rng.normal(0.5, 1, 15)).round(3)
data = [f"let x : Array[Double] = {arr(x)}", f"let y : Array[Double] = {arr(y)}", f"let z : Array[Double] = {arr(z)}"]

for alt in ALT:
    r = st.ttest_1samp(x, 4.5, alternative=alt)
    ci = r.confidence_interval(0.9)
    block(f"one-sample t ({alt})", data + [
        f"let r = @infer.t_test_one(x[:], mu=4.5, alternative={ALT[alt]}, level=0.9, normality=normality())",
    ] + checks("r", r.statistic, r.pvalue, r.df, ci=(ci.low, ci.high)))
    for eq, var in ((False, "Welch"), (True, "Pooled")):
        r = st.ttest_ind(x, y, equal_var=eq, alternative=alt)
        ci = r.confidence_interval(0.95)
        extra = ", equal_variance=@infer.Assumed::postulate()" if eq else ""
        block(f"two-sample t {var} ({alt})", data + [
            f"let r = @infer.t_test_two(x[:], y[:], variances={var}, alternative={ALT[alt]}, level=0.95, normality=normality(){extra})",
        ] + checks("r", r.statistic, r.pvalue, r.df, ci=(ci.low, ci.high)))
    r = st.ttest_rel(x, z, alternative=alt)
    ci = r.confidence_interval(0.95)
    block(f"paired t ({alt})", data + [
        f"let r = @infer.t_test_paired(x[:], z[:], alternative={ALT[alt]}, level=0.95, normality=normality())",
    ] + checks("r", r.statistic, r.pvalue, r.df, ci=(ci.low, ci.high)))

r = st.ttest_ind_from_stats(5.1, 1.9, 20, 6.3, 2.7, 25, equal_var=False)
block("t from summaries (Welch)", [
    "let r = @infer.t_test_summary({ n: 20, mean: 5.1, sd: 1.9 }, { n: 25, mean: 6.3, sd: 2.7 }, variances=Welch, alternative=TwoSided, level=0.95, normality=normality())",
] + checks("r", r.statistic, r.pvalue))

block("pooled t without the certificate is an error", data + [
    "let refused = try @infer.t_test_two(x[:], y[:], variances=Pooled, alternative=TwoSided, level=0.95, normality=normality()) catch {",
    "  MissingAssumption(name~) => name == \"EqualVariance\"",
    "  _ => false",
    "} noraise {",
    "  _ => false",
    "}",
    "assert_true(refused)",
])

for k, n, p0 in ((7, 20, 0.3), (2, 25, 0.35), (40, 50, 0.6)):
    for alt in ALT:
        r = st.binomtest(k, n, p0, alternative=alt)
        ci = r.proportion_ci(0.95, method="exact")
        block(f"binomial k={k} n={n} p={p0} ({alt})", [
            f"let r = @infer.binomial_test({{ successes: {k}, trials: {n} }}, p0={lit(p0)}, alternative={ALT[alt]}, two_sided=MinLikelihood, level=0.95)",
        ] + checks("r", p=r.pvalue, ci=(ci.low, ci.high)))

obs = [18, 22, 30, 30]
r = st.chisquare(obs)
block("chi-square goodness of fit, uniform", [
    f"let r = @infer.chi_square_goodness_of_fit({{ counts: [{', '.join(f'{c}L' for c in obs)}], null_probs: None }}, fitted_params=0)",
] + checks("r", r.statistic, r.pvalue))
probs = [0.1, 0.2, 0.3, 0.4]
r = st.chisquare(obs, f_exp=np.array(probs) * sum(obs), ddof=1)
block("chi-square goodness of fit, probabilities, one fitted", [
    f"let r = @infer.chi_square_goodness_of_fit({{ counts: [{', '.join(f'{c}L' for c in obs)}], null_probs: Some({arr(probs)}) }}, fitted_params=1)",
] + checks("r", r.statistic, r.pvalue))

for table, corr in (([[12, 5], [7, 16]], True), ([[12, 5], [7, 16]], False), ([[10, 20, 30], [15, 15, 25]], True)):
    res = st.chi2_contingency(table, correction=corr)
    flat = [c for row in table for c in row]
    block(f"chi-square independence {len(table)}x{len(table[0])} correction={corr}", [
        f"let t = @data.Table::new({len(table)}, {len(table[0])}, [{', '.join(f'{c}L' for c in flat)}][:])",
        f"let r = @infer.chi_square_independence(t, continuity={'true' if corr else 'false'})",
    ] + checks("r", res.statistic, res.pvalue, res.dof))

for table in ([[8, 2], [1, 5]], [[3, 9], [12, 4]], [[20, 15], [11, 30]]):
    for alt in ALT:
        res = st.fisher_exact(table, alternative=alt)
        flat = [c for row in table for c in row]
        block(f"fisher exact {flat} ({alt})", [
            f"let t = @data.Table::new(2, 2, [{', '.join(f'{c}L' for c in flat)}][:])",
            f"let r = @infer.fisher_exact(t, alternative={ALT[alt]}, two_sided=MinLikelihood)",
        ] + checks("r", res.statistic, res.pvalue))

g3 = rng.normal(5.5, 2.5, 10).round(3)
data3 = data + [f"let w : Array[Double] = {arr(g3)}"]
r = st.f_oneway(x, y, g3)
block("one-way ANOVA", data3 + [
    "let r = @infer.anova_one_way([x[:], y[:], w[:]][:], variances=Pooled, normality=normality())",
] + checks("r", r.statistic, r.pvalue))
r = st.f_oneway(x, y, g3, equal_var=False)
block("Welch's ANOVA", data3 + [
    "let r = @infer.anova_one_way([x[:], y[:], w[:]][:], variances=Welch, normality=normality())",
] + checks("r", r.statistic, r.pvalue))
for center, mb in (("mean", "Mean"), ("median", "Median"), ("trimmed", "Trimmed(0.1)")):
    r = st.levene(x, y, g3, center=center, proportiontocut=0.1)
    block(f"Levene ({center})", data3 + [
        f"let r = @infer.levene([x[:], y[:], w[:]][:], center={mb})",
    ] + checks("r", r.statistic, r.pvalue))
r = st.bartlett(x, y, g3)
block("Bartlett", data3 + [
    "let r = @infer.bartlett([x[:], y[:], w[:]][:], normality=normality())",
] + checks("r", r.statistic, r.pvalue))
r = st.kruskal(x, y, g3)
block("Kruskal-Wallis", data3 + [
    "let r = @infer.kruskal_wallis([x[:], y[:], w[:]][:])",
] + checks("r", r.statistic, r.pvalue))

for alt in ALT:
    r = st.pearsonr(x, z, alternative=alt)
    ci = r.confidence_interval(0.95)
    block(f"Pearson ({alt})", data + [
        "let xy = @data.Bivariate::new(x[:], z[:])",
        f"let r = @infer.correlation_test(xy, Pearson, alternative={ALT[alt]}, level=0.95)",
        f"assert_close(r.estimate.unwrap(), {lit(r.statistic)}, \"r\")",
        f"assert_close(r.p_value, {lit(r.pvalue)}, \"p\")",
        "let iv = r.interval.unwrap()",
        f"assert_close(iv.lo, {lit(ci.low)}, \"ci lo\")" if math.isfinite(ci.low) else "assert_true(iv.lo == -1.0 || iv.lo == -1.0 / 0.0 || iv.lo <= -0.9999999999)",
        f"assert_close(iv.hi, {lit(ci.high)}, \"ci hi\")" if math.isfinite(ci.high) else "assert_true(iv.hi >= 0.9999999999)",
    ])
    r = st.spearmanr(x, z, alternative=alt)
    block(f"Spearman ({alt})", data + [
        "let xy = @data.Bivariate::new(x[:], z[:])",
        f"let r = @infer.correlation_test(xy, Spearman, alternative={ALT[alt]}, level=0.95)",
        f"assert_close(r.estimate.unwrap(), {lit(r.statistic)}, \"rho\")",
        f"assert_close(r.p_value, {lit(r.pvalue)}, \"p\")",
    ])

a = rng.normal(0, 1, 7).round(4)
b = rng.normal(0.8, 1, 6).round(4)
ab = [f"let a : Array[Double] = {arr(a)}", f"let b : Array[Double] = {arr(b)}"]
for alt in ALT:
    r = st.mannwhitneyu(a, b, alternative=alt, method="exact")
    block(f"Mann-Whitney exact ({alt})", ab + [
        f"let r = @infer.mann_whitney(a[:], b[:], alternative={ALT[alt]}, exactness=Exact, continuity=false)",
    ] + checks("r", r.statistic, r.pvalue))
    r = st.mannwhitneyu(x, y, alternative=alt, method="asymptotic", use_continuity=True)
    block(f"Mann-Whitney asymptotic ({alt})", data + [
        f"let r = @infer.mann_whitney(x[:], y[:], alternative={ALT[alt]}, exactness=Asymptotic, continuity=true)",
    ] + checks("r", r.statistic, r.pvalue))

r = st.permutation_test((a, b), lambda u, v: np.mean(u) - np.mean(v), permutation_type="independent", n_resamples=np.inf, vectorized=False, alternative="two-sided")
block("exact permutation test of a mean difference", ab + [
    "let q = @infer.Question::new(MeanDifference, Two(a[:], b[:]))",
    "let rng = @num.Xoshiro256::seed(1)",
    "let stat = (s : @infer.Sample) => match s {",
    "  Two(u, v) => {",
    "    let mut su = 0.0",
    "    let mut sv = 0.0",
    "    for t in u { su = su + t }",
    "    for t in v { sv = sv + t }",
    "    su / u.length().to_double() - sv / v.length().to_double()",
    "  }",
    "  _ => 0.0",
    "}",
    "ignore(q)",
    "let m = @infer.permutation_test(Two(a[:], b[:]), stat, scheme=Independent, reps=None, ties={ relative: 1.0e-12 }, rng~)",
    f"assert_close(m.p_value(TwoSided), {lit(r.pvalue)}, \"p\")",
])

block("likelihood ratio", [
    "let r = @infer.likelihood_ratio(-120.5, -115.25, 2)",
] + checks("r", 10.5, st.chi2.sf(10.5, 2)))

header = f"// Generated by scripts/gen_infer_tests.py from SciPy {scipy.__version__}; do not edit.\n\n"
(ROOT / "src/infer/infer_oracle_test.mbt").write_text(header + HELPERS + "\n" + "\n\n".join(TESTS) + "\n")
print(len(TESTS), "tests")

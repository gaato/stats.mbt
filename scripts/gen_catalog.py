"""Generate the distribution catalogue from SciPy's metadata.

Writes src/dist/catalog.mbt (types, constructors, families, stubs) and
src/dist/catalog_test.mbt (constructor checks against SciPy's _argcheck,
and SciPy oracle values for the families implemented by hand in
src/dist/closed_forms.mbt).

Run from the repository root:
    uv run --no-project --with scipy --with mpmath python -I scripts/gen_catalog.py
"""

import math
import pathlib

import numpy as np
import scipy
import scipy.stats as st
from scipy.stats._distr_params import distcont, distdiscrete

ROOT = pathlib.Path(__file__).resolve().parent.parent

# SciPy names not generated.
EXISTING = set()
# Families that keep SciPy's location and scale as their own parameters,
# because the capabilities written on top of them (exponential-family
# fitting, conjugate updates, gamma by rate) need them; the rest take
# location and scale from `Affine`.
BUILTIN = {"norm": ("loc", "scale"), "gamma": ("scale",)}
# Aliases of another SciPy name.
ALIASES = {"reciprocal": "loguniform"}

NAMES = {
    "norm": "Normal", "t": "StudentT", "gamma": "Gamma", "beta": "Beta",
    "chi2": "ChiSquared", "f": "FisherF", "nct": "NoncentralT",
    "ncx2": "NoncentralChiSquared", "ncf": "NoncentralF", "binom": "Binomial",
    "poisson": "Poisson", "hypergeom": "Hypergeometric",
    "alpha": "Alpha", "anglit": "Anglit", "arcsine": "Arcsine", "argus": "Argus",
    "betaprime": "BetaPrime", "bradford": "Bradford", "burr": "Burr",
    "burr12": "Burr12", "cauchy": "Cauchy", "chi": "Chi", "cosine": "Cosine",
    "crystalball": "CrystalBall", "dgamma": "DoubleGamma",
    "dpareto_lognorm": "DoubleParetoLognormal", "dweibull": "DoubleWeibull",
    "erlang": "Erlang", "expon": "Exponential",
    "exponnorm": "ExponentiallyModifiedNormal", "exponpow": "ExponentialPower",
    "exponweib": "ExponentiatedWeibull", "fatiguelife": "FatigueLife",
    "fisk": "Fisk", "foldcauchy": "FoldedCauchy", "foldnorm": "FoldedNormal",
    "gausshyper": "GaussHypergeometric", "genexpon": "GeneralizedExponential",
    "genextreme": "GeneralizedExtremeValue", "gengamma": "GeneralizedGamma",
    "genhalflogistic": "GeneralizedHalfLogistic",
    "genhyperbolic": "GeneralizedHyperbolic",
    "geninvgauss": "GeneralizedInverseGaussian",
    "genlogistic": "GeneralizedLogistic", "gennorm": "GeneralizedNormal",
    "halfgennorm": "HalfGeneralizedNormal", "genpareto": "GeneralizedPareto",
    "gibrat": "Gibrat", "gompertz": "Gompertz", "gumbel_l": "GumbelLeft",
    "gumbel_r": "GumbelRight", "halfcauchy": "HalfCauchy",
    "halflogistic": "HalfLogistic", "halfnorm": "HalfNormal",
    "hypsecant": "HyperbolicSecant", "invgamma": "InverseGamma",
    "invgauss": "InverseGaussian", "invweibull": "InverseWeibull",
    "irwinhall": "IrwinHall", "jf_skew_t": "JonesFaddySkewT",
    "johnsonsb": "JohnsonSB", "johnsonsu": "JohnsonSU", "kappa4": "Kappa4",
    "kappa3": "Kappa3", "ksone": "KolmogorovSmirnovOne",
    "kstwo": "KolmogorovSmirnovTwo", "kstwobign": "Kolmogorov",
    "landau": "Landau", "laplace": "Laplace",
    "laplace_asymmetric": "AsymmetricLaplace", "levy": "Levy",
    "levy_l": "LevyLeft", "levy_stable": "Stable", "loggamma": "LogGamma",
    "logistic": "Logistic", "loglaplace": "LogLaplace", "lognorm": "LogNormal",
    "loguniform": "LogUniform", "lomax": "Lomax", "maxwell": "Maxwell",
    "mielke": "Mielke", "moyal": "Moyal", "nakagami": "Nakagami",
    "norminvgauss": "NormalInverseGaussian", "pareto": "Pareto",
    "pearson3": "PearsonIII", "powerlaw": "PowerLaw",
    "powerlognorm": "PowerLogNormal", "powernorm": "PowerNormal",
    "rayleigh": "Rayleigh", "rdist": "RDist",
    "recipinvgauss": "ReciprocalInverseGaussian",
    "rel_breitwigner": "RelativisticBreitWigner", "rice": "Rice",
    "semicircular": "Semicircular", "skewcauchy": "SkewCauchy",
    "skewnorm": "SkewNormal", "studentized_range": "StudentizedRange",
    "trapezoid": "Trapezoid", "triang": "Triangular",
    "truncexpon": "TruncatedExponential", "truncnorm": "TruncatedNormal",
    "truncpareto": "TruncatedPareto", "truncweibull_min": "TruncatedWeibull",
    "tukeylambda": "TukeyLambda", "uniform": "Uniform", "vonmises": "VonMises",
    "vonmises_line": "VonMisesLine", "wald": "Wald", "weibull_max": "WeibullMax",
    "weibull_min": "Weibull", "wrapcauchy": "WrappedCauchy",
    # discrete
    "bernoulli": "Bernoulli", "betabinom": "BetaBinomial",
    "betanbinom": "BetaNegativeBinomial", "boltzmann": "Boltzmann",
    "dlaplace": "DiscreteLaplace", "geom": "Geometric",
    "nchypergeom_fisher": "FisherNoncentralHypergeometric",
    "nchypergeom_wallenius": "WalleniusNoncentralHypergeometric",
    "logser": "LogSeries", "nbinom": "NegativeBinomial", "planck": "Planck",
    "randint": "DiscreteUniform", "skellam": "Skellam", "zipf": "Zipf",
    "zipfian": "Zipfian", "yulesimon": "YuleSimon",
    "nhypergeom": "NegativeHypergeometric",
}

# Field names where SciPy's are uppercase or ambiguous.
FIELD_RENAMES = {
    ("nchypergeom_fisher", "M"): "population",
    ("nchypergeom_fisher", "n"): "successes",
    ("nchypergeom_fisher", "N"): "draws",
    ("nchypergeom_wallenius", "M"): "population",
    ("nchypergeom_wallenius", "n"): "successes",
    ("nchypergeom_wallenius", "N"): "draws",
    ("hypergeom", "M"): "population",
    ("hypergeom", "n"): "successes",
    ("hypergeom", "N"): "draws",
    ("nhypergeom", "M"): "population",
    ("nhypergeom", "n"): "successes",
    ("nhypergeom", "r"): "failures",
    ("boltzmann", "N"): "n",
    ("exponnorm", "K"): "k",
}

# Cross-parameter constraints SciPy checks in _argcheck beyond the boxes.
EXTRA = {
    "genhyperbolic": ("a > b.abs()", "a", "need |b| < a"),
    "loguniform": ("a < b", "b", "need a < b"),
    "norminvgauss": ("a > b.abs()", "a", "need |b| < a"),
    "trapezoid": ("c <= d", "d", "need c <= d"),
    "truncnorm": ("a < b", "b", "need a < b"),
    "truncweibull_min": ("a < b", "b", "need a < b"),
    "boltzmann": ("n > 0", "n", "need n > 0"),
    "hypergeom": (
        "successes <= population && draws <= population", "successes",
        "need successes, draws <= population",
    ),
    "nchypergeom_fisher": (
        "successes <= population && draws <= population", "successes",
        "need successes, draws <= population",
    ),
    "nchypergeom_wallenius": (
        "successes <= population && draws <= population", "successes",
        "need successes, draws <= population",
    ),
    "nbinom": ("n > 0 && p > 0.0", "n", "need n > 0 and p > 0"),
    "randint": ("low < high", "high", "need low < high"),
    "zipfian": ("n > 0", "n", "need n > 0"),
    "nhypergeom": (
        "successes <= population && failures <= population - successes",
        "successes", "need successes <= population, failures <= population - successes",
    ),
}

# Implemented by hand in src/dist/closed_forms.mbt.
HAND_CONTINUOUS = {
    "uniform", "expon", "laplace", "logistic", "cauchy", "gumbel_r", "gumbel_l",
    "weibull_min", "weibull_max", "pareto", "lomax", "rayleigh", "triang",
    "arcsine", "halfcauchy", "halflogistic", "hypsecant", "genlogistic",
    "genpareto", "genextreme", "burr", "burr12", "fisk", "powerlaw",
    "bradford", "loguniform", "anglit", "truncexpon", "exponpow", "exponweib",
    "gompertz", "kappa3", "dweibull", "truncpareto", "loglaplace",
    "laplace_asymmetric", "skewcauchy",
}
HAND_DISCRETE = {
    "bernoulli", "geom", "randint", "dlaplace", "planck", "boltzmann", "zipfian",
}

# Required methods written by hand (closed_forms.mbt, extras.mbt); the
# generator emits a stub only for the others.
HAND_METHODS = {name: {"pdf", "cdf", "quantile"} for name in HAND_CONTINUOUS}
HAND_METHODS.update({name: {"pmf", "cdf_int", "quantile_int"} for name in HAND_DISCRETE})
HAND_METHODS.update({
    "norm": {"pdf", "cdf", "quantile"}, "t": {"pdf", "cdf", "quantile"},
    "gamma": {"pdf", "cdf", "quantile"}, "chi2": {"pdf", "cdf", "quantile"},
    "f": {"pdf", "cdf", "quantile"}, "beta": {"pdf", "cdf", "quantile"},
    "nct": {"cdf"}, "ncx2": {"cdf"}, "ncf": {"cdf"},
    "binom": {"pmf", "cdf_int", "quantile_int"},
    "poisson": {"pmf", "cdf_int", "quantile_int"},
    "hypergeom": {"pmf", "cdf_int", "quantile_int"},
})


def field(dist, shape):
    return FIELD_RENAMES.get((dist, shape), shape.lower().rstrip("_"))


def num(x):
    """A MoonBit Double literal."""
    r = repr(float(x))
    if "e" in r and "." not in r.split("e")[0]:
        m, e = r.split("e")
        r = m + ".0e" + e
    elif "." not in r and "e" not in r and "inf" not in r and "nan" not in r:
        r += ".0"
    return r


def shape_infos(d):
    try:
        return list(d._shape_info())
    except Exception:
        return []


def box_checks(dist, infos):
    lines = []
    for s in infos:
        f = field(dist, s.name)
        lo, hi = float(s.domain[0]), float(s.domain[1])
        il, ih = s.inclusive
        conds, words = [], []
        if s.integrality:
            if math.isfinite(lo):
                conds.append(f"{f} >= {int(math.ceil(lo)) if il else int(math.floor(lo)) + 1}")
            if math.isfinite(hi):
                conds.append(f"{f} <= {int(math.floor(hi)) if ih else int(math.ceil(hi)) - 1}")
            value = f"{f}.to_double()"
        else:
            if math.isfinite(lo):
                conds.append(f"{f} {'>=' if il else '>'} {num(lo)}")
            if math.isfinite(hi):
                conds.append(f"{f} {'<=' if ih else '<'} {num(hi)}")
            if not conds:
                conds.append(f"!{f}.is_nan()")
            value = f
        if math.isfinite(lo) or math.isfinite(hi):
            words.append(
                ("[" if il else "(") + (repr(lo) if math.isfinite(lo) else "-inf")
                + ", " + (repr(hi) if math.isfinite(hi) else "inf") + ("]" if ih else ")")
            )
        else:
            words.append("not NaN")
        if conds:
            lines.append(
                f"  if !({' && '.join(conds)}) {{\n"
                f"    raise InvalidParam(name=\"{f}\", value={value}, reason=\"must be in {words[0]}\")\n"
                f"  }}"
            )
    return lines


def constraint(s):
    lo, hi = float(s.domain[0]), float(s.domain[1])
    # SciPy writes "strictly positive" as a lower bound of the smallest
    # subnormal in some families
    if 0.0 < lo < 1e-300:
        lo = 0.0
    if not math.isfinite(lo) and not math.isfinite(hi):
        return "Real"
    if lo == 0.0 and not math.isfinite(hi):
        return "Positive"
    if lo == 0.0 and hi == 1.0:
        return "UnitInterval"
    if math.isfinite(lo) and math.isfinite(hi):
        return f"Interval(lo={num(lo)}, hi={num(hi)})"
    if math.isfinite(lo):
        return f"Above({num(lo)})"
    return f"Below({num(hi)})"


def emit_family(name, kind):
    d = getattr(st, name)
    tname = NAMES[name]
    infos = shape_infos(d)
    fields = [(field(name, s.name), "Int" if s.integrality else "Double", s) for s in infos]
    builtin = BUILTIN.get(name, ())
    out = []
    if builtin:
        doc = [f"/// `scipy.stats.{name}` with {' and '.join(builtin)} as its own parameters."]
    else:
        doc = [
            f"/// `scipy.stats.{name}` in standard form"
            + (" (loc 0, scale 1); wrap in `Affine` for location and scale." if kind == "c"
               else " (loc 0)."),
        ]
    if name in EXTRA:
        doc.append(f"/// Parameters: {EXTRA[name][2]}.")
    out.append("///|\n" + "\n".join(doc))
    if fields or builtin:
        out.append(
            f"pub struct {tname} {{\n"
            + "\n".join([f"  priv {f} : {t}" for f, t, _ in fields] + [f"  priv {b} : Double" for b in builtin])
            + "\n}"
        )
    else:
        out.append(f"pub struct {tname} {{}}")
    # constructor
    args = ", ".join(
        [f"{f} : {t}" for f, t, _ in fields]
        + [f"{b}? : Double = {'0.0' if b == 'loc' else '1.0'}" for b in builtin]
    )
    out.append("///|")
    if fields or builtin:
        body = box_checks(name, infos)
        for b in builtin:
            cond = "!loc.is_nan() && !loc.is_inf()" if b == "loc" else "scale > 0.0 && !scale.is_inf()"
            why = "must be finite" if b == "loc" else "must be in (0, inf)"
            body.append(f"  if !({cond}) {{\n    raise InvalidParam(name=\"{b}\", value={b}, reason=\"{why}\")\n  }}")
        if name in EXTRA:
            cond, f0, why = EXTRA[name]
            t0 = dict((f, t) for f, t, _ in fields)[f0]
            v0 = f"{f0}.to_double()" if t0 == "Int" else f0
            body.append(
                f"  if !({cond}) {{\n    raise InvalidParam(name=\"{f0}\", value={v0}, reason=\"{why}\")\n  }}"
            )
        ctor = ", ".join([f for f, _, _ in fields] + list(builtin))
        if len(fields) + len(builtin) == 1:
            ctor += ","
        out.append(
            f"pub fn {tname}::new({args}) -> {tname} raise ParamError {{\n"
            + "\n".join(body) + f"\n  {{ {ctor} }}\n}}"
        )
    else:
        out.append(f"pub fn {tname}::new() -> {tname} {{\n  {tname}::{{  }}\n}}")
    # family: continuous parameters free, integer parameters fixed
    free = [(f, constraint(s)) for f, t, s in fields if t == "Double"]
    free += [(b, "Real" if b == "loc" else "Positive") for b in builtin]
    fixed = [(f, t) for f, t, _ in fields if t == "Int"]
    fixed_args = ", ".join(f"{f}~ : Int" for f, _ in fixed)
    specs = ", ".join(f"{{ name: \"{f}\", constraint: {c} }}" for f, c in free)
    make_args = []
    k = 0
    for f, t, _ in fields:
        if t == "Int":
            make_args.append(f)
        else:
            make_args.append(f"p[{k}]")
            k += 1
    for b in builtin:
        make_args.append(f"{b}=p[{k}]")
        k += 1
    make = f"{tname}::new({', '.join(make_args)})"
    params = "[" + ", ".join(f"d.{f}" for f, _ in free) + "]"
    out.append("///|")
    out.append(
        f"pub fn {tname}::family({fixed_args}) -> Family[{tname}] {{\n"
        f"  {{\n    name: \"{name}\",\n    specs: [{specs}],\n"
        f"    make: {'p' if free else '_'} => {make},\n"
        f"    params: {'d' if free else '_'} => {params if free else '[]'},\n  }}\n}}"
    )
    hand = HAND_METHODS.get(name, set())
    if kind == "c":
        for trait, sig in (
            ("Density", "pdf(self, x)"), ("Cdf", "cdf(self, x)"),
            ("Quantile", "quantile(self, p)"),
        ):
            if sig.split("(")[0] not in hand:
                a = sig[sig.index("(") + 1:-1]
                out.append(f"///|\npub impl {trait} for {tname} with {sig} {{\n  ignore(({a}))\n  ...\n}}")
        out.append(
            f"///|\n/// Inversion of the quantile function.\n"
            f"pub impl Sample for {tname} with sample(self, rng) {{\n"
            f"  Quantile::quantile(self, @num.uniform_open(rng))\n}}"
        )
        out.append(f"///|\npub impl Univariate for {tname}")
    else:
        for trait, sig in (
            ("Mass", "pmf(self, k)"), ("IntCdf", "cdf_int(self, k)"),
            ("IntQuantile", "quantile_int(self, p)"),
        ):
            if sig.split("(")[0] not in hand:
                a = sig[sig.index("(") + 1:-1]
                out.append(f"///|\npub impl {trait} for {tname} with {sig} {{\n  ignore(({a}))\n  ...\n}}")
        out.append(
            f"///|\n/// Inversion of the quantile function.\n"
            f"pub impl SampleInt for {tname} with sample_int(self, rng) {{\n"
            f"  IntQuantile::quantile_int(self, @num.uniform_open(rng))\n}}"
        )
        out.append(
            f"///|\npub impl Cdf for {tname} with cdf(self, x) {{\n"
            f"  IntCdf::cdf_int(self, @math.floor(x).to_int())\n}}"
        )
        out.append(
            f"///|\npub impl Cdf for {tname} with sf(self, x) {{\n"
            f"  IntCdf::sf_int(self, @math.floor(x).to_int())\n}}"
        )
        out.append(
            f"///|\npub impl Cdf for {tname} with sf_at_least(self, x) {{\n"
            f"  IntCdf::sf_int(self, @math.ceil(x).to_int() - 1)\n}}"
        )
        out.append(
            f"///|\npub impl Quantile for {tname} with quantile(self, p) {{\n"
            f"  IntQuantile::quantile_int(self, p).to_double()\n}}"
        )
        out.append(
            f"///|\npub impl Sample for {tname} with sample(self, rng) {{\n"
            f"  SampleInt::sample_int(self, rng).to_double()\n}}"
        )
        out.append(f"///|\npub impl Univariate for {tname}")
    return "\n\n".join(out)


def ctor_call(name, params, infos):
    tname = NAMES[name]
    vals = []
    for p, s in zip(params, infos):
        vals.append(str(int(round(p))) if s.integrality else num(p))
    return f"@dist.{tname}::new({', '.join(vals)})"


def argcheck_tests(name, kind, rng):
    d = getattr(st, name)
    infos = shape_infos(d)
    if not infos:
        return ""
    accept, reject = [], []
    for _ in range(400):
        ps = []
        for s in infos:
            lo, hi = max(float(s.domain[0]), -20.0), min(float(s.domain[1]), 20.0)
            v = rng.uniform(lo, hi)
            if s.integrality:
                v = float(round(v))
            ps.append(v)
        # stay inside the boxes, which the generated checks enforce exactly
        inside = all(
            (s.domain[0] <= v if s.inclusive[0] else s.domain[0] < v)
            and (v <= s.domain[1] if s.inclusive[1] else v < s.domain[1])
            for v, s in zip(ps, infos)
        )
        if not inside:
            continue
        ok = bool(np.all(d._argcheck(*ps)))
        (accept if ok else reject).append(ps)
    accept, reject = accept[:6], reject[:6]
    lines = [f"test \"catalog: {name} parameters agree with SciPy\" {{"]
    for ps in accept:
        lines.append(f"  ignore({ctor_call(name, ps, infos)})")
    for ps in reject:
        lines.append(
            f"  assert_true(\n    (try {ctor_call(name, ps, infos)} catch {{ _ => true }} noraise {{ _ => false }}),\n  )"
        )
    lines.append("}")
    return "///|\n" + "\n".join(lines)


def floats(xs):
    return "[" + ", ".join(num(x) for x in xs) + "]"


# The defining distribution functions, evaluated at high precision with
# mpmath. They are the oracle for the hand-written closed forms: densities
# come from differentiating these, quantiles from root finding on them, so
# the oracle shares no formula with the MoonBit code. Each is also checked
# against SciPy below, which catches a mismatch in parameter conventions.
def mp_cdfs():
    from mpmath import mp, mpf, exp, log, atan, asin, sin, tanh, pi, sqrt, power, ncdf, gammainc, betainc, binomial, inf
    one = mpf(1)

    def t_cdf(df):
        # near 0 use I_z(1/2, ν/2) with z = x²/(ν + x²) small; in the tails
        # I_t(ν/2, 1/2) with t = ν/(ν + x²) small, so neither form sits next
        # to an endpoint of betainc, where mpmath's derivative loses digits
        def F(x):
            if abs(x) < 1:
                h = betainc(one / 2, df / 2, 0, x * x / (df + x * x), regularized=True) / 2
                return one / 2 + h if x >= 0 else one / 2 - h
            lower = betainc(df / 2, one / 2, 0, df / (df + x * x), regularized=True) / 2
            return lower if x < 0 else 1 - lower
        return F

    def hyper(M, n, N):
        def pmf(k):
            return binomial(n, k) * binomial(M - n, N - k) / binomial(M, N)
        lo = max(0, N - (M - n))
        return lambda k: sum((pmf(j) for j in range(lo, int(k) + 1)), mpf(0)) if k >= lo else mpf(0)

    return {
        "norm": lambda: lambda x: ncdf(x),
        "t": t_cdf,
        "chi2": lambda df: lambda x: gammainc(df / 2, 0, x / 2, regularized=True),
        "gamma": lambda a: lambda x: gammainc(a, 0, x, regularized=True),
        "beta": lambda a, b: lambda x: betainc(a, b, 0, x, regularized=True),
        "f": lambda m, n: lambda x: betainc(m / 2, n / 2, 0, m * x / (m * x + n), regularized=True),
        "binom": lambda n, p: lambda k: mpf(0) if k < 0 else (one if k >= n else betainc(n - k, k + 1, 0, 1 - p, regularized=True)),
        "poisson": lambda mu: lambda k: mpf(0) if k < 0 else gammainc(k + 1, mu, inf, regularized=True),
        "hypergeom": hyper,
        "uniform": lambda: lambda x: x,
        "expon": lambda: lambda x: 1 - exp(-x),
        "laplace": lambda: lambda x: exp(x) / 2 if x < 0 else 1 - exp(-x) / 2,
        "logistic": lambda: lambda x: 1 / (1 + exp(-x)),
        "cauchy": lambda: lambda x: one / 2 + atan(x) / pi,
        "gumbel_r": lambda: lambda x: exp(-exp(-x)),
        "gumbel_l": lambda: lambda x: 1 - exp(-exp(x)),
        "weibull_min": lambda c: lambda x: 1 - exp(-power(x, c)),
        "weibull_max": lambda c: lambda x: exp(-power(-x, c)),
        "pareto": lambda b: lambda x: 1 - power(x, -b),
        "lomax": lambda c: lambda x: 1 - power(1 + x, -c),
        "rayleigh": lambda: lambda x: 1 - exp(-x * x / 2),
        "triang": lambda c: lambda x: x * x / c if x < c else 1 - (1 - x) ** 2 / (1 - c),
        "arcsine": lambda: lambda x: 2 / pi * asin(sqrt(x)),
        "halfcauchy": lambda: lambda x: 2 / pi * atan(x),
        "halflogistic": lambda: lambda x: tanh(x / 2),
        "hypsecant": lambda: lambda x: 2 / pi * atan(exp(x)),
        "genlogistic": lambda c: lambda x: power(1 + exp(-x), -c),
        "genpareto": lambda c: lambda x: 1 - power(1 + c * x, -1 / c),
        "genextreme": lambda c: lambda x: exp(-power(1 - c * x, 1 / c)),
        "burr": lambda c, d: lambda x: power(1 + power(x, -c), -d),
        "burr12": lambda c, d: lambda x: 1 - power(1 + power(x, c), -d),
        "fisk": lambda c: lambda x: 1 / (1 + power(x, -c)),
        "powerlaw": lambda a: lambda x: power(x, a),
        "bradford": lambda c: lambda x: log(1 + c * x) / log(1 + c),
        "loguniform": lambda a, b: lambda x: log(x / a) / log(b / a),
        "anglit": lambda: lambda x: sin(x + pi / 4) ** 2,
        "truncexpon": lambda b: lambda x: (1 - exp(-x)) / (1 - exp(-b)),
        "exponpow": lambda b: lambda x: 1 - exp(1 - exp(power(x, b))),
        "exponweib": lambda a, c: lambda x: power(1 - exp(-power(x, c)), a),
        "gompertz": lambda c: lambda x: 1 - exp(-c * (exp(x) - 1)),
        "kappa3": lambda a: lambda x: x * power(a + power(x, a), -1 / a),
        "dweibull": lambda c: lambda x: exp(-power(-x, c)) / 2 if x < 0 else 1 - exp(-power(x, c)) / 2,
        "truncpareto": lambda b, c: lambda x: (1 - power(x, -b)) / (1 - power(c, -b)),
        "loglaplace": lambda c: lambda x: power(x, c) / 2 if x < 1 else 1 - power(x, -c) / 2,
        "laplace_asymmetric": lambda k: lambda x: 1 - exp(-x * k) / (1 + k * k) if x >= 0 else k * k * exp(x / k) / (1 + k * k),
        "skewcauchy": lambda a: lambda x: (1 - a) / 2 + (1 - a) / pi * atan(x / (1 - a)) if x <= 0 else (1 - a) / 2 + (1 + a) / pi * atan(x / (1 + a)),
        # discrete: P(X <= k)
        "bernoulli": lambda p: lambda k: 0 if k < 0 else (1 - p if k < 1 else one),
        "geom": lambda p: lambda k: 0 if k < 1 else 1 - power(1 - p, k),
        "randint": lambda lo, hi: lambda k: 0 if k < lo else (one if k >= hi - 1 else (k - lo + 1) / mpf(hi - lo)),
        "dlaplace": lambda a: lambda k: 1 - exp(-a * (k + 1)) / (1 + exp(-a)) if k >= 0 else exp(a * k) / (1 + exp(-a)),
        "planck": lambda lam: lambda k: 0 if k < 0 else 1 - exp(-lam * (k + 1)),
        "boltzmann": lambda lam, n: lambda k: 0 if k < 0 else (one if k >= n - 1 else (1 - exp(-lam * (k + 1))) / (1 - exp(-lam * n))),
        "zipfian": lambda a, n: lambda k: 0 if k < 1 else (one if k >= n else sum(power(j, -a) for j in range(1, int(k) + 1)) / sum(power(j, -a) for j in range(1, int(n) + 1))),
    }


def mp_invert(F, target, x0, a, b):
    """The x in [a, b] with F(x) = target for nondecreasing F, by bisection
    from a bracket grown around the starting point."""
    from mpmath import mpf, fabs
    x0 = mpf(x0)
    step = (fabs(x0) + mpf(10) ** -12) * mpf(10) ** -9
    lo, hi = x0 - step, x0 + step
    lo = max(lo, mpf(a)) if math.isfinite(a) else lo
    hi = min(hi, mpf(b)) if math.isfinite(b) else hi
    for _ in range(400):
        if F(lo) <= target or (math.isfinite(a) and lo == a):
            break
        step *= 2
        lo = x0 - step
        if math.isfinite(a) and lo < a:
            lo = mpf(a)
    step = (fabs(x0) + mpf(10) ** -12) * mpf(10) ** -9
    for _ in range(400):
        if F(hi) >= target or (math.isfinite(b) and hi == b):
            break
        step *= 2
        hi = x0 + step
        if math.isfinite(b) and hi > b:
            hi = mpf(b)
    for _ in range(400):
        mid = (lo + hi) / 2
        if F(mid) < target:
            lo = mid
        else:
            hi = mid
        if hi - lo <= (fabs(mid) + mpf(10) ** -300) * mpf(10) ** -45:
            break
    return (lo + hi) / 2


def oracle_test(name, kind, params):
    from mpmath import mp, mpf, diff, log
    mp.dps = 60
    d = getattr(st, name)
    infos = shape_infos(d)
    frozen = d(*params)
    call = ctor_call(name, params, infos)
    mparams = [mpf(p) if not s.integrality else int(round(p)) for p, s in zip(params, infos)]
    F = mp_cdfs()[name](*mparams)
    ps = [1e-6, 0.01, 0.25, 0.5, 0.75, 0.99, 1 - 1e-6]

    def agree(ours, theirs, what):
        ours, theirs = float(ours), float(theirs)
        if abs(ours - theirs) > 1e-7 * abs(theirs) + 1e-12:
            raise SystemExit(f"{name}: oracle and SciPy disagree on {what}: {ours} vs {theirs}")

    if kind == "c":
        rows = []
        for p in ps:
            x = float(frozen.ppf(p))
            if not math.isfinite(x) or frozen.pdf(x) <= 0:
                continue
            X = mpf(x)
            c = F(X)
            dens = diff(F, X)
            rows.append((x, float(dens), float(log(dens)), float(c), float(1 - c)))
            agree(c, frozen.cdf(x), f"cdf({x})")
            agree(dens, frozen.pdf(x), f"pdf({x})")
        qs = []
        for p in ps:
            x0, y0 = float(frozen.ppf(p)), float(frozen.isf(p))
            if not (math.isfinite(x0) and math.isfinite(y0)):
                continue
            P = mpf(p)
            a, b = (float(v) for v in frozen.support())
            xq = mp_invert(F, P, x0, a, b)
            yq = mp_invert(F, 1 - P, y0, a, b)
            agree(xq, x0, f"ppf({p})")
            qs.append((p, float(xq), float(yq)))
        return (
            f"///|\ntest \"catalog: {name} matches the high-precision oracle\" {{\n"
            f"  let d = {call}\n"
            f"  check_continuous(\n    d,\n"
            f"    x={floats([r[0] for r in rows])},\n"
            f"    pdf={floats([r[1] for r in rows])},\n"
            f"    logpdf={floats([r[2] for r in rows])},\n"
            f"    cdf={floats([r[3] for r in rows])},\n"
            f"    sf={floats([r[4] for r in rows])},\n"
            f"    p={floats([q[0] for q in qs])},\n"
            f"    ppf={floats([q[1] for q in qs])},\n"
            f"    isf={floats([q[2] for q in qs])},\n  )\n}}"
        )
    ks = sorted({int(frozen.ppf(p)) for p in ps if math.isfinite(frozen.ppf(p))})
    rows = []
    for k in ks:
        c = F(k)
        m = c - F(k - 1)
        agree(c, frozen.cdf(k), f"cdf({k})")
        rows.append((k, float(m), float(c), float(1 - c)))
    qs = []
    for p in ps:
        k = int(frozen.ppf(p))
        P = mpf(p)
        while F(k - 1) >= P:
            k -= 1
        while F(k) < P:
            k += 1
        qs.append((p, k))
    return (
        f"///|\ntest \"catalog: {name} matches the high-precision oracle\" {{\n"
        f"  let d = {call}\n"
        f"  check_discrete(\n    d,\n"
        f"    k=[{', '.join(str(r[0]) for r in rows)}],\n"
        f"    pmf={floats([r[1] for r in rows])},\n"
        f"    cdf={floats([r[2] for r in rows])},\n"
        f"    sf={floats([r[3] for r in rows])},\n"
        f"    p={floats([q[0] for q in qs])},\n"
        f"    ppf=[{', '.join(str(q[1]) for q in qs)}],\n  )\n}}"
    )


HELPERS = '''///|
fn catalog_close(got : Double, want : Double, rel : Double) -> Bool {
  if want == got {
    return true
  }
  (got - want).abs() <= rel * want.abs() + 1.0e-15 * rel / 1.0e-13
}

///|
fn[D : @dist.Density + @dist.Cdf + @dist.Quantile] check_continuous(
  d : D,
  x~ : Array[Double],
  pdf~ : Array[Double],
  logpdf~ : Array[Double],
  cdf~ : Array[Double],
  sf~ : Array[Double],
  p~ : Array[Double],
  ppf~ : Array[Double],
  isf~ : Array[Double],
) -> Unit raise {
  for i, xi in x {
    assert_true(catalog_close(@dist.Density::pdf(d, xi), pdf[i], 1.0e-13), msg="pdf at \\{xi}: \\{@dist.Density::pdf(d, xi)} vs \\{pdf[i]}")
    assert_true(catalog_close(@dist.Density::logpdf(d, xi), logpdf[i], 1.0e-13), msg="logpdf at \\{xi}: \\{@dist.Density::logpdf(d, xi)} vs \\{logpdf[i]}")
    assert_true(catalog_close(@dist.Cdf::cdf(d, xi), cdf[i], 1.0e-13), msg="cdf at \\{xi}: \\{@dist.Cdf::cdf(d, xi)} vs \\{cdf[i]}")
    assert_true(catalog_close(@dist.Cdf::sf(d, xi), sf[i], 1.0e-13), msg="sf at \\{xi}: \\{@dist.Cdf::sf(d, xi)} vs \\{sf[i]}")
  }
  for i, pi in p {
    assert_true(catalog_close(@dist.Quantile::quantile(d, pi), ppf[i], 1.0e-13), msg="ppf at \\{pi}: \\{@dist.Quantile::quantile(d, pi)} vs \\{ppf[i]}")
    assert_true(catalog_close(@dist.Quantile::isf(d, pi), isf[i], 1.0e-13), msg="isf at \\{pi}: \\{@dist.Quantile::isf(d, pi)} vs \\{isf[i]}")
  }
}

///|
fn[D : @dist.Mass + @dist.IntCdf + @dist.IntQuantile] check_discrete(
  d : D,
  k~ : Array[Int],
  pmf~ : Array[Double],
  cdf~ : Array[Double],
  sf~ : Array[Double],
  p~ : Array[Double],
  ppf~ : Array[Int],
) -> Unit raise {
  for i, ki in k {
    assert_true(catalog_close(@dist.Mass::pmf(d, ki), pmf[i], 1.0e-13), msg="pmf at \\{ki}: \\{@dist.Mass::pmf(d, ki)} vs \\{pmf[i]}")
    assert_true(catalog_close(@dist.IntCdf::cdf_int(d, ki), cdf[i], 1.0e-13), msg="cdf at \\{ki}: \\{@dist.IntCdf::cdf_int(d, ki)} vs \\{cdf[i]}")
    assert_true(catalog_close(@dist.IntCdf::sf_int(d, ki), sf[i], 1.0e-13), msg="sf at \\{ki}: \\{@dist.IntCdf::sf_int(d, ki)} vs \\{sf[i]}")
  }
  for i, pi in p {
    assert_eq(@dist.IntQuantile::quantile_int(d, pi), ppf[i], msg="ppf at \\{pi}")
  }
}
'''


def main():
    rng = np.random.default_rng(20261119)
    header = (
        "// Generated by scripts/gen_catalog.py from SciPy "
        f"{scipy.__version__}; do not edit by hand.\n"
    )
    code, tests = [], []
    for kind, lst in (("c", distcont), ("d", distdiscrete)):
        seen = set()
        for name, params in lst:
            if not isinstance(name, str) or name in seen:
                continue
            seen.add(name)
            if name in EXISTING or name in ALIASES or name not in NAMES:
                continue
            code.append(emit_family(name, kind))
            tests.append(argcheck_tests(name, kind, rng))
            if name in mp_cdfs():
                tests.append(oracle_test(name, kind, params))
    text = header + "\n" + "\n\n".join(code) + "\n"
    (ROOT / "src/dist/catalog.mbt").write_text(text.replace("///|\n\n", "///|\n"))
    (ROOT / "src/dist/catalog_test.mbt").write_text(
        header + "\n" + HELPERS + "\n" + "\n\n".join(t for t in tests if t) + "\n"
    )
    print(f"{len(code)} families, {len(tests)} test blocks")


main()

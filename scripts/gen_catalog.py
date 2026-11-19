"""Generate the distribution catalogue from SciPy's metadata.

Writes src/dist/catalog.mbt and src/dist/catalog_test.mbt.

For every family the generator emits the type (SciPy's shape parameters,
plus `loc` and `scale` for a continuous one), its constructor, its
`Family` record, and the trait implementations. The formulas themselves
are written by hand, elsewhere in src/dist, as private methods the
generator finds by name and wraps:

- a continuous family's `T::std_pdf`, `std_cdf`, `std_sf`, `std_quantile`
  and `std_isf` (required) and `std_logpdf`, `std_logcdf`, `std_logsf`,
  `std_mean` with `std_variance`, `std_entropy`, `std_mgf`, `std_cf`,
  `std_support`, `std_mode` (optional) describe the standard form, at
  loc 0 and scale 1; location and scale are applied here, in one place,
  so no family can forget a tail;
- a discrete family's `T::lattice_pmf`, `lattice_cdf`, `lattice_sf`,
  `lattice_quantile` and `lattice_support` take the point as an
  integer-valued Double, and both the `Double` and the `Int64` interface
  are built on them.

A required method not written yet becomes a typed `...` stub. Traits a
family implements directly by hand (`pub impl Moments for Binomial`) are
found too and exposed with `pub extend` like the generated ones.

The tests check the constructors against SciPy's _argcheck and the
hand-written formulas against a high-precision oracle.

Run from the repository root:
    uv run --no-project --with scipy --with mpmath python -I scripts/gen_catalog.py
"""

import math
import re
import pathlib

import numpy as np
import scipy
import scipy.stats as st
from scipy.stats._distr_params import distcont, distdiscrete
from scipy.stats._distn_infrastructure import rv_generic

ROOT = pathlib.Path(__file__).resolve().parent.parent

# SciPy names not generated.
EXISTING = set()
# Extra lines for a family's doc comment, where SciPy's convention is easy
# to get wrong.
NOTES = {
    "randint": ["The integers `low ..< high`: `high` itself is not included."],
    "hypergeom": [
        "`population` items of which `successes` are marked; `draws` are taken",
        "without replacement, and the variable counts the marked ones drawn.",
    ],
}
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

def hand_written():
    """What is written by hand anywhere in src/dist, so that the generator
    wraps it and never emits a second copy: the private standard-form and
    lattice methods as (type, method), and the trait implementations as
    {(type, trait): methods}."""
    methods, impls = set(), {}
    for path in sorted((ROOT / "src/dist").glob("*.mbt")):
        if path.name.startswith("catalog") or path.name.endswith("_test.mbt"):
            continue
        text = path.read_text()
        for m in re.finditer(r"^fn (\w+)::((?:std|lattice)_\w+)\(", text, re.M):
            methods.add((m.group(1), m.group(2)))
        for m in re.finditer(r"^pub impl (\w+) for (\w+) with fn (\w+)\(", text, re.M):
            impls.setdefault((m.group(2), m.group(1)), set()).add(m.group(3))
    return methods, impls


HAND, HAND_IMPLS = hand_written()

# The methods `pub extend` exposes for each trait, so that callers write
# `d.cdf(x)`: an `impl` is not reachable with dot syntax otherwise.
TRAIT_METHODS = {
    "Density": ["pdf", "logpdf"],
    "Cdf": ["cdf", "sf", "logcdf", "logsf", "sf_at_least"],
    "Quantile": ["quantile", "isf"],
    "Sample": ["sample", "sample_into"],
    "Mass": ["pmf", "logpmf"],
    "IntCdf": ["cdf_int", "sf_int"],
    "IntQuantile": ["quantile_int"],
    "SampleInt": ["sample_int"],
    "Moments": ["mean", "variance"],
    "Shape": ["skewness", "excess_kurtosis"],
    "Entropy": ["entropy"],
    "Mgf": ["mgf", "cgf"],
    "Characteristic": ["cf"],
    "Support": ["support_lower", "support_upper"],
    "Mode": ["mode"],
    "Fit": ["fit"],
}
EXTENDED = {m for ms in TRAIT_METHODS.values() for m in ms}


def field(dist, shape):
    return FIELD_RENAMES.get((dist, shape), shape.lower().rstrip("_"))


def num(x):
    """A MoonBit Double literal."""
    x = float(x)
    if math.isinf(x):
        return "1.0 / 0.0" if x > 0 else "-1.0 / 0.0"
    r = repr(x)
    if "e" in r and "." not in r.split("e")[0]:
        m, e = r.split("e")
        r = m + ".0e" + e
    elif "." not in r and "e" not in r and "nan" not in r:
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


def constant_support(d):
    """SciPy's support when it does not depend on the shape parameters."""
    if type(d)._get_support is not rv_generic._get_support:
        return None
    return float(d.a), float(d.b)


def stub(tname, method, params, ret):
    names = ", ".join(p.split(" : ")[0] for p in ["self"] + params)
    args = ", ".join([f"self : {tname}"] + params)
    use = f"ignore(({names}))" if params else "ignore(self)"
    return f"///|\nfn {tname}::{method}({args}) -> {ret} {{\n  {use}\n  ...\n}}"


def emit_family(name, kind):
    d = getattr(st, name)
    tname = NAMES[name]
    infos = shape_infos(d)
    fields = [(field(name, s.name), "Int64" if s.integrality else "Double", s) for s in infos]
    clash = {f for f, _, _ in fields} & (EXTENDED | {"loc", "scale"})
    if clash:
        raise SystemExit(f"{name}: parameter names {clash} collide with method names")
    hand = {m for t, m in HAND if t == tname}
    labelled = len(fields) >= 2
    out = []
    if kind == "c":
        doc = [
            f"/// `scipy.stats.{name}`, with SciPy's shape parameters, location `loc`",
            "/// and scale `scale`: the standard form at (x - loc) / scale.",
        ]
    else:
        doc = [f"/// `scipy.stats.{name}` (loc 0)."]
    if name in EXTRA:
        doc.append(f"/// Parameters: {EXTRA[name][2]}.")
    if name in NOTES:
        doc.extend(f"/// {line}" for line in NOTES[name])
    own = [f"  {f} : {t}" for f, t, _ in fields]
    if kind == "c":
        own += ["  loc : Double", "  scale : Double"]
    out.append("///|\n" + "\n".join(doc) + "\n#derive(Eq, Debug)")
    out[-1] += f"\npub struct {tname} {{\n" + "\n".join(own) + "\n}" if own else f"\npub struct {tname} {{}}"
    # constructor: one shape parameter is positional, several are labelled,
    # so that two parameters of the same type cannot be swapped silently
    args = [f"{f}{'~' if labelled else ''} : {t}" for f, t, _ in fields]
    body = box_checks(name, infos)
    if kind == "c":
        args += ["loc~ : Double = 0.0", "scale~ : Double = 1.0"]
        body.append(
            "  if !(!loc.is_nan() && !loc.is_inf()) {\n"
            "    raise InvalidParam(name=\"loc\", value=loc, reason=\"must be finite\")\n  }"
        )
        body.append(
            "  if !(scale > 0.0 && !scale.is_inf()) {\n"
            "    raise InvalidParam(name=\"scale\", value=scale, reason=\"must be in (0, inf)\")\n  }"
        )
    if name in EXTRA:
        cond, f0, why = EXTRA[name]
        t0 = dict((f, t) for f, t, _ in fields)[f0]
        v0 = f"{f0}.to_double()" if t0 == "Int64" else f0
        body.append(
            f"  if !({cond}) {{\n    raise InvalidParam(name=\"{f0}\", value={v0}, reason=\"{why}\")\n  }}"
        )
    members = [f for f, _, _ in fields] + (["loc", "scale"] if kind == "c" else [])
    out.append("///|")
    if kind == "c":
        out.append(
            "/// The distribution with these shapes, location `loc` and scale\n"
            "/// `scale`; raises `ParamError` outside SciPy's parameter domain."
        )
    elif members:
        out.append(
            "/// The distribution with these parameters; raises `ParamError`\n"
            "/// outside SciPy's parameter domain."
        )
    else:
        out.append("/// The distribution, which has no parameters.")
    if members:
        out.append(
            f"pub fn {tname}::new({', '.join(args)}) -> {tname} raise ParamError {{\n"
            + "\n".join(body) + f"\n  {{ {', '.join(members)}{',' if len(members) == 1 else ''} }}\n}}"
        )
    else:
        out.append(f"pub fn {tname}::new() -> {tname} {{\n  {tname}::{{  }}\n}}")
    # family: continuous parameters free, integer parameters fixed
    free = [(f, constraint(s)) for f, t, s in fields if t == "Double"]
    if kind == "c":
        free += [("loc", "Real"), ("scale", "Positive")]
    fixed = [(f, t) for f, t, _ in fields if t == "Int64"]
    fixed_args = ", ".join(f"{f}~ : Int64" for f, _ in fixed)
    specs = ", ".join(f"{{ name: \"{f}\", constraint: {c} }}" for f, c in free)
    make_args = []
    k = 0
    for f, t, _ in fields:
        value = f if t == "Int64" else f"p[{k}]"
        if t == "Double":
            k += 1
        make_args.append(f"{f}={value}" if labelled else value)
    if kind == "c":
        make_args += [f"loc=p[{k}]", f"scale=p[{k + 1}]"]
    make = f"{tname}::new({', '.join(make_args)})"
    params = "[" + ", ".join(f"d.{f}" for f, _ in free) + "]"
    out.append("///|")
    if fixed:
        out.append(
            "/// The family for fitting, over the continuous parameters, with the\n"
            "/// integer ones fixed at the values given."
        )
    else:
        out.append("/// The family for fitting, over every parameter.")
    out.append(
        f"pub fn {tname}::family({fixed_args}) -> Family[{tname}] {{\n"
        f"  {{\n    name: \"{name}\",\n    specs: [{specs}],\n"
        f"    make: {'p' if free else '_'} => {make},\n"
        f"    params: {'d' if free else '_'} => {params if free else '[]'},\n  }}\n}}"
    )
    traits = {}

    def impl(trait, method, sig, body):
        traits.setdefault(trait, []).append(method)
        out.append(f"///|\npub impl {trait} for {tname} with fn {method}({sig}) {{\n  {body}\n}}")

    if kind == "c":
        # standard forms not written by hand yet are typed stubs; `std_sf`
        # and `std_isf` are never filled in from the other tail
        for method, param in (
            ("std_pdf", "x"), ("std_cdf", "x"), ("std_sf", "x"),
            ("std_quantile", "p"), ("std_isf", "q"),
        ):
            if method not in hand:
                out.append(stub(tname, method, [f"{param} : Double"], "Double"))
        support = constant_support(d)
        if "std_support" not in hand and support is not None:
            out.append(
                f"///|\nfn {tname}::std_support(_self : {tname}) -> (Double, Double) {{\n"
                f"  ({num(support[0])}, {num(support[1])})\n}}"
            )
            hand.add("std_support")
        z = "(x - self.loc) / self.scale"
        impl("Density", "pdf", "self, x", f"self.std_pdf({z}) / self.scale")
        if "std_logpdf" in hand:
            impl("Density", "logpdf", "self, x", f"self.std_logpdf({z}) - @math.ln(self.scale)")
        else:
            impl("Density", "logpdf", "self, x", f"@math.ln(self.std_pdf({z})) - @math.ln(self.scale)")
        impl("Cdf", "cdf", "self, x", f"self.std_cdf({z})")
        impl("Cdf", "sf", "self, x", f"self.std_sf({z})")
        for m in ("logcdf", "logsf"):
            if f"std_{m}" in hand:
                impl("Cdf", m, "self, x", f"self.std_{m}({z})")
        impl("Quantile", "quantile", "self, p", "self.loc + self.scale * self.std_quantile(p)")
        impl("Quantile", "isf", "self, q", "self.loc + self.scale * self.std_isf(q)")
        out.append("///|\n/// Inversion of the quantile function.")
        out[-1] = out[-1]  # the doc line attaches to the impl below
        traits.setdefault("Sample", []).append("sample")
        out[-1] += (
            f"\npub impl Sample for {tname} with fn sample(self, rng) {{\n"
            f"  Quantile::quantile(self, @num.uniform_open(rng))\n}}"
        )
        if "std_mean" in hand and "std_variance" in hand:
            impl(
                "Moments", "mean", "self",
                "match self.std_mean() {\n    Some(m) => Some(self.loc + self.scale * m)\n    None => None\n  }",
            )
            impl(
                "Moments", "variance", "self",
                "match self.std_variance() {\n    Some(v) => Some(self.scale * self.scale * v)\n    None => None\n  }",
            )
        if "std_entropy" in hand:
            impl("Entropy", "entropy", "self", "self.std_entropy() + @math.ln(self.scale)")
        if "std_mgf" in hand:
            impl(
                "Mgf", "mgf", "self, t",
                "match self.std_mgf(self.scale * t) {\n    Some(m) => Some(@math.exp(self.loc * t) * m)\n    None => None\n  }",
            )
        if "std_cf" in hand:
            impl(
                "Characteristic", "cf", "self, t",
                "let (re, im) = self.std_cf(self.scale * t)\n"
                "  let c = @math.cos(self.loc * t)\n  let s = @math.sin(self.loc * t)\n"
                "  (re * c - im * s, re * s + im * c)",
            )
        if "std_support" in hand:
            impl("Support", "support_lower", "self", "self.loc + self.scale * self.std_support().0")
            impl("Support", "support_upper", "self", "self.loc + self.scale * self.std_support().1")
        if "std_mode" in hand:
            impl("Mode", "mode", "self", "self.loc + self.scale * self.std_mode()")
    else:
        for method, params, ret in (
            ("lattice_pmf", ["k : Double"], "Double"),
            ("lattice_cdf", ["k : Double"], "Double"),
            ("lattice_sf", ["k : Double"], "Double"),
            ("lattice_quantile", ["q : Double"], "Double"),
        ):
            if method not in hand:
                out.append(stub(tname, method, params, ret))
        support = constant_support(d)
        if "lattice_support" not in hand:
            if support is None:
                out.append(stub(tname, "lattice_support", [], "(Double, Double)"))
            else:
                out.append(
                    f"///|\nfn {tname}::lattice_support(_self : {tname}) -> (Double, Double) {{\n"
                    f"  ({num(support[0])}, {num(support[1])})\n}}"
                )
        impl("Mass", "pmf", "self, k", "self.lattice_pmf(k.to_double())")
        impl("IntCdf", "cdf_int", "self, k", "self.lattice_cdf(k.to_double())")
        impl("IntCdf", "sf_int", "self, k", "self.lattice_sf(k.to_double())")
        impl("IntQuantile", "quantile_int", "self, p", "saturate_int64(self.lattice_quantile(p))")
        impl(
            "SampleInt", "sample_int", "self, rng",
            "saturate_int64(self.lattice_quantile(@num.uniform_open(rng)))",
        )
        # the Double interface takes floor(x) straight to the lattice
        # methods, without an integer type in between
        impl("Cdf", "cdf", "self, x", "self.lattice_cdf(@math.floor(x))")
        impl("Cdf", "sf", "self, x", "self.lattice_sf(@math.floor(x))")
        impl("Cdf", "sf_at_least", "self, x", "self.lattice_sf(@math.ceil(x) - 1.0)")
        impl("Quantile", "quantile", "self, p", "self.lattice_quantile(p)")
        impl(
            "Quantile", "isf", "self, q",
            "let (lo, hi) = self.lattice_support()\n  lattice_isf(k => self.lattice_sf(k), q, lo, hi)",
        )
        impl("Sample", "sample", "self, rng", "self.lattice_quantile(@num.uniform_open(rng))")
        impl("Support", "support_lower", "self", "self.lattice_support().0")
        impl("Support", "support_upper", "self", "self.lattice_support().1")
    out.append(f"///|\npub impl Univariate for {tname}")
    # dot syntax for everything the type implements, generated or by hand
    for (t, trait), methods in sorted(HAND_IMPLS.items()):
        if t == tname:
            traits.setdefault(trait, [])
    for trait in sorted(traits):
        out.append(
            f"///|\npub extend {tname} with {trait}::{{{', '.join(TRAIT_METHODS[trait])}}}"
        )
    return "\n\n".join(out)


def ctor_call(name, params, infos, extra=()):
    tname = NAMES[name]
    vals = []
    for p, s in zip(params, infos):
        v = str(int(round(p))) if s.integrality else num(p)
        vals.append(f"{field(name, s.name)}={v}" if len(infos) >= 2 else v)
    return f"@dist.{tname}::new({', '.join(vals + list(extra))})"


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
    tname = NAMES[name]
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
            rows.append((x, dens, c))
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
            qs.append((p, xq, yq))
        # the log tails where the tail itself is at most a half: a default
        # ln(tail) is accurate there, and a closed form must agree with it.
        # Next to a finite end of the support the tail is ill-conditioned in
        # x (one ulp of x moves ln(tail) by far more than the tolerance), so
        # a point is used only where the condition number
        # |x f(x) / (tail ln(tail))| is small.
        def well_conditioned(x, dens, tail):
            return abs(x * dens / (tail * log(tail))) < 100

        lower = [(x, log(c)) for x, f, c in rows if c <= 0.5 and well_conditioned(x, f, c)]
        upper = [(x, log(1 - c)) for x, f, c in rows if 1 - c <= 0.5 and well_conditioned(x, f, 1 - c)]
        # far beyond where the tail underflows, for the families that state
        # their log tails in closed form
        hand = {m for t, m in HAND if t == tname}
        mp.dps = 1200
        if "std_logsf" in hand:
            x = float(frozen.isf(1e-300))
            if math.isfinite(x):
                x += 0.3 * abs(x)
                upper.append((x, log(1 - F(mpf(x)))))
        if "std_logcdf" in hand:
            x = float(frozen.ppf(1e-300))
            if math.isfinite(x):
                x -= 0.3 * abs(x)
                lower.append((x, log(F(mpf(x)))))
        mp.dps = 60

        def block(title, ctor, k):
            # k is the scale: x and the quantiles are multiplied by it, the
            # density divided; a power of two keeps every value exact
            return (
                f"///|\ntest \"catalog: {name} {title}\" {{\n"
                f"  let d = {ctor}\n"
                f"  check_continuous(\n    d,\n"
                f"    x={floats([k * r[0] for r in rows])},\n"
                f"    pdf={floats([r[1] / k for r in rows])},\n"
                f"    logpdf={floats([log(r[1]) - log(mpf(k)) for r in rows])},\n"
                f"    cdf={floats([r[2] for r in rows])},\n"
                f"    sf={floats([1 - r[2] for r in rows])},\n"
                f"    p={floats([q[0] for q in qs])},\n"
                f"    ppf={floats([k * q[1] for q in qs])},\n"
                f"    isf={floats([k * q[2] for q in qs])},\n"
                f"    lower_x={floats([k * r[0] for r in lower])},\n"
                f"    logcdf={floats([r[1] for r in lower])},\n"
                f"    upper_x={floats([k * r[0] for r in upper])},\n"
                f"    logsf={floats([r[1] for r in upper])},\n  )\n}}"
            )

        return (
            block("matches the high-precision oracle", call, 1)
            + "\n\n"
            + block(
                "keeps its accuracy under a scale",
                ctor_call(name, params, infos, extra=("scale=2.0",)), 2,
            )
        )
    ks = sorted({int(frozen.ppf(p)) for p in ps if math.isfinite(frozen.ppf(p))})
    rows = []
    for k in ks:
        c = F(k)
        m = c - F(k - 1)
        agree(c, frozen.cdf(k), f"cdf({k})")
        rows.append((k, float(m), float(c), float(1 - c)))
    lo, hi = (float(v) for v in frozen.support())
    qs = []
    for p in ps:
        k = int(frozen.ppf(p))
        P = mpf(p)
        while F(k - 1) >= P:
            k -= 1
        while F(k) < P:
            k += 1
        # the smallest j with P(X > j) <= p
        j = int(frozen.isf(p))
        while j - 1 >= lo and 1 - F(j - 1) <= P:
            j -= 1
        while 1 - F(j) > P:
            j += 1
        qs.append((p, k, j))
    return (
        f"///|\ntest \"catalog: {name} matches the high-precision oracle\" {{\n"
        f"  let d = {call}\n"
        f"  check_discrete(\n    d,\n"
        f"    k=[{', '.join(str(r[0]) for r in rows)}],\n"
        f"    pmf={floats([r[1] for r in rows])},\n"
        f"    cdf={floats([r[2] for r in rows])},\n"
        f"    sf={floats([r[3] for r in rows])},\n"
        f"    p={floats([q[0] for q in qs])},\n"
        f"    ppf=[{', '.join(str(q[1]) for q in qs)}],\n"
        f"    isf=[{', '.join(str(q[2]) for q in qs)}],\n  )\n}}"
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
  lower_x~ : Array[Double],
  logcdf~ : Array[Double],
  upper_x~ : Array[Double],
  logsf~ : Array[Double],
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
  for i, xi in lower_x {
    assert_true(catalog_close(@dist.Cdf::logcdf(d, xi), logcdf[i], 1.0e-13), msg="logcdf at \\{xi}: \\{@dist.Cdf::logcdf(d, xi)} vs \\{logcdf[i]}")
  }
  for i, xi in upper_x {
    assert_true(catalog_close(@dist.Cdf::logsf(d, xi), logsf[i], 1.0e-13), msg="logsf at \\{xi}: \\{@dist.Cdf::logsf(d, xi)} vs \\{logsf[i]}")
  }
}

///|
/// Both interfaces of a distribution on the integers: the `Int64` one and
/// the `Double` one, which must agree at every integer.
fn[D : @dist.Mass + @dist.IntCdf + @dist.IntQuantile + @dist.Cdf + @dist.Quantile] check_discrete(
  d : D,
  k~ : Array[Int64],
  pmf~ : Array[Double],
  cdf~ : Array[Double],
  sf~ : Array[Double],
  p~ : Array[Double],
  ppf~ : Array[Int64],
  isf~ : Array[Int64],
) -> Unit raise {
  for i, ki in k {
    assert_true(catalog_close(@dist.Mass::pmf(d, ki), pmf[i], 1.0e-13), msg="pmf at \\{ki}: \\{@dist.Mass::pmf(d, ki)} vs \\{pmf[i]}")
    assert_true(catalog_close(@dist.IntCdf::cdf_int(d, ki), cdf[i], 1.0e-13), msg="cdf at \\{ki}: \\{@dist.IntCdf::cdf_int(d, ki)} vs \\{cdf[i]}")
    assert_true(catalog_close(@dist.IntCdf::sf_int(d, ki), sf[i], 1.0e-13), msg="sf at \\{ki}: \\{@dist.IntCdf::sf_int(d, ki)} vs \\{sf[i]}")
    let x = ki.to_double()
    assert_eq(@dist.Cdf::cdf(d, x + 0.5), @dist.IntCdf::cdf_int(d, ki))
    assert_eq(@dist.Cdf::sf(d, x + 0.5), @dist.IntCdf::sf_int(d, ki))
    assert_eq(@dist.Cdf::sf_at_least(d, x), @dist.IntCdf::sf_int(d, ki - 1))
  }
  for i, pi in p {
    assert_eq(@dist.IntQuantile::quantile_int(d, pi), ppf[i], msg="ppf at \\{pi}")
    assert_eq(@dist.Quantile::quantile(d, pi), ppf[i].to_double(), msg="quantile at \\{pi}")
    assert_eq(@dist.Quantile::isf(d, pi), isf[i].to_double(), msg="isf at \\{pi}")
  }
}
'''


def main():
    rng = np.random.default_rng(20261119)
    header = (
        "// Generated by scripts/gen_catalog.py from SciPy "
        f"{scipy.__version__}; DO NOT EDIT.\n"

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

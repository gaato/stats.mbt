"""Generate src/special/special_test.mbt: special functions against mpmath.

Points are drawn from each function's domain with extra weight on the hard
places (tails, the zeros of ln Γ and ψ, large parameters, p near 0 and 1);
the expected values are computed by mpmath at 50 digits, inverses by
bisection on the forward function in mpmath.

    uv run --no-project --with mpmath python -I scripts/gen_special_tests.py [--report]

--report sets every tolerance to 0, so each test fails and prints the
largest relative error it saw; that is how the tolerances below were set.
"""

import math
import pathlib
import random
import sys

from mpmath import mp, mpf, erf, erfc, exp, log, ncdf, loggamma, digamma, psi, gammainc, betainc, quad, pi, nsum, inf, binomial

mp.dps = 50
ROOT = pathlib.Path(__file__).resolve().parent.parent
REPORT = "--report" in sys.argv
rnd = random.Random(20261119)


def lit(x):
    r = repr(float(x))
    if r in ("inf", "-inf", "nan"):
        raise ValueError(r)
    if "e" in r and "." not in r.split("e")[0]:
        m, e = r.split("e")
        r = m + ".0e" + e
    elif "." not in r and "e" not in r:
        r += ".0"
    return r


def logu(lo, hi):
    return math.exp(rnd.uniform(math.log(lo), math.log(hi)))


def bisect(f, target, lo, hi, iters=400):
    lo, hi = mpf(lo), mpf(hi)
    for _ in range(iters):
        mid = (lo + hi) / 2
        if f(mid) < target:
            lo = mid
        else:
            hi = mid
    return (lo + hi) / 2


def ln_bisect(f, target, lo, hi, iters=400):
    """For targets far below 1: bisect on log f against log target."""
    lt = log(target)
    return bisect(lambda x: log(f(x)) if f(x) > 0 else mpf(-10) ** 9, lt, lo, hi, iters)


def logx_bisect(f, target, lo=mpf("1e-330"), hi=mpf(1e6), iters=500, log_target=False):
    """For increasing f on positive x: bisect in u = ln x, so solutions at
    any scale (1e-300 as well as 1e5) are resolved to ~1e-100 relative."""
    lo, hi = log(mpf(lo)), log(mpf(hi))
    lt = log(target) if log_target else target
    for _ in range(iters):
        mid = (lo + hi) / 2
        v = f(exp(mid))
        v = (log(v) if v > 0 else mpf(-10) ** 9) if log_target else v
        if v < lt:
            lo = mid
        else:
            hi = mid
    return exp((lo + hi) / 2)


# Tolerances: about twice the worst relative error measured with --report.
# The incomplete gamma and beta functions and their inverses get 1e-13 to
# 2e-13 because their worst cases are limited by conditioning, not by the
# method: near x ~ 1e-287, say, d ln x / d ln a is |ln x| ~ 660, so the
# rounding of a alone is 7e-14 (computed per case with mpmath when these
# were set).
TOL = {
    "erf": 1e-15, "erfc": 5e-15, "erfcx": 3e-15, "ndtr": 5e-15,
    "log_ndtr": 3e-15, "ndtri": 1e-15, "erfinv": 1e-15, "lgamma": 1e-15,
    "digamma": 5e-15, "trigamma": 1e-15, "lbeta": 2e-15, "lchoose": 2e-15,
    "gammainc": 1e-13, "gammaincc": 2e-13, "gammaincinv": 1e-15,
    "gammainccinv": 1e-13, "betainc": 2e-13, "betaincc": 2e-13,
    "betaincinv": 2e-13, "owens_t": 3e-15, "kolmogorov_sf": 1e-14,
}

TESTS = []


def emit(name, call, args, expected, tol, floor=0.0):
    """call: MoonBit expression over the argument names a0, a1, ..."""
    rows = [(a, e) for a, e in zip(args, expected)]
    n = len(args[0]) if args and isinstance(args[0], tuple) else 1
    cols = list(zip(*[a if isinstance(a, tuple) else (a,) for a, _ in rows]))
    lines = [f"///|\ntest \"special: {name} against mpmath\" {{"]
    for i, col in enumerate(cols):
        lines.append(f"  let a{i} : Array[Double] = [{', '.join(lit(v) for v in col)}]")
    lines.append(f"  let want : Array[Double] = [{', '.join(lit(e) for _, e in rows)}]")
    lines.append("  let mut worst = 0.0")
    lines.append("  let mut at = 0")
    lines.append("  for i in 0..<want.length() {")
    lines.append(f"    let got = {call}")
    lines.append(f"    let err = (got - want[i]).abs() / (want[i].abs() + {lit(floor)})")
    lines.append("    if !(err <= worst) {")
    lines.append("      worst = err")
    lines.append("      at = i")
    lines.append("    }")
    lines.append("  }")
    t = 0.0 if REPORT else TOL[name]
    lines.append(f"  assert_true(worst <= {lit(t)}, msg=\"{name}: worst relative error \\{{worst}} at case \\{{at}}\")")
    lines.append("}")
    TESTS.append("\n".join(lines))


def dbl(x):
    return float(x)


# --- error functions -----------------------------------------------------
xs = [rnd.uniform(-6, 6) for _ in range(60)] + [logu(1e-300, 0.5) * rnd.choice([-1, 1]) for _ in range(20)]
emit("erf", "@special.erf(a0[i])", xs, [erf(mpf(x)) for x in xs], 2e-16)
xs = [rnd.uniform(-6, 27) for _ in range(80)] + [rnd.uniform(0.9, 1.1) for _ in range(10)]
emit("erfc", "@special.erfc(a0[i])", xs, [erfc(mpf(x)) for x in xs], 2e-16)
xs = [logu(1e-3, 1e7) for _ in range(50)]
emit("erfcx", "@special.erfcx(a0[i])", xs, [erfc(mpf(x)) * exp(mpf(x) ** 2) for x in xs], 2e-16)
zs = [rnd.uniform(-38, 9) for _ in range(80)]
emit("ndtr", "@special.ndtr(a0[i])", zs, [ncdf(mpf(z)) for z in zs], 2e-16)
zs = [-logu(1e-3, 1e5) for _ in range(40)] + [rnd.uniform(-3, 8) for _ in range(30)]
emit("log_ndtr", "@special.log_ndtr(a0[i])", zs, [log(ncdf(mpf(z))) for z in zs], 2e-16)
ps = [logu(1e-300, 0.5) for _ in range(50)] + [1 - logu(1e-16, 0.5) for _ in range(30)]
emit(
    "ndtri", "@special.ndtri(a0[i])", ps,
    [ln_bisect(ncdf, mpf(p), -40, 40) if p < 0.5 else bisect(lambda z: -(1 - ncdf(z)), -(1 - mpf(p)), -40, 40) for p in ps],
    2e-16,
)
ys = [rnd.uniform(-0.99, 0.99) for _ in range(30)] + [logu(1e-300, 0.1) for _ in range(15)] + [1 - logu(1e-16, 0.1) for _ in range(20)]
emit(
    "erfinv", "@special.erfinv(a0[i])", ys,
    [math.copysign(1, y) * (logx_bisect(erf, abs(mpf(y)), hi=10, log_target=True) if abs(y) < 0.5 else logx_bisect(lambda x: -erfc(x), -(1 - abs(mpf(y))), hi=30)) for y in ys],
    2e-16,
)

# --- gamma family ----------------------------------------------------------
xs = (
    [logu(1e-300, 1e10) for _ in range(60)]
    + [1 + rnd.choice([-1, 1]) * logu(1e-15, 0.4) for _ in range(15)]
    + [2 + rnd.choice([-1, 1]) * logu(1e-15, 0.4) for _ in range(15)]
    + [-rnd.uniform(0.01, 50) for _ in range(20)]
)
xs = [x for x in xs if x != round(x)]
emit("lgamma", "@special.lgamma(a0[i])", xs, [mp.re(loggamma(mpf(x))) if x > 0 else log(abs(mp.gamma(mpf(x)))) for x in xs], 2e-16)
root = mp.findroot(digamma, 1.46)
xs = [logu(1e-10, 1e10) for _ in range(60)] + [float(root) + rnd.uniform(-0.3, 0.3) for _ in range(20)] + [-rnd.uniform(0.01, 30) for _ in range(15)]
emit("digamma", "@special.digamma(a0[i])", xs, [digamma(mpf(x)) for x in xs], 2e-16)
xs = [logu(1e-10, 1e10) for _ in range(60)] + [-rnd.uniform(0.01, 30) for _ in range(15)]
emit("trigamma", "@special.trigamma(a0[i])", xs, [psi(1, mpf(x)) for x in xs], 2e-16)
pairs = [(logu(1e-3, 1e7), logu(1e-3, 1e7)) for _ in range(80)]
emit("lbeta", "@special.lbeta(a0[i], a1[i])", pairs, [loggamma(mpf(a)) + loggamma(mpf(b)) - loggamma(mpf(a) + mpf(b)) for a, b in pairs], 2e-16)
pairs = []
for _ in range(60):
    n = float(round(logu(2, 1e7)))
    k = float(rnd.randint(1, int(n) - 1))
    pairs.append((n, k))
emit("lchoose", "@special.lchoose(a0[i], a1[i])", pairs, [log(binomial(mpf(n), mpf(k))) for n, k in pairs], 2e-16)

# --- incomplete gamma -------------------------------------------------------
cases = []
for _ in range(120):
    a = logu(1e-3, 1e4)
    x = a * logu(1e-2, 1e2) if rnd.random() < 0.7 else a + rnd.uniform(-3, 3) * math.sqrt(a)
    if x > 0:
        cases.append((a, x))
emit("gammainc", "@special.gammainc(a0[i], a1[i])", cases, [gammainc(mpf(a), 0, mpf(x), regularized=True) for a, x in cases], 2e-16, floor=1e-300)
emit("gammaincc", "@special.gammaincc(a0[i], a1[i])", cases, [gammainc(mpf(a), mpf(x), inf, regularized=True) for a, x in cases], 2e-16, floor=1e-300)
cases = []
for _ in range(60):
    a = logu(1e-2, 1e4)
    p = logu(1e-200, 0.5) if rnd.random() < 0.6 else 1 - logu(1e-14, 0.5)
    cases.append((a, p))


def P(a):
    return lambda x: gammainc(a, 0, x, regularized=True)


def Q(a):
    return lambda x: gammainc(a, x, inf, regularized=True)


emit(
    "gammaincinv", "@special.gammaincinv(a0[i], a1[i])", cases,
    [logx_bisect(P(mpf(a)), mpf(p), log_target=True) if p < 0.5 else logx_bisect(lambda x: -Q(mpf(a))(x), -(1 - mpf(p))) for a, p in cases],
    2e-16,
)
emit(
    "gammainccinv", "@special.gammainccinv(a0[i], a1[i])", cases,
    [logx_bisect(lambda x: -Q(mpf(a))(x), -mpf(p)) if p < 0.5 else logx_bisect(P(mpf(a)), 1 - mpf(p), log_target=True) for a, p in cases],
    2e-16,
)

# --- incomplete beta --------------------------------------------------------
cases = []
for _ in range(120):
    a, b = logu(1e-2, 1e4), logu(1e-2, 1e4)
    m = a / (a + b)
    x = rnd.random() if rnd.random() < 0.5 else min(max(m + rnd.uniform(-3, 3) * math.sqrt(a * b / (a + b) ** 3), 1e-12), 1 - 1e-12)
    cases.append((a, b, x))
emit("betainc", "@special.betainc(a0[i], a1[i], a2[i])", cases, [betainc(mpf(a), mpf(b), 0, mpf(x), regularized=True) for a, b, x in cases], 2e-16, floor=1e-300)
# 1 - I_x(a, b) as I_{1-x}(b, a): mpmath's betainc(a, b, x, 1) subtracts and
# returns 0 below about 1e-50
emit("betaincc", "@special.betaincc(a0[i], a1[i], a2[i])", cases, [betainc(mpf(b), mpf(a), 0, 1 - mpf(x), regularized=True) for a, b, x in cases], 2e-16, floor=1e-300)
cases = []
for _ in range(60):
    a, b = logu(1e-1, 1e3), logu(1e-1, 1e3)
    p = logu(1e-100, 0.5) if rnd.random() < 0.6 else 1 - logu(1e-14, 0.5)
    cases.append((a, b, p))


def I(a, b):
    return lambda x: betainc(a, b, 0, x, regularized=True)


emit(
    "betaincinv", "@special.betaincinv(a0[i], a1[i], a2[i])", cases,
    [logx_bisect(I(mpf(a), mpf(b)), mpf(p), hi=1, log_target=True) if p < 0.5 else 1 - logx_bisect(I(mpf(b), mpf(a)), 1 - mpf(p), hi=1, log_target=True) for a, b, p in cases],
    2e-16,
)

# --- the rest ---------------------------------------------------------------
cases = [(rnd.uniform(-8, 8), rnd.uniform(-5, 5)) for _ in range(40)]
emit(
    "owens_t", "@special.owens_t(a0[i], a1[i])", cases,
    [quad(lambda t: exp(-mpf(h) ** 2 * (1 + t * t) / 2) / (1 + t * t), [0, mpf(a)]) / (2 * pi) for h, a in cases],
    2e-16, floor=1e-300,
)
xs = [rnd.uniform(0.2, 6) for _ in range(40)]
emit(
    "kolmogorov_sf", "@special.kolmogorov_sf(a0[i])", xs,
    [2 * nsum(lambda k: (-1) ** (k - 1) * exp(-2 * k * k * mpf(x) ** 2), [1, inf]) for x in xs],
    2e-16, floor=1e-300,
)

header = "// Generated by scripts/gen_special_tests.py with mpmath; do not edit.\n\n"
(ROOT / "src/special/special_test.mbt").write_text(header + "\n\n".join(TESTS) + "\n")
print(len(TESTS), "tests", "(report mode)" if REPORT else "")

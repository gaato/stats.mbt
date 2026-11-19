///|
/// Standard forms written by hand for the catalogue's reference
/// distributions, through the special functions: each `std_` method is the
/// family at loc 0 and scale 1, and the generated code in catalog.mbt
/// applies location and scale to it and exposes the result as the trait
/// implementations. Field names follow SciPy (generated in catalog.mbt).

// ---------------------------------------------------------------------------
// Normal
// ---------------------------------------------------------------------------

///|
fn Normal::std_pdf(self : Normal, x : Double) -> Double {
  @math.exp(self.std_logpdf(x))
}

///|
/// Closed form, so the far tails stay finite on the log scale.
fn Normal::std_logpdf(_self : Normal, x : Double) -> Double {
  -0.5 * x * x - 0.5 * @math.ln(2.0 * @math.PI)
}

///|
fn Normal::std_cdf(_self : Normal, x : Double) -> Double {
  @special.ndtr(x)
}

///|
/// Symmetry gives the upper tail without cancellation.
fn Normal::std_sf(_self : Normal, x : Double) -> Double {
  @special.ndtr(-x)
}

///|
/// ln Φ from its own expansion, finite where Φ underflows.
fn Normal::std_logcdf(_self : Normal, x : Double) -> Double {
  @special.log_ndtr(x)
}

///|
fn Normal::std_logsf(_self : Normal, x : Double) -> Double {
  @special.log_ndtr(-x)
}

///|
fn Normal::std_quantile(_self : Normal, p : Double) -> Double {
  @special.ndtri(p)
}

///|
/// By symmetry, never through 1 - q.
fn Normal::std_isf(_self : Normal, q : Double) -> Double {
  -@special.ndtri(q)
}

///|
fn Normal::std_mean(_self : Normal) -> Double? {
  Some(0.0)
}

///|
fn Normal::std_variance(_self : Normal) -> Double? {
  Some(1.0)
}

///|
fn Normal::std_mode(_self : Normal) -> Double {
  0.0
}

///|
/// φ(t) = exp(-t²/2).
fn Normal::std_cf(_self : Normal, t : Double) -> (Double, Double) {
  (@math.exp(-0.5 * t * t), 0.0)
}

///|
/// M(t) = exp(t²/2), finite for every t.
fn Normal::std_mgf(_self : Normal, t : Double) -> Double? {
  Some(@math.exp(0.5 * t * t))
}

///|
/// Maximum likelihood: the sample mean and the divisor-n standard
/// deviation. Needs two points, since one gives scale = 0.
pub impl Fit for Normal with fn fit(data) {
  if data.length() < 2 {
    raise TooFewPoints(needed=2, got=data.length())
  }
  ...
}

// ---------------------------------------------------------------------------
// Student t
// ---------------------------------------------------------------------------

///|
fn StudentT::std_mean(self : StudentT) -> Double? {
  if self.df > 1.0 {
    Some(0.0)
  } else {
    None
  }
}

///|
fn StudentT::std_variance(self : StudentT) -> Double? {
  if self.df > 2.0 {
    Some(self.df / (self.df - 2.0))
  } else if self.df > 1.0 {
    Some(1.0 / 0.0)
  } else {
    None
  }
}

///|
fn StudentT::std_pdf(self : StudentT, x : Double) -> Double {
  @math.exp(self.std_logpdf(x))
}

///|
/// -ln √ν - ln B(1/2, ν/2) - (ν+1)/2 · ln(1 + x²/ν).
fn StudentT::std_logpdf(self : StudentT, x : Double) -> Double {
  -0.5 * @math.ln(self.df) -
  @special.lbeta(0.5, 0.5 * self.df) -
  0.5 * (self.df + 1.0) * @math.ln_1p(x * x / self.df)
}

///|
/// Half the lower tail, I_t(ν/2, 1/2)/2 with t = ν/(ν + x²); both t and
/// 1 - t = x²/(ν + x²) are formed directly, so neither end cancels.
fn StudentT::half_tail(self : StudentT, x : Double) -> Double {
  let d = self.df + x * x
  let (i, _) = @special.betainc_pair(0.5 * self.df, 0.5, self.df / d, x * x / d)
  0.5 * i
}

///|
fn StudentT::std_cdf(self : StudentT, x : Double) -> Double {
  if x < 0.0 {
    self.half_tail(x)
  } else {
    1.0 - self.half_tail(x)
  }
}

///|
/// The p-value path of every t test, from the incomplete beta directly.
fn StudentT::std_sf(self : StudentT, x : Double) -> Double {
  if x > 0.0 {
    self.half_tail(x)
  } else {
    1.0 - self.half_tail(x)
  }
}

///|
/// Inverts the tail it is in: I_t(ν/2, 1/2) = 2q for the far tail, and the
/// complementary I_z(1/2, ν/2) = 1 - 2q near the centre, so x stays
/// accurate both at 1e-300 and next to the median.
fn StudentT::std_quantile(self : StudentT, p : Double) -> Double {
  if p == 0.5 {
    return 0.0
  }
  let q = if p < 0.5 { p } else { 1.0 - p }
  let target = 2.0 * q
  let x = if target <= 0.5 {
    let t = @special.betaincinv(0.5 * self.df, 0.5, target)
    (self.df * (1.0 - t) / t).sqrt()
  } else {
    let z = @special.betaincinv(0.5, 0.5 * self.df, 1.0 - target)
    (self.df * z / (1.0 - z)).sqrt()
  }
  if p < 0.5 {
    -x
  } else {
    x
  }
}

///|
fn StudentT::std_isf(self : StudentT, q : Double) -> Double {
  -self.std_quantile(q)
}

///|
/// No `Mgf`: E[exp(tX)] diverges for every t != 0. The characteristic
/// function is a modified Bessel function of the second kind.
fn StudentT::std_cf(self : StudentT, t : Double) -> (Double, Double) {
  ignore((self, t))
  ...
}

// ---------------------------------------------------------------------------
// Gamma: shape a, with the scale every family has; the rate form is an
// alternative constructor.
// ---------------------------------------------------------------------------

///|
/// The same distribution by shape and rate (scale = 1/rate, loc = 0).
pub fn Gamma::with_rate(a : Double, rate~ : Double) -> Gamma raise ParamError {
  if !(rate > 0.0) || rate.is_inf() {
    raise InvalidParam(name="rate", value=rate, reason="must be in (0, inf)")
  }
  Gamma::new(a, scale=1.0 / rate)
}

///|
/// Gamma by shape and rate: a separate family value, so a fit reports its
/// parameters in the convention asked for. The location stays at 0.
pub fn Gamma::family_rate() -> Family[Gamma] {
  {
    name: "gamma (a, rate)",
    specs: [
      { name: "a", constraint: Positive, },
      { name: "rate", constraint: Positive, },
    ],
    make: p => Gamma::with_rate(p[0], rate=p[1]),
    params: d => [d.a, 1.0 / d.scale],
  }
}

///|
fn Gamma::std_cdf(self : Gamma, x : Double) -> Double {
  if x <= 0.0 {
    0.0
  } else {
    @special.gammainc(self.a, x)
  }
}

///|
fn Gamma::std_sf(self : Gamma, x : Double) -> Double {
  if x <= 0.0 {
    1.0
  } else {
    @special.gammaincc(self.a, x)
  }
}

///|
fn Gamma::std_pdf(self : Gamma, x : Double) -> Double {
  if x < 0.0 {
    0.0
  } else {
    @math.exp(self.std_logpdf(x))
  }
}

///|
fn Gamma::std_logpdf(self : Gamma, x : Double) -> Double {
  if x < 0.0 {
    -1.0 / 0.0
  } else {
    (self.a - 1.0) * @math.ln(x) - x - @special.lgamma(self.a)
  }
}

///|
fn Gamma::std_quantile(self : Gamma, p : Double) -> Double {
  @special.gammaincinv(self.a, p)
}

///|
fn Gamma::std_isf(self : Gamma, q : Double) -> Double {
  @special.gammainccinv(self.a, q)
}

///|
fn Gamma::std_mean(self : Gamma) -> Double? {
  Some(self.a)
}

///|
fn Gamma::std_variance(self : Gamma) -> Double? {
  Some(self.a)
}

///|
/// M(t) = (1 - t)^(-a) for t < 1; it diverges from there on.
fn Gamma::std_mgf(self : Gamma, t : Double) -> Double? {
  if t < 1.0 {
    Some(@math.pow(1.0 - t, -self.a))
  } else {
    None
  }
}

// ---------------------------------------------------------------------------
// Reference distributions for tests and power
// ---------------------------------------------------------------------------

///|
/// P(X <= x) = P(df/2, x/2), the lower regularised incomplete gamma.
fn ChiSquared::std_cdf(self : ChiSquared, x : Double) -> Double {
  if x <= 0.0 {
    0.0
  } else {
    @special.gammainc(0.5 * self.df, 0.5 * x)
  }
}

///|
/// P(X > x) = Q(df/2, x/2), computed directly for small p-values.
fn ChiSquared::std_sf(self : ChiSquared, x : Double) -> Double {
  if x <= 0.0 {
    1.0
  } else {
    @special.gammaincc(0.5 * self.df, 0.5 * x)
  }
}

///|
fn ChiSquared::std_pdf(self : ChiSquared, x : Double) -> Double {
  if x < 0.0 {
    0.0
  } else {
    @math.exp(self.std_logpdf(x))
  }
}

///|
fn ChiSquared::std_logpdf(self : ChiSquared, x : Double) -> Double {
  if x < 0.0 {
    -1.0 / 0.0
  } else {
    let k = 0.5 * self.df
    (k - 1.0) * @math.ln(x) - 0.5 * x - k * @math.ln(2.0) - @special.lgamma(k)
  }
}

///|
fn ChiSquared::std_quantile(self : ChiSquared, p : Double) -> Double {
  2.0 * @special.gammaincinv(0.5 * self.df, p)
}

///|
fn ChiSquared::std_isf(self : ChiSquared, q : Double) -> Double {
  2.0 * @special.gammainccinv(0.5 * self.df, q)
}

///|
/// (P(X <= x), P(X > x)) = I_u(dfn/2, dfd/2) and its complement, with
/// u = dfn x/(dfn x + dfd) and 1 - u = dfd/(dfn x + dfd) both formed directly.
fn FisherF::tails(self : FisherF, x : Double) -> (Double, Double) {
  let d = self.dfn * x + self.dfd
  @special.betainc_pair(
    0.5 * self.dfn,
    0.5 * self.dfd,
    self.dfn * x / d,
    self.dfd / d,
  )
}

///|
fn FisherF::std_cdf(self : FisherF, x : Double) -> Double {
  if x <= 0.0 {
    0.0
  } else {
    self.tails(x).0
  }
}

///|
fn FisherF::std_sf(self : FisherF, x : Double) -> Double {
  if x <= 0.0 {
    1.0
  } else {
    self.tails(x).1
  }
}

///|
fn FisherF::std_pdf(self : FisherF, x : Double) -> Double {
  if x < 0.0 {
    0.0
  } else {
    @math.exp(self.std_logpdf(x))
  }
}

///|
fn FisherF::std_logpdf(self : FisherF, x : Double) -> Double {
  if x < 0.0 {
    -1.0 / 0.0
  } else {
    let m = self.dfn
    let n = self.dfd
    0.5 * m * @math.ln(m / n) +
    (0.5 * m - 1.0) * @math.ln(x) -
    0.5 * (m + n) * @math.ln_1p(m * x / n) -
    @special.lbeta(0.5 * m, 0.5 * n)
  }
}

///|
/// The lower tail inverts I_u(dfn/2, dfd/2); above the median the upper
/// tail is inverted instead, so a small 1 - u is never formed.
fn FisherF::std_quantile(self : FisherF, p : Double) -> Double {
  if p > 0.5 {
    return self.std_isf(1.0 - p)
  }
  let u = @special.betaincinv(0.5 * self.dfn, 0.5 * self.dfd, p)
  self.dfd * u / (self.dfn * (1.0 - u))
}

///|
fn FisherF::std_isf(self : FisherF, q : Double) -> Double {
  if q > 0.5 {
    return self.std_quantile(1.0 - q)
  }
  let v = @special.betaincinv(0.5 * self.dfd, 0.5 * self.dfn, q)
  self.dfd * (1.0 - v) / (self.dfn * v)
}

///|
fn Beta::std_pdf(self : Beta, x : Double) -> Double {
  if x < 0.0 || x > 1.0 {
    0.0
  } else {
    @math.exp(self.std_logpdf(x))
  }
}

///|
fn Beta::std_logpdf(self : Beta, x : Double) -> Double {
  if x < 0.0 || x > 1.0 {
    -1.0 / 0.0
  } else {
    (self.a - 1.0) * @math.ln(x) +
    (self.b - 1.0) * @math.ln_1p(-x) -
    @special.lbeta(self.a, self.b)
  }
}

///|
fn Beta::std_sf(self : Beta, x : Double) -> Double {
  if x <= 0.0 {
    1.0
  } else if x >= 1.0 {
    0.0
  } else {
    @special.betaincc(self.a, self.b, x)
  }
}

///|
/// 1 - I_x(a, b) = q is I_{1-x}(b, a) = q.
fn Beta::std_isf(self : Beta, q : Double) -> Double {
  1.0 - @special.betaincinv(self.b, self.a, q)
}

///|
fn Beta::std_cdf(self : Beta, x : Double) -> Double {
  if x <= 0.0 {
    0.0
  } else if x >= 1.0 {
    1.0
  } else {
    @special.betainc(self.a, self.b, x)
  }
}

///|
/// Clopper–Pearson intervals for a proportion are quantiles of this.
fn Beta::std_quantile(self : Beta, p : Double) -> Double {
  @special.betaincinv(self.a, self.b, p)
}

///|

///|
/// More accumulators. Each is a monoid like `MeanVar`, so it can be
/// computed in parallel or over a stream and merged.

///|
/// Count, means and the co-moment of (x, y): covariance and Pearson
/// correlation in one pass. Values are accumulated relative to a shift (the
/// first point seen), so a large common offset such as 1e9 does not cost
/// digits; plain Welford loses about 1e-9 relative there.
pub struct Comoment {
  priv n : Int64
  priv kx : Double
  priv ky : Double
  priv mx : Double
  priv my : Double
  priv m2x : Double
  priv m2y : Double
  priv cxy : Double
}

///|
pub fn Comoment::empty() -> Comoment {
  { n: 0, kx: 0.0, ky: 0.0, mx: 0.0, my: 0.0, m2x: 0.0, m2y: 0.0, cxy: 0.0 }
}

///|
pub fn Comoment::push(self : Comoment, x : Double, y : Double) -> Comoment {
  let (kx, ky) = if self.n == 0 { (x, y) } else { (self.kx, self.ky) }
  let u = x - kx
  let v = y - ky
  let n = self.n + 1
  let fn_ = n.to_double()
  let du = u - self.mx
  let dv = v - self.my
  let mx = self.mx + du / fn_
  let my = self.my + dv / fn_
  {
    n,
    kx,
    ky,
    mx,
    my,
    m2x: self.m2x + du * (u - mx),
    m2y: self.m2y + dv * (v - my),
    cxy: self.cxy + du * (v - my),
  }
}

///|
/// Chan's pairwise update, with b's means moved onto a's shift first.
pub fn Comoment::merge(a : Comoment, b : Comoment) -> Comoment {
  if a.n == 0 {
    return b
  }
  if b.n == 0 {
    return a
  }
  let n = a.n + b.n
  let fa = a.n.to_double()
  let fb = b.n.to_double()
  let fn_ = n.to_double()
  // difference of the true means, b - a, computed from the shifts and the
  // relative means separately
  let dx = b.kx - a.kx + (b.mx - a.mx)
  let dy = b.ky - a.ky + (b.my - a.my)
  {
    n,
    kx: a.kx,
    ky: a.ky,
    mx: a.mx + dx * fb / fn_,
    my: a.my + dy * fb / fn_,
    m2x: a.m2x + b.m2x + dx * dx * fa * fb / fn_,
    m2y: a.m2y + b.m2y + dy * dy * fa * fb / fn_,
    cxy: a.cxy + b.cxy + dx * dy * fa * fb / fn_,
  }
}

///|
pub fn Comoment::mean_x(self : Comoment) -> Double? {
  if self.n == 0 {
    None
  } else {
    Some(self.kx + self.mx)
  }
}

///|
pub fn Comoment::covariance(self : Comoment, ddof~ : Int) -> Double? {
  let d = self.n - ddof.to_int64()
  if d <= 0 {
    None
  } else {
    Some(self.cxy / d.to_double())
  }
}

///|
/// Pearson correlation; `None` when either variable is constant.
pub fn Comoment::correlation(self : Comoment) -> Double? {
  if self.n < 2 || self.m2x <= 0.0 || self.m2y <= 0.0 {
    None
  } else {
    Some(self.cxy / (self.m2x * self.m2y).sqrt())
  }
}

///|
/// Mean and variance under weights whose meaning is stated: frequency
/// weights give n = Σw, reliability weights the effective n
/// (Σw)²/Σw², sampling weights a design-based variance.
pub struct WeightedMeanVar {
  priv kind : @data.WeightKind
  priv w : Double
  priv w2 : Double
  priv mean : Double
  priv s : Double
}

///|
pub fn WeightedMeanVar::empty(kind : @data.WeightKind) -> WeightedMeanVar {
  { kind, w: 0.0, w2: 0.0, mean: 0.0, s: 0.0 }
}

///|
/// West's weighted update.
pub fn WeightedMeanVar::push(
  self : WeightedMeanVar,
  x : Double,
  weight : Double,
) -> WeightedMeanVar {
  if weight == 0.0 {
    return self
  }
  let w = self.w + weight
  let d = x - self.mean
  let mean = self.mean + d * weight / w
  {
    kind: self.kind,
    w,
    w2: self.w2 + weight * weight,
    mean,
    s: self.s + weight * d * (x - mean),
  }
}

///|
pub fn WeightedMeanVar::mean(self : WeightedMeanVar) -> Double? {
  if self.w == 0.0 {
    None
  } else {
    Some(self.mean)
  }
}

///|
/// The unbiased variance for the weight kind.
pub fn WeightedMeanVar::variance(self : WeightedMeanVar) -> Double? {
  ignore(self)
  ...
}

///|
/// Third and fourth moments as well. The bias-corrected sample skewness
/// and excess kurtosis differ from the plain moment ratios, and libraries
/// disagree on the default, so the caller chooses.
pub struct Moments4 {
  priv n : Int64
  priv mean : Double
  priv m2 : Double
  priv m3 : Double
  priv m4 : Double
}

///|
pub fn Moments4::empty() -> Moments4 {
  { n: 0, mean: 0.0, m2: 0.0, m3: 0.0, m4: 0.0 }
}

///|
pub fn Moments4::push(self : Moments4, x : Double) -> Moments4 {
  ignore((self, x))
  ...
}

///|
pub fn Moments4::merge(a : Moments4, b : Moments4) -> Moments4 {
  ignore((a, b))
  ...
}

///|
pub fn Moments4::skewness(self : Moments4, bias_corrected~ : Bool) -> Double? {
  ignore((self, bias_corrected))
  ...
}

///|
/// Excess kurtosis (0 for the normal).
pub fn Moments4::excess_kurtosis(
  self : Moments4,
  bias_corrected~ : Bool,
) -> Double? {
  ignore((self, bias_corrected))
  ...
}

///|
/// Spearman's rank correlation: Pearson on ranks with average ties.
pub fn spearman(xy : @data.Bivariate) -> Double? {
  ignore(xy)
  ...
}

///|
/// Kendall's tau-b, with the tie adjustment.
pub fn kendall_tau_b(xy : @data.Bivariate) -> Double? {
  ignore(xy)
  ...
}

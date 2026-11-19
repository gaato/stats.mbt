///|
/// The gamma function family: ln Γ, ψ, ψ', ln B and ln C(n, k).
///
/// ln Γ is accurate in relative terms also next to its zeros at 1 and 2,
/// where a series in z = x - 2 with the coefficients ζ(k) - 1 replaces the
/// subtraction of nearly equal logarithms; ψ likewise uses its Taylor
/// series around its positive zero. Constants come from constants.mbt.

///|
/// ln Γ(2 + z) for |z| <= 1/2:
/// z(1 - γ) + Σ_{k>=2} (-1)^k (ζ(k) - 1) z^k / k.
fn lgamma_2p(z : Double) -> Double {
  let mut sum = z * (1.0 - euler_gamma)
  // zk runs through (-1)^k z^k, starting before k = 2
  let mut zk = -z
  for i in 0..<zeta_minus_one.length() {
    let k = i + 2
    zk = zk * -z
    let term = zeta_minus_one[i] * zk / k.to_double()
    sum = sum + term
    if term.abs() <= 1.0e-17 * sum.abs() {
      break
    }
  }
  sum
}

///|
/// ln Γ(1 + z) for |z| <= 1/2, taking z itself so that a tiny z keeps its
/// digits (1 + z would round them away).
fn lgamma_1p(z : Double) -> Double {
  lgamma_2p(z) - @math.ln_1p(z)
}

///|
/// The remainder of Stirling's series, ln Γ(x) - [(x - 1/2) ln x - x +
/// ln √(2π)], for x >= 8.
fn stirling_delta(x : Double) -> Double {
  let r = 1.0 / (x * x)
  let mut s = 0.0
  for i = stirling.length() - 1; i >= 0; i = i - 1 {
    s = s * r + stirling[i]
  }
  s / x
}

///|
/// sin(πx) with the argument reduced exactly first.
fn sinpi(x : Double) -> Double {
  let r = x - 2.0 * @math.round(0.5 * x)
  // r in [-1, 1]; fold to [-1/2, 1/2] by symmetry
  if r > 0.5 {
    @math.sin(@math.PI * (1.0 - r))
  } else if r < -0.5 {
    -@math.sin(@math.PI * (1.0 + r))
  } else {
    @math.sin(@math.PI * r)
  }
}

///|
/// ln |Γ(x)|. +inf at the poles (0 and the negative integers).
pub fn lgamma(x : Double) -> Double {
  if x.is_nan() {
    return x
  }
  if x.is_inf() {
    return 1.0 / 0.0
  }
  if x <= 0.0 {
    if x == @math.floor(x) {
      return 1.0 / 0.0
    }
    // reflection: Γ(x) Γ(1 - x) = π / sin(πx)
    return @math.ln(@math.PI / sinpi(x).abs()) - lgamma(1.0 - x)
  }
  if x < 0.5 {
    lgamma_1p(x) - @math.ln(x)
  } else if x < 1.5 {
    lgamma_1p(x - 1.0)
  } else if x < 2.5 {
    lgamma_2p(x - 2.0)
  } else if x < 8.0 {
    // move down into [1.5, 2.5); the logarithms of the factors are added,
    // never subtracted
    let m = @math.floor(x - 1.5)
    let y = x - m
    let mut prod = 1.0
    let mut i = 0.0
    while i < m {
      prod = prod * (y + i)
      i = i + 1.0
    }
    lgamma_2p(y - 2.0) + @math.ln(prod)
  } else {
    (x - 0.5) * @math.ln(x) - x + ln_sqrt_2pi + stirling_delta(x)
  }
}

///|
/// ψ(x) = Γ'(x)/Γ(x). NaN at the poles.
pub fn digamma(x : Double) -> Double {
  if x.is_nan() {
    return x
  }
  if x == 1.0 / 0.0 {
    return x
  }
  if x <= 0.0 {
    if x == @math.floor(x) {
      return 0.0 / 0.0
    }
    // reflection: ψ(1 - x) - ψ(x) = π cot(πx)
    let r = x - @math.round(x)
    return digamma(1.0 - x) - @math.PI / @math.tan(@math.PI * r)
  }
  let z = x - digamma_root_hi - digamma_root_lo
  if z.abs() < 0.25 {
    let mut s = 0.0
    for i = digamma_root_taylor.length() - 1; i >= 0; i = i - 1 {
      s = s * z + digamma_root_taylor[i]
    }
    return s * z
  }
  let mut acc = 0.0
  let mut y = x
  while y < 10.0 {
    acc = acc + 1.0 / y
    y = y + 1.0
  }
  let r = 1.0 / (y * y)
  let mut s = 0.0
  for i = digamma_asym.length() - 1; i >= 0; i = i - 1 {
    s = s * r + digamma_asym[i]
  }
  @math.ln(y) - 0.5 / y - s * r - acc
}

///|
/// ψ'(x), the trigamma function.
pub fn trigamma(x : Double) -> Double {
  if x.is_nan() {
    return x
  }
  if x <= 0.0 {
    if x == @math.floor(x) {
      return 1.0 / 0.0
    }
    // reflection: ψ'(1 - x) + ψ'(x) = π² / sin²(πx)
    let s = sinpi(x)
    return @math.PI * @math.PI / (s * s) - trigamma(1.0 - x)
  }
  let mut acc = 0.0
  let mut y = x
  while y < 10.0 {
    acc = acc + 1.0 / (y * y)
    y = y + 1.0
  }
  let r = 1.0 / (y * y)
  let mut s = 0.0
  for i = trigamma_asym.length() - 1; i >= 0; i = i - 1 {
    s = s * r + trigamma_asym[i]
  }
  acc + 1.0 / y + 0.5 * r + s * r / y
}

///|
/// ln B(a, b) for a, b > 0. When an argument is large the Stirling forms
/// are combined before taking differences, so ln B(10⁶, 0.5) keeps its
/// digits instead of being the difference of two numbers near 1.3e7.
pub fn lbeta(a : Double, b : Double) -> Double {
  if a.is_nan() || b.is_nan() || a <= 0.0 || b <= 0.0 {
    return 0.0 / 0.0
  }
  let (big, small) = if a >= b { (a, b) } else { (b, a) }
  if small >= 8.0 {
    // ln √(2π) - ln √small + (big - 1/2) ln(big/(a+b)) + small ln(small/(a+b))
    ln_sqrt_2pi -
    0.5 * @math.ln(small) -
    (big - 0.5) * @math.ln_1p(small / big) -
    small * @math.ln_1p(big / small) +
    stirling_delta(big) +
    stirling_delta(small) -
    stirling_delta(big + small)
  } else if big >= 8.0 {
    // ln Γ(small) + [ln Γ(big) - ln Γ(big + small)]
    lgamma(small) -
    small * @math.ln(big) -
    (big + small - 0.5) * @math.ln_1p(small / big) +
    small +
    stirling_delta(big) -
    stirling_delta(big + small)
  } else {
    lgamma(a) + lgamma(b) - lgamma(a + b)
  }
}

///|
/// ln C(n, k) for real 0 <= k <= n.
pub fn lchoose(n : Double, k : Double) -> Double {
  if n.is_nan() || k.is_nan() || k < 0.0 || k > n {
    return 0.0 / 0.0
  }
  if k == 0.0 || k == n {
    return 0.0
  }
  -@math.ln_1p(n) - lbeta(n - k + 1.0, k + 1.0)
}

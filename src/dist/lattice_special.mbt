///|
/// The discrete reference distributions, through the special functions:
/// binomial and Poisson tails are incomplete beta and gamma functions, each
/// tail from its own one, and the hypergeometric tails are summed on their
/// own sides. Written on integer-valued Doubles like lattice.mbt.

// ---------------------------------------------------------------------------
// Binomial
// ---------------------------------------------------------------------------

///|
fn Binomial::lattice_pmf(self : Binomial, k : Double) -> Double {
  let n = self.n.to_double()
  if k < 0.0 || k > n {
    0.0
  } else if self.p == 0.0 {
    if k == 0.0 {
      1.0
    } else {
      0.0
    }
  } else if self.p == 1.0 {
    if k == n {
      1.0
    } else {
      0.0
    }
  } else {
    @math.exp(
      @special.lchoose(n, k) +
      k * @math.ln(self.p) +
      (n - k) * @math.ln_1p(-self.p),
    )
  }
}

///|
/// P(X <= k) = I_{1-p}(n - k, k + 1).
fn Binomial::lattice_cdf(self : Binomial, k : Double) -> Double {
  let n = self.n.to_double()
  if k < 0.0 {
    0.0
  } else if k >= n {
    1.0
  } else {
    @special.betainc(n - k, k + 1.0, 1.0 - self.p)
  }
}

///|
/// P(X > k) = I_p(k + 1, n - k).
fn Binomial::lattice_sf(self : Binomial, k : Double) -> Double {
  let n = self.n.to_double()
  if k < 0.0 {
    1.0
  } else if k >= n {
    0.0
  } else {
    @special.betainc(k + 1.0, n - k, self.p)
  }
}

///|
/// A normal-approximation start, corrected by exact comparisons.
fn Binomial::lattice_quantile(self : Binomial, q : Double) -> Double {
  let n = self.n.to_double()
  let guess = n * self.p +
    @special.ndtri(q) * (n * self.p * (1.0 - self.p)).sqrt()
  refine_quantile(k => self.lattice_cdf(k), q, guess, 0.0, n)
}

///|
fn Binomial::lattice_support(self : Binomial) -> (Double, Double) {
  (0.0, self.n.to_double())
}

///|
pub impl Moments for Binomial with fn mean(self) {
  Some(self.n.to_double() * self.p)
}

///|
pub impl Moments for Binomial with fn variance(self) {
  Some(self.n.to_double() * self.p * (1.0 - self.p))
}

///|
/// floor((n + 1) p); when (n + 1) p is an integer the mass is the same one
/// step below, and the smaller point is returned.
pub impl Mode for Binomial with fn mode(self) {
  let n = self.n.to_double()
  let m = (n + 1.0) * self.p
  if self.p == 1.0 {
    n
  } else if m == @math.floor(m) && m > 0.0 {
    m - 1.0
  } else {
    @math.floor(m)
  }
}

///|
/// φ(t) = (1 - p + p e^{it})^n, raised to the n-th power in polar form.
pub impl Characteristic for Binomial with fn cf(self, t) {
  let re = 1.0 - self.p + self.p * @math.cos(t)
  let im = self.p * @math.sin(t)
  let n = self.n.to_double()
  let r = @math.pow(@math.hypot(re, im), n)
  let theta = n * @math.atan2(im, re)
  (r * @math.cos(theta), r * @math.sin(theta))
}

// ---------------------------------------------------------------------------
// Poisson
// ---------------------------------------------------------------------------

///|
fn Poisson::lattice_pmf(self : Poisson, k : Double) -> Double {
  if k < 0.0 {
    0.0
  } else if self.mu == 0.0 {
    if k == 0.0 {
      1.0
    } else {
      0.0
    }
  } else {
    @math.exp(k * @math.ln(self.mu) - self.mu - @special.lgamma(k + 1.0))
  }
}

///|
/// P(X <= k) = Q(k + 1, mu).
fn Poisson::lattice_cdf(self : Poisson, k : Double) -> Double {
  if k < 0.0 {
    0.0
  } else {
    @special.gammaincc(k + 1.0, self.mu)
  }
}

///|
/// P(X > k) = P(k + 1, mu).
fn Poisson::lattice_sf(self : Poisson, k : Double) -> Double {
  if k < 0.0 {
    1.0
  } else {
    @special.gammainc(k + 1.0, self.mu)
  }
}

///|
fn Poisson::lattice_quantile(self : Poisson, q : Double) -> Double {
  let guess = self.mu + @special.ndtri(q) * self.mu.sqrt()
  refine_quantile(k => self.lattice_cdf(k), q, guess, 0.0, 1.0 / 0.0)
}

///|
fn Poisson::lattice_support(_self : Poisson) -> (Double, Double) {
  (0.0, 1.0 / 0.0)
}

///|
pub impl Moments for Poisson with fn mean(self) {
  Some(self.mu)
}

///|
pub impl Moments for Poisson with fn variance(self) {
  Some(self.mu)
}

///|
/// floor(mu); for an integer mu the mass at mu - 1 is the same, and the
/// smaller point is returned.
pub impl Mode for Poisson with fn mode(self) {
  let m = @math.floor(self.mu)
  if m == self.mu && m > 0.0 {
    m - 1.0
  } else {
    m
  }
}

// ---------------------------------------------------------------------------
// Hypergeometric: the null distribution of Fisher's exact test for a
// 2 × 2 table.
// ---------------------------------------------------------------------------

///|
fn Hypergeometric::lowest(self : Hypergeometric) -> Double {
  let k = (self.draws - (self.population - self.successes)).to_double()
  if k > 0.0 {
    k
  } else {
    0.0
  }
}

///|
fn Hypergeometric::highest(self : Hypergeometric) -> Double {
  if self.successes < self.draws {
    self.successes.to_double()
  } else {
    self.draws.to_double()
  }
}

///|
fn Hypergeometric::lattice_pmf(self : Hypergeometric, k : Double) -> Double {
  if k < self.lowest() || k > self.highest() {
    0.0
  } else {
    let population = self.population.to_double()
    let successes = self.successes.to_double()
    let draws = self.draws.to_double()
    @math.exp(
      @special.lchoose(successes, k) +
      @special.lchoose(population - successes, draws - k) -
      @special.lchoose(population, draws),
    )
  }
}

///|
/// Σ pmf over the integers from..=to.
fn Hypergeometric::partial(
  self : Hypergeometric,
  from : Double,
  to : Double,
) -> Double {
  let mut s = @num.CompensatedSum::zero()
  let mut j = from
  while j <= to {
    s = s.add(self.lattice_pmf(j))
    j = j + 1.0
  }
  s.value()
}

///|
/// Each tail is summed on its own side, so neither is 1 minus the other.
fn Hypergeometric::lattice_cdf(self : Hypergeometric, k : Double) -> Double {
  if k < self.lowest() {
    0.0
  } else if k >= self.highest() {
    1.0
  } else {
    self.partial(self.lowest(), k)
  }
}

///|
fn Hypergeometric::lattice_sf(self : Hypergeometric, k : Double) -> Double {
  if k < self.lowest() {
    1.0
  } else if k >= self.highest() {
    0.0
  } else {
    self.partial(k + 1.0, self.highest())
  }
}

///|
fn Hypergeometric::lattice_quantile(
  self : Hypergeometric,
  q : Double,
) -> Double {
  let guess = self.draws.to_double() *
    self.successes.to_double() /
    self.population.to_double()
  refine_quantile(
    k => self.lattice_cdf(k),
    q,
    guess,
    self.lowest(),
    self.highest(),
  )
}

///|
fn Hypergeometric::lattice_support(self : Hypergeometric) -> (Double, Double) {
  (self.lowest(), self.highest())
}

///|
pub impl Moments for Hypergeometric with fn mean(self) {
  Some(
    self.draws.to_double() *
    self.successes.to_double() /
    self.population.to_double(),
  )
}

///|
pub impl Moments for Hypergeometric with fn variance(self) {
  let population = self.population.to_double()
  if population <= 1.0 {
    return Some(0.0)
  }
  let share = self.successes.to_double() / population
  Some(
    self.draws.to_double() *
    share *
    (1.0 - share) *
    (population - self.draws.to_double()) /
    (population - 1.0),
  )
}

///|
/// floor((draws + 1)(successes + 1)/(population + 2)), kept inside the
/// support; when that ratio is an integer the mass is the same one step
/// below, and the smaller point is returned.
pub impl Mode for Hypergeometric with fn mode(self) {
  let r = (self.draws.to_double() + 1.0) *
    (self.successes.to_double() + 1.0) /
    (self.population.to_double() + 2.0)
  let m = if r == @math.floor(r) { r - 1.0 } else { @math.floor(r) }
  if m < self.lowest() {
    self.lowest()
  } else if m > self.highest() {
    self.highest()
  } else {
    m
  }
}

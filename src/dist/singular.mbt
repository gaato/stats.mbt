///|
/// The Cantor distribution: no density and no atoms, only a singular
/// part. It has a distribution function, a quantile function and a
/// sampler, and nothing else, which makes it the strictest test that the
/// capability traits are split correctly.
pub struct Cantor {}

///|
pub fn Cantor::new() -> Cantor {
  Cantor::{ }
}

///|
/// The Cantor function, evaluated digit by digit in base 3.
pub impl Cdf for Cantor with fn cdf(self, x) {
  ignore((self, x))
  ...
}

///|
pub impl Quantile for Cantor with fn quantile(self, p) {
  ignore((self, p))
  ...
}

///|
/// Random base-3 digits from {0, 2}.
pub impl Sample for Cantor with fn sample(self, rng) {
  ignore((self, rng))
  ...
}

///|
pub impl Univariate for Cantor

///|
/// Mean 1/2 and variance 1/8, although there is no density to integrate.
pub impl Moments for Cantor with fn mean(_self) {
  Some(0.5)
}

///|
pub impl Moments for Cantor with fn variance(_self) {
  Some(0.125)
}

///|
/// φ(t) = exp(it/2) Π cos(t / 3^k), k >= 1.
pub impl Characteristic for Cantor with fn cf(self, t) {
  ignore((self, t))
  ...
}

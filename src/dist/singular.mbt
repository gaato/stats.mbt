///|
/// The Cantor distribution: no density and no atoms, only a singular
/// part. It has a distribution function, a quantile function and a
/// sampler, and nothing else, which makes it the strictest test that the
/// capability traits are split correctly.
pub struct Cantor {}

///|
/// The Cantor distribution on [0, 1].
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
/// The Cantor function from the right: digits 0 and 2 exchanged.
pub impl Cdf for Cantor with fn sf(self, x) {
  ignore((self, x))
  ...
}

///|
/// The smallest x with P(X ≤ x) ≥ p.
pub impl Quantile for Cantor with fn quantile(self, p) {
  ignore((self, p))
  ...
}

///|
/// The smallest x with P(X > x) ≤ q.
pub impl Quantile for Cantor with fn isf(self, q) {
  ignore((self, q))
  ...
}

///|
/// Random base-3 digits from {0, 2}.
pub impl Sample for Cantor with fn sample(self, rng) {
  ignore((self, rng))
  ...
}

///|
/// Usable wherever a `&Univariate` is taken.
pub impl Univariate for Cantor

///|
/// The lower end of the support.
pub impl Support for Cantor with fn support_lower(_self) {
  0.0
}

///|
/// The upper end of the support.
pub impl Support for Cantor with fn support_upper(_self) {
  1.0
}

///|
/// Mean 1/2 and variance 1/8, although there is no density to integrate.
pub impl Moments for Cantor with fn mean(_self) {
  Some(0.5)
}

///|
/// The variance, or `None` when it does not exist.
pub impl Moments for Cantor with fn variance(_self) {
  Some(0.125)
}

///|
/// φ(t) = exp(it/2) Π cos(t / 3^k), k >= 1.
pub impl Characteristic for Cantor with fn cf(self, t) {
  ignore((self, t))
  ...
}

///|
/// Dot syntax for `Cdf`: `d.cdf`, `d.sf`, `d.logcdf`, `d.logsf`, `d.sf_at_least`.
pub extend Cantor with Cdf::{cdf, sf, logcdf, logsf, sf_at_least}

///|
/// Dot syntax for `Quantile`: `d.quantile`, `d.isf`.
pub extend Cantor with Quantile::{quantile, isf}

///|
/// Dot syntax for `Sample`: `d.sample`, `d.sample_into`.
pub extend Cantor with Sample::{sample, sample_into}

///|
/// Dot syntax for `Moments`: `d.mean`, `d.variance`.
pub extend Cantor with Moments::{mean, variance}

///|
/// Dot syntax for `Support`: `d.support_lower`, `d.support_upper`.
pub extend Cantor with Support::{support_lower, support_upper}

///|
/// Dot syntax for `Characteristic`: `d.cf`.
pub extend Cantor with Characteristic::{cf}

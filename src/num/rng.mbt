///|
/// Random numbers, below the distributions.
///
/// Reproducibility contract: for a given generator, seed and variate
/// method, the stream of values never changes between library versions.
/// A faster or better algorithm is added as a new enum case, never swapped
/// in under an existing one. (NumPy had to keep its legacy RandomState for
/// exactly this reason.)

///|
/// A source of uniform 64-bit words. Distributions take it as a trait
/// object, so no global random state exists anywhere in the library.
pub(open) trait RandomSource {
  fn next_u64(Self) -> UInt64
}

///|
/// xoshiro256** (Blackman and Vigna), seeded through splitmix64.
pub struct Xoshiro256 {
  priv mut s0 : UInt64
  priv mut s1 : UInt64
  priv mut s2 : UInt64
  priv mut s3 : UInt64
}

///|
fn rotl(x : UInt64, k : Int32) -> UInt64 {
  (x << k) | (x >> (64 - k))
}

///|
fn splitmix64(state : UInt64) -> (UInt64, UInt64) {
  let s = state + 0x9e3779b97f4a7c15UL
  let mut z = s
  z = (z ^ (z >> 30)) * 0xbf58476d1ce4e5b9UL
  z = (z ^ (z >> 27)) * 0x94d049bb133111ebUL
  (s, z ^ (z >> 31))
}

///|
pub fn Xoshiro256::seed(seed : UInt64) -> Xoshiro256 {
  let (a, s0) = splitmix64(seed)
  let (b, s1) = splitmix64(a)
  let (c, s2) = splitmix64(b)
  let (_, s3) = splitmix64(c)
  { s0, s1, s2, s3, }
}

///|
/// Raw state, for checking against the reference implementation.
pub fn Xoshiro256::from_state(
  s0 : UInt64,
  s1 : UInt64,
  s2 : UInt64,
  s3 : UInt64,
) -> Xoshiro256 {
  { s0, s1, s2, s3, }
}

///|
pub impl RandomSource for Xoshiro256 with fn next_u64(self) {
  let result = rotl(self.s1 * 5, 7) * 9
  let t = self.s1 << 17
  self.s2 = self.s2 ^ self.s0
  self.s3 = self.s3 ^ self.s1
  self.s1 = self.s1 ^ self.s2
  self.s0 = self.s0 ^ self.s3
  self.s2 = self.s2 ^ t
  self.s3 = rotl(self.s3, 45)
  result
}

///|
/// Uniform on [0, 1) with 53 random bits.
pub fn uniform(rng : &RandomSource) -> Double {
  (rng.next_u64() >> 11).to_double() * 1.1102230246251565e-16
}

///|
/// Uniform on (0, 1), never 0, for inversion methods that take a log.
pub fn uniform_open(rng : &RandomSource) -> Double {
  ((rng.next_u64() >> 12).to_double() + 0.5) * 2.220446049250313e-16
}

///|
/// Uniform integer in [0, n), n > 0, unbiased: words below 2^64 mod n are
/// rejected so that every residue is equally likely.
pub fn uniform_below(rng : &RandomSource, n : UInt64) -> UInt64 {
  let threshold = (0UL - n) % n
  while true {
    let x = rng.next_u64()
    if x >= threshold {
      return x % n
    }
  }
  0
}

///|
/// Methods for the standard normal. Each is frozen once released.
#derive(Eq, Debug)
pub(all) enum NormalMethod {
  /// Inverse of the normal distribution function on one uniform.
  Inversion
  /// The Marsaglia–Tsang ziggurat.
  Ziggurat
}

///|
pub fn standard_normal(rng : &RandomSource, how : NormalMethod) -> Double {
  ignore((rng, how))
  ...
}

///|
#derive(Eq, Debug)
pub(all) enum ExponentialMethod {
  Inversion
  Ziggurat
}

///|
pub fn standard_exponential(
  rng : &RandomSource,
  how : ExponentialMethod,
) -> Double {
  match how {
    Inversion => -@math.ln(uniform_open(rng))
    Ziggurat => ...
  }
}

///|
#derive(Eq, Debug)
pub(all) enum GammaMethod {
  /// Marsaglia and Tsang (2000), with the shape < 1 boost.
  MarsagliaTsang
}

///|
pub fn standard_gamma(
  rng : &RandomSource,
  shape : Double,
  how : GammaMethod,
  normal~ : NormalMethod,
) -> Double {
  ignore((rng, shape, how, normal))
  ...
}

///|
#derive(Eq, Debug)
pub(all) enum DiscreteMethod {
  /// Inversion by sequential search; exact, slow for large means.
  Inversion
  /// BTPE for the binomial (Kachitvichyanukul and Schmeiser), PTRS for the
  /// Poisson (Hörmann).
  Rejection
}

///|
/// Shuffle in place (Fisher–Yates), the core of permutation tests.
pub fn[T] shuffle(rng : &RandomSource, xs : FixedArray[T]) -> Unit {
  for i = xs.length() - 1; i > 0; i = i - 1 {
    let j = uniform_below(rng, (i + 1).to_uint64()).to_int()
    let t = xs[i]
    xs[i] = xs[j]
    xs[j] = t
  }
}

///|
/// `rng.next_u64()` on a concrete generator, for callers outside this package.
pub extend Xoshiro256 with RandomSource::{next_u64}

///|
/// Independent streams from one seed, for parallel resampling: `jump`
/// advances the state by 2^128 draws (xoshiro256's jump polynomial), so
/// each worker takes the next jumped copy.
pub(open) trait Jump {
  fn jump(Self) -> Unit
}

///|
pub impl Jump for Xoshiro256 with fn jump(self) {
  ignore(self)
  ...
}

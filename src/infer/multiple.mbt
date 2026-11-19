///|
/// Multiple testing, effect sizes and power.

///|
/// p-value adjustments. FWER methods control the chance of any false
/// rejection; FDR methods the expected share of false rejections.
#derive(Eq, Debug)
pub(all) enum Adjust {
  Bonferroni
  /// Holm's step-down (uniformly better than Bonferroni).
  Holm
  /// Hochberg's step-up; needs independence or positive dependence.
  Hochberg
  Hommel
  /// Benjamini–Hochberg FDR; independence or positive dependence.
  BenjaminiHochberg
  /// Benjamini–Yekutieli FDR under any dependence.
  BenjaminiYekutieli
}

///|
/// Adjusted p-values in the input order, as R's p.adjust returns them.
pub fn adjust(p : ArrayView[Double], how : Adjust) -> FixedArray[Double] {
  let n = p.length()
  let out = FixedArray::make(n, 0.0)
  if n == 0 {
    return out
  }
  let fn_ = n.to_double()
  let order = FixedArray::makei(n, i => i)
  order.sort_by((a, b) => p[a].compare(p[b]))
  let clamp = (x : Double) => if x > 1.0 { 1.0 } else { x }
  match how {
    Bonferroni => for i in 0..<n {
      out[i] = clamp(fn_ * p[i])
    }
    Holm => {
      // step-down: running maximum of (n - k) p_(k) over ascending p
      let mut run = 0.0
      for k in 0..<n {
        let v = clamp((fn_ - k.to_double()) * p[order[k]])
        if v > run {
          run = v
        }
        out[order[k]] = run
      }
    }
    Hochberg => {
      // step-up: running minimum of (n - k) p_(k) over descending p
      let mut run = 1.0
      for k = n - 1; k >= 0; k = k - 1 {
        let v = clamp((fn_ - k.to_double()) * p[order[k]])
        if v < run {
          run = v
        }
        out[order[k]] = run
      }
    }
    BenjaminiHochberg | BenjaminiYekutieli => {
      let c = match how {
        BenjaminiYekutieli => {
          let mut h = 0.0
          for i in 1..=n {
            h = h + 1.0 / i.to_double()
          }
          h
        }
        _ => 1.0
      }
      let mut run = 1.0
      for k = n - 1; k >= 0; k = k - 1 {
        let v = clamp(c * fn_ / (k + 1).to_double() * p[order[k]])
        if v < run {
          run = v
        }
        out[order[k]] = run
      }
    }
    Hommel => ...
  }
  out
}

///|
/// Which hypotheses are rejected at level `alpha` after adjustment.
pub fn reject(
  p : ArrayView[Double],
  alpha : Double,
  how : Adjust,
) -> FixedArray[Bool] {
  let adj = adjust(p, how)
  FixedArray::makei(adj.length(), i => adj[i] <= alpha)
}

///|
#derive(Eq, Debug)
pub(all) enum EffectKind {
  /// Standardised mean difference with the pooled SD.
  CohenD
  /// Cohen's d with the small-sample correction.
  HedgesG
  /// Mean difference over the control group's SD.
  GlassDelta
  /// Correlation as an effect size.
  R
  EtaSquared
  OmegaSquared
  CramersV
  OddsRatio
  RiskRatio
  /// For rank tests: P(X > Y) - P(X < Y).
  RankBiserial
}

///|
pub struct Effect {
  kind : EffectKind
  value : Double
  interval : Interval?
}

///|
/// Power analysis solves for exactly one unknown: give three of sample
/// size, effect, significance level and power, leave the fourth `None`.
pub struct PowerQuery {
  n : Double?
  effect : Double?
  alpha : Double?
  power : Double?
}

///|
#derive(Debug)
pub(all) suberror PowerError {
  /// Not exactly one of the four quantities was left open.
  UnknownCount(Int)
  /// No value of the unknown achieves the others.
  Unreachable
}

///|
pub fn PowerQuery::new(
  n? : Double,
  effect? : Double,
  alpha? : Double,
  power? : Double,
) -> PowerQuery raise PowerError {
  let open = [n, effect, alpha, power].iter().filter(v => v is None).count()
  if open != 1 {
    raise UnknownCount(open)
  }
  { n, effect, alpha, power }
}

///|
#derive(Eq, Debug)
pub(all) enum PowerDesign {
  OneSampleT
  TwoSampleT
  PairedT
  OneProportion
  TwoProportions
  OneWayAnova(groups~ : Int)
  ChiSquare(df~ : Int)
}

///|
/// The solved value of the open quantity. Uses the noncentral t, χ² or F
/// distribution, not a normal approximation.
pub fn solve_power(
  q : PowerQuery,
  design : PowerDesign,
  alternative~ : Alternative,
) -> Double raise PowerError {
  ignore((q, design, alternative))
  ...
}

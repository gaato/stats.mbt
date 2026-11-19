///|
/// Order statistics: quantiles, ranks, and the robust summaries built on
/// them. Quantiles and medians are the most reimplemented primitive on
/// mooncakes (31 and 23 modules), always with a hard-coded definition;
/// here the definition is an argument.

///|
/// The p-quantile of data already sorted ascending, by one of the nine
/// Hyndman–Fan definitions. Index arithmetic follows the paper with
/// 1-based positions; positions outside the data are clamped to the ends.
/// NaN for empty data and for p outside [0, 1]: this is the unchecked
/// inner step, and `quantile` raises for the same inputs.
pub fn quantile_sorted(
  sorted : ArrayView[Double],
  p : Double,
  how : @data.QuantileMethod,
) -> Double {
  let n = sorted.length()
  if n == 0 || !(p >= 0.0 && p <= 1.0) {
    return 0.0 / 0.0
  }
  // x(k) for 1-based k, clamped to [1, n]
  let at = (k : Int32) => {
    sorted[(if k < 1 { 1 } else if k > n { n } else { k }) - 1]
  }
  let np = n.to_double() * p
  match how {
    Type1 => {
      let j = @math.floor(np).to_int()
      if np - j.to_double() > 0.0 {
        at(j + 1)
      } else {
        at(j)
      }
    }
    Type2 => {
      let j = @math.floor(np).to_int()
      if np - j.to_double() > 0.0 {
        at(j + 1)
      } else {
        0.5 * (at(j) + at(j + 1))
      }
    }
    Type3 => {
      let h = np - 0.5
      let j = @math.floor(h).to_int()
      if h - j.to_double() == 0.0 && j % 2 == 0 {
        at(j)
      } else {
        at(j + 1)
      }
    }
    _ => {
      // continuous types: position h = np + m with the type's m
      let m = match how {
        Type4 => 0.0
        Type5 => 0.5
        Type6 => p
        Type7 => 1.0 - p
        Type8 => (p + 1.0) / 3.0
        _ => p / 4.0 + 3.0 / 8.0
      }
      let h = np + m
      let j = @math.floor(h).to_int()
      let g = h - j.to_double()
      if g == 0.0 {
        at(j)
      } else {
        (1.0 - g) * at(j) + g * at(j + 1)
      }
    }
  }
}

///|
/// Sorts a copy, then `quantile_sorted`. Rejects empty input, NaN, and a
/// p outside [0, 1].
pub fn quantile(
  data : ArrayView[Double],
  p : Double,
  how : @data.QuantileMethod,
) -> Double raise @data.DataError {
  if !(p >= 0.0 && p <= 1.0) {
    raise @data.InvalidProbability(p)
  }
  let (xs, _) = @data.clean(data, Reject)
  if xs.length() == 0 {
    raise @data.Empty
  }
  xs.sort()
  quantile_sorted(xs[:], p, how)
}

///|
/// Ranks starting at 1, with the tie rule given by the caller.
pub fn ranks(
  data : ArrayView[Double],
  ties : @data.RankTies,
) -> FixedArray[Double] {
  let n = data.length()
  let order = FixedArray::makei(n, i => i)
  order.sort_by((a, b) => {
    let c = data[a].compare(data[b])
    if c != 0 {
      c
    } else {
      a.compare(b)
    }
  })
  let out = FixedArray::make(n, 0.0)
  let mut i = 0
  let mut dense = 0
  while i < n {
    let mut j = i
    while j + 1 < n && data[order[j + 1]] == data[order[i]] {
      j = j + 1
    }
    dense = dense + 1
    for k in i..<=j {
      out[order[k]] = match ties {
        Average => 0.5 * ((i + 1).to_double() + (j + 1).to_double())
        Min => (i + 1).to_double()
        Max => (j + 1).to_double()
        Dense => dense.to_double()
        Ordinal => (k + 1).to_double()
      }
    }
    i = j + 1
  }
  out
}

///|
/// The scale factor of the median absolute deviation: `Raw` is the plain
/// median of |x - median|; `Normal` multiplies by 1.4826 so that it
/// estimates σ for normal data (R's default, not SciPy's).
#derive(Eq, Debug)
pub(all) enum MadScale {
  Raw
  Normal
}

///|
pub fn mad(
  data : ArrayView[Double],
  scale : MadScale,
) -> Double raise @data.DataError {
  let m = quantile(data, 0.5, Type7)
  let dev = FixedArray::from_iter(data.iter().map(x => (x - m).abs()))
  let raw = quantile(dev[:], 0.5, Type7)
  match scale {
    Raw => raw
    // 1/Φ⁻¹(3/4), so that the MAD estimates σ for normal data
    Normal => 1.482602218505602 * raw
  }
}

///|
/// Mean of the data with `proportion` cut from each end.
/// floor(proportion · n) values are cut from each end, as SciPy's trim_mean.
/// The proportion is in [0, 1/2): from a half on nothing would be left.
pub fn trimmed_mean(
  data : ArrayView[Double],
  proportion : Double,
) -> Double raise @data.DataError {
  if !(proportion >= 0.0 && proportion < 0.5) {
    raise @data.InvalidProbability(proportion)
  }
  let (xs, _) = @data.clean(data, Reject)
  if xs.length() == 0 {
    raise @data.Empty
  }
  xs.sort()
  let cut = @math.floor(proportion * xs.length().to_double()).to_int()
  let mut s = @num.CompensatedSum::zero()
  for i in cut..<(xs.length() - cut) {
    s = s.add(xs[i])
  }
  s.value() / (xs.length() - 2 * cut).to_double()
}

///|
/// Values beyond the given quantiles replaced by those quantiles.
pub fn winsorize(
  data : ArrayView[Double],
  lower : Double,
  upper : Double,
  how : @data.QuantileMethod,
) -> FixedArray[Double] raise @data.DataError {
  ignore((data, lower, upper, how))
  ...
}

///|
/// How histogram bins are chosen.
pub(all) enum Bins {
  /// Explicit edges, ascending; values on an inner edge go to the right
  /// bin, the last edge is closed.
  Edges(FixedArray[Double])
  /// k equal-width bins over [lo, hi].
  Uniform(k~ : Int32, lo~ : Double, hi~ : Double)
  Sturges
  Scott
  FreedmanDiaconis
}

///|
pub struct Histogram {
  edges : FixedArray[Double]
  counts : FixedArray[Int64]
  /// Values outside the edges, counted rather than dropped silently.
  below : Int64
  above : Int64
}

///|
pub fn histogram(
  data : ArrayView[Double],
  bins : Bins,
) -> Histogram raise @data.DataError {
  ignore((data, bins))
  ...
}

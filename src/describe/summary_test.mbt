///|
test "plain summaries of an array" {
  let xs : Array[Double] = [2.0, 4.0, 4.0, 4.0, 5.0, 5.0, 7.0, 9.0]
  assert_eq(@describe.sum(xs), 40.0)
  assert_eq(@describe.mean(xs), 5.0)
  assert_eq(@describe.variance(xs, ddof=0), 4.0)
  assert_eq(@describe.variance(xs, ddof=1), 32.0 / 7.0)
  assert_eq(@describe.std_dev(xs, ddof=0), 2.0)
  assert_eq(@describe.median(xs), 4.5)
  assert_eq(@describe.median([3.0, 1.0, 2.0]), 2.0)
  // the sum keeps what plain addition loses
  assert_eq(@describe.sum([1.0e16, 1.0, -1.0e16, 1.0]), 2.0)
}

///|
test "summaries refuse what they cannot summarise" {
  let empty : Array[Double] = []
  let is_empty = try @describe.mean(empty) catch {
    @data.Empty => true
    _ => false
  } noraise {
    _ => false
  }
  assert_true(is_empty)
  let nan = try @describe.sum([1.0, 0.0 / 0.0]) catch {
    @data.ContainsNaN(index~) => index == 1
    _ => false
  } noraise {
    _ => false
  }
  assert_true(nan)
  let few = try @describe.variance([1.0], ddof=1) catch {
    @data.TooFewValues(needed~, got~) => needed == 2 && got == 1
    _ => false
  } noraise {
    _ => false
  }
  assert_true(few)
}

///|
/// A probability outside [0, 1] used to return an end of the data.
test "quantile and trimmed mean check their proportion" {
  let xs : Array[Double] = [1.0, 2.0, 3.0, 4.0]
  for p in [2.0, -1.0, 0.0 / 0.0] {
    let bad = try @describe.quantile(xs, p, Type7) catch {
      @data.InvalidProbability(_) => true
      _ => false
    } noraise {
      _ => false
    }
    assert_true(bad)
    assert_true(@describe.quantile_sorted(xs, p, Type1).is_nan())
  }
  assert_true(@describe.quantile_sorted([], 0.5, Type7).is_nan())
  assert_eq(@describe.quantile(xs, 0.0, Type7), 1.0)
  assert_eq(@describe.quantile(xs, 1.0, Type7), 4.0)
  let half = try @describe.trimmed_mean(xs, 0.5) catch {
    @data.InvalidProbability(_) => true
    _ => false
  } noraise {
    _ => false
  }
  assert_true(half)
  assert_eq(@describe.trimmed_mean(xs, 0.25), 2.5)
}

///|
test "accumulators start, merge and fold the same way" {
  let xs : Array[Double] = [3.0, -1.0, 4.0, 1.0, 5.0]
  let mv = @describe.MeanVar::empty().push(3.0).push(-1.0)
  assert_eq(mv.count(), 2L)
  assert_eq(mv.mean(), Some(1.0))
  let all : @describe.MeanVar = @describe.accumulate(xs)
  assert_true((all.mean().unwrap() - 2.4).abs() < 1.0e-15)
  let e : @describe.Extremes = @describe.accumulate(xs)
  assert_eq((e.count(), e.min(), e.max()), (5L, Some(-1.0), Some(5.0)))
  let none = @describe.Extremes::empty()
  assert_eq((none.count(), none.min(), none.max()), (0L, None, None))
  let left : @describe.Extremes = @describe.accumulate(xs[0:2])
  let right : @describe.Extremes = @describe.accumulate(xs[2:])
  let both = left.merge(right)
  assert_eq((both.count(), both.min(), both.max()), (5L, Some(-1.0), Some(5.0)))
}

///|
test "co-moments report both variables" {
  let mut c = @describe.Comoment::empty()
  let xs = [1.0, 2.0, 3.0, 4.0]
  let ys = [2.0, 4.0, 7.0, 7.0]
  for i in 0..<4 {
    c = c.push(xs[i], ys[i])
  }
  assert_eq(c.count(), 4L)
  assert_eq(c.mean_x(), Some(2.5))
  assert_eq(c.mean_y(), Some(5.0))
  assert_eq(c.variance_x(ddof=1), Some(5.0 / 3.0))
  assert_eq(c.variance_y(ddof=1), Some(6.0))
  assert_eq(c.variance_y(ddof=4), None)
}

///|
test "weighted means merge" {
  let a = @describe.WeightedMeanVar::empty(Frequency)
    .push(1.0, 2.0)
    .push(4.0, 1.0)
  let b = @describe.WeightedMeanVar::empty(Frequency).push(10.0, 3.0)
  let whole = @describe.WeightedMeanVar::empty(Frequency)
    .push(1.0, 2.0)
    .push(4.0, 1.0)
    .push(10.0, 3.0)
  let merged = a.merge(b)
  assert_eq(merged.sum_weights(), 6.0)
  assert_true((merged.mean().unwrap() - whole.mean().unwrap()).abs() < 1.0e-15)
  assert_eq(merged.mean(), Some(6.0))
  let none = @describe.WeightedMeanVar::empty(Frequency)
  assert_eq(none.merge(b).mean(), Some(10.0))
  assert_eq(b.merge(none).mean(), Some(10.0))
}

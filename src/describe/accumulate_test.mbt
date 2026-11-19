///|
/// The monoid laws, checked on data whose two halves have very different
/// means: merging the halves must agree with one pass over everything.
test "merge agrees with a single pass" {
  let xs = [1.0e9 + 4.0, 1.0e9 + 7.0, 1.0e9 + 13.0, 1.0e9 + 16.0]
  let ys = [2.0, 3.0, 5.0, 7.0, 11.0]
  let empty : MeanVar = Accumulator::empty()
  let one_pass = (xs + ys).fold(init=empty, (acc, x) => acc.push(x))
  let left = xs.fold(init=empty, (acc, x) => acc.push(x))
  let right = ys.fold(init=empty, (acc, x) => acc.push(x))
  let merged = left.merge(right)
  assert_eq(merged.count(), one_pass.count())
  let rel = (a : Double, b : Double) => (a - b).abs() / b.abs()
  assert_true(rel(merged.mean().unwrap(), one_pass.mean().unwrap()) < 1.0e-15)
  assert_true(
    rel(merged.variance(ddof=1).unwrap(), one_pass.variance(ddof=1).unwrap()) <
    1.0e-12,
  )
  // Identity on both sides.
  assert_eq(empty.merge(left).mean(), left.mean())
  assert_eq(left.merge(empty).mean(), left.mean())
  // Undefined moments are None, not NaN.
  assert_eq(empty.mean(), None)
  assert_eq(empty.push(1.0).variance(ddof=1), None)
}

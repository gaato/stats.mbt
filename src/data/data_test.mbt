///|
test "missing-value policies" {
  let xs = [1.0, 0.0 / 0.0, 3.0]
  let (kept, dropped) = @data.clean(xs[:], Omit)
  assert_eq(kept, ([1.0, 3.0] : FixedArray[Double]))
  assert_eq(dropped, 1)
  let rejected = try @data.clean(xs[:], Reject) catch {
    ContainsNaN(index~) => index == 1
    _ => false
  } noraise {
    _ => false
  }
  assert_true(rejected)
  let (all, _) = @data.clean(xs[:], Propagate)
  assert_true(all[1].is_nan())
}

///|
test "contingency table" {
  let t = @data.Table::new(2, 2, [3L, 1L, 1L, 3L][:])
  assert_eq(t.at(0, 1), 1L)
  assert_eq(t.total(), 8L)
  let bad = try @data.Table::new(2, 2, [1L, 2L, 3L][:]) catch {
    LengthMismatch(expected~, got~) => expected == 4 && got == 3
    _ => false
  } noraise {
    _ => false
  }
  assert_true(bad)
}

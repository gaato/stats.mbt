///|
test "a design checks its per-row inputs against the matrix" {
  let x = @data.DataMatrix::from_rows([[1.0, 2.0], [3.0, 4.0], [5.0, 6.0]])
  let d = @regress.Design::new(x, intercept=true)
  assert_eq(d.names, ["x0", "x1"])
  assert_eq(d.offset, None)
  let named = @regress.Design::new(x, intercept=false, names=["age", "dose"])
  assert_eq(named.names, ["age", "dose"])
  let w = @data.Weights::new([1.0, 1.0], kind=Frequency)
  let short = try @regress.Design::new(x, intercept=true, weights=w) catch {
    LengthMismatch(expected~, got~) => expected == 3 && got == 2
    _ => false
  } noraise {
    _ => false
  }
  assert_true(short)
  let offset = try @regress.Design::new(x, intercept=true, offset=[0.0]) catch {
    LengthMismatch(expected~, got~) => expected == 3 && got == 1
    _ => false
  } noraise {
    _ => false
  }
  assert_true(offset)
}

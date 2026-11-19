///|
test "draws made elsewhere become a chain" {
  let c = @bayes.Chain::from_draws(@data.DataMatrix::new(2, 1, [0.5, 0.7]))
  assert_eq(c.acceptance, None)
  assert_eq(c.draws.rows, 2)
}

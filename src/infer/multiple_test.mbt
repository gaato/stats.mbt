///|
/// Expected values from statsmodels' multipletests for the same p-values.
test "p-value adjustments agree with statsmodels" {
  let p = [0.01, 0.04, 0.03, 0.005, 0.2, 0.04]
  let cases : Array[(@infer.Adjust, Array[Double])] = [
    (Bonferroni, [0.06, 0.24, 0.18, 0.03, 1.0, 0.24]),
    (Holm, [0.05, 0.12, 0.12, 0.03, 0.2, 0.12]),
    (Hochberg, [0.05, 0.08, 0.08, 0.03, 0.2, 0.08]),
    (BenjaminiHochberg, [0.03, 0.048, 0.048, 0.03, 0.2, 0.048]),
    (BenjaminiYekutieli, [0.0735, 0.1176, 0.1176, 0.0735, 0.49, 0.1176]),
  ]
  for c in cases {
    let (how, want) = c
    let got = @infer.adjust(p[:], how)
    for i in 0..<p.length() {
      assert_true(
        (got[i] - want[i]).abs() <= 1.0e-12,
        msg="case \{i}: got \{got[i]}, want \{want[i]}",
      )
    }
  }
  assert_eq(
    @infer.reject(p[:], alpha=0.05, BenjaminiHochberg),
    FixedArray::from_array([true, true, true, true, false, true]),
  )
}

///|
/// The unknown is the case: a query names one and supplies the other three.
test "power queries name their one unknown" {
  let q = @infer.PowerQuery::SampleSize(effect=0.5, alpha=0.05, power=0.8)
  guard q is SampleSize(power~, ..) else {
    fail("expected a sample-size query")
  }
  assert_eq(power, 0.8)
}

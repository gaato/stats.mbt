///|
/// Shape parameters far beyond the generated test points. Next to x = a
/// the series and the continued fraction need about 9√a terms; they are
/// allowed that many, and where even that is out of reach the result is
/// NaN, never a partial sum. References are 50-digit mpmath values.
test "incomplete gamma with a large shape parameter" {
  let close = (got : Double, want : Double, tol : Double) => {
    (got - want).abs() <= tol * want.abs()
  }
  // 9e5 terms of the series: the sum itself costs a few 1e-13
  assert_true(
    close(
      @special.gammaincc(1.0e10 + 1.0, 1.0e10),
      0.50000265961520264,
      2.0e-12,
    ),
  )
  assert_true(
    close(@special.gammaincc(1.0e8, 1.0e8), 0.49998670192398588, 4.0e-13),
  )
  // ten standard deviations out, by the continued fraction
  assert_true(
    close(
      @special.gammaincc(1.0e8, 1.0005e8),
      2.8784296868527810513e-7,
      1.0e-13,
    ),
  )
  assert_true(
    close(
      @special.gammaincc(1.0e10 - 1.0e6 + 1.0, 1.0e10),
      7.607550991662096336e-24,
      1.0e-13,
    ),
  )
  // beyond the term limit: not computed, and said so
  assert_true(@special.gammainc(1.0e14, 1.0e14).is_nan())
  assert_true(@special.gammaincc(1.0e14, 1.0e14).is_nan())
}

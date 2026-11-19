///|
/// The published reference outputs of xoshiro256** from state {1, 2, 3, 4},
/// and splitmix64 seeding checked against an independent implementation.
test "xoshiro256** matches the reference" {
  let r = @num.Xoshiro256::from_state(1, 2, 3, 4)
  let got = [r.next_u64(), r.next_u64(), r.next_u64(), r.next_u64()]
  assert_eq(got, [11520UL, 0UL, 1509978240UL, 1215971899390074240UL])
  let s = @num.Xoshiro256::seed(42)
  assert_eq([s.next_u64(), s.next_u64(), s.next_u64()], [
    1546998764402558742UL, 6990951692964543102UL, 12544586762248559009UL,
  ])
}

///|
test "uniform conversions stay in range" {
  let r = @num.Xoshiro256::seed(7)
  for _ in 0..<10000 {
    let u = @num.uniform(r)
    assert_true(u >= 0.0 && u < 1.0)
    let v = @num.uniform_open(r)
    assert_true(v > 0.0 && v < 1.0)
  }
  // the extreme words map to the ends of the ranges
  let lo = @num.Xoshiro256::from_state(0, 0, 0, 0)
  assert_eq(@num.uniform(lo), 0.0)
}

///|
test "log-space arithmetic survives where the naive forms do not" {
  // exp(1000) overflows, the stable sum does not
  assert_eq(@num.log_add_exp(1000.0, 1000.0), 1000.0 + @math.ln(2.0))
  assert_eq(@num.log_add_exp(-1.0 / 0.0, -3.0), -3.0)
  // log(1 - exp(-1e-20)) = log(1e-20) to full precision
  assert_true((@num.log1m_exp(-1.0e-20) - @math.ln(1.0e-20)).abs() < 1.0e-12)
  assert_eq(@num.log1p_exp(800.0), 800.0)
  assert_eq(@num.log_sum_exp([-1000.0, -1000.0][:]), -1000.0 + @math.ln(2.0))
  assert_eq(@num.log1p(1.0e-300), 1.0e-300)
}

///|
test "compensated sum keeps what plain summation loses" {
  let xs = [1.0e16, 1.0, -1.0e16, 1.0]
  let mut plain = 0.0
  let mut acc = @num.CompensatedSum::zero()
  for x in xs {
    plain = plain + x
    acc = acc.add(x)
  }
  assert_eq(plain, 1.0)
  assert_eq(acc.value(), 2.0)
}

///|
test "brent and bracketed newton find x = cos x" {
  let r = @num.find_root(x => x - @math.cos(x), 0.0, 1.0)
  assert_eq(r.status, Converged)
  assert_true((r.value - 0.7390851332151607).abs() < 2.0e-16)
  let n = @num.find_root_newton(x => (x - @math.cos(x), 1.0 + @math.sin(x)), 0.0, 1.0, 0.5)
  assert_eq(n.status, Converged)
  assert_true((n.value - 0.7390851332151607).abs() < 2.0e-16)
  assert_eq(@num.find_root(x => x * x + 1.0, -1.0, 1.0).status, NoBracket)
}

///|
test "series and continued fraction" {
  // e = Σ 1/n!
  let s = @num.sum_series(n => {
    let mut f = 1.0
    for i in 1..=n {
      f = f * i.to_double()
    }
    1.0 / f
  })
  assert_true((s.value - 2.718281828459045).abs() < 5.0e-16)
  // √2 = 1 + 1/(2 + 1/(2 + ...))
  let cf = @num.continued_fraction(1.0, _ => (1.0, 2.0))
  assert_eq(cf.status, Converged)
  // Lentz multiplies one rounded ratio per term, so the error grows with the
  // number of terms: 22 terms here, about 5 ulp
  assert_true((cf.value - 1.4142135623730951).abs() < 1.0e-15, msg="\{cf.value} after \{cf.iterations}")
}

///|
test "tanh-sinh handles endpoint singularities" {
  let a = @num.integrate(x => @math.ln(x), 0.0, 1.0)
  assert_true((a.value + 1.0).abs() < 1.0e-14, msg="\{a.value}")
  let b = @num.integrate(x => x.sqrt(), 0.0, 1.0)
  assert_true((b.value - 2.0 / 3.0).abs() < 1.0e-15, msg="\{b.value}")
  let c = @num.integrate(x => 1.0 / (1.0 + x * x), 0.0, 1.0)
  assert_true((c.value - 0.25 * @math.PI).abs() < 4.0e-16, msg="\{c.value}")
}

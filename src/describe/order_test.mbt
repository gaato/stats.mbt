///|
/// Expected values from NumPy's quantile, whose nine methods are the
/// Hyndman–Fan types in order.
test "nine quantile definitions agree with NumPy" {
  let xs = [3.1, 1.2, 7.7, 4.4, 9.0, 2.5, 6.3, 5.8, 8.1, 0.4, 4.4]
  let methods : Array[@data.QuantileMethod] = [
    Type1,
    Type2,
    Type3,
    Type4,
    Type5,
    Type6,
    Type7,
    Type8,
    Type9,
  ]
  let expected : Array[(Double, Array[Double])] = [
    (0.0, [0.4, 0.4, 0.4, 0.4, 0.4, 0.4, 0.4, 0.4, 0.4]),
    (
      0.1,
      [
        1.2, 1.2, 0.4, 0.4800000000000001, 0.8800000000000001, 0.5600000000000002,
        1.2, 0.7733333333333334, 0.8,
      ],
    ),
    (0.25, [2.5, 2.5, 2.5, 2.175, 2.65, 2.5, 2.8, 2.6, 2.6125]),
    (0.5, [4.4, 4.4, 4.4, 4.4, 4.4, 4.4, 4.4, 4.4, 4.4]),
    (
      0.73,
      [
        7.7, 7.7, 6.3, 6.341999999999999, 7.041999999999999, 7.364, 6.72, 7.149333333333332,
        7.122499999999998,
      ],
    ),
    (0.99, [9.0, 9.0, 9.0, 8.901, 9.0, 9.0, 8.91, 9.0, 9.0]),
    (1.0, [9.0, 9.0, 9.0, 9.0, 9.0, 9.0, 9.0, 9.0, 9.0]),
  ]
  for row in expected {
    let (p, want) = row
    for k, m in methods {
      let got = @describe.quantile(xs[:], p, m)
      assert_true(
        (got - want[k]).abs() <= 4.0e-15 * want[k].abs(),
        msg="p=\{p} type \{k + 1}: got \{got}, want \{want[k]}",
      )
    }
  }
}

///|
/// Expected values from SciPy's rankdata.
test "rank tie rules agree with SciPy" {
  let xs = [3.1, 1.2, 7.7, 4.4, 9.0, 2.5, 6.3, 5.8, 8.1, 0.4, 4.4, 1.2, 4.4]
  let cases : Array[(@data.RankTies, Array[Double])] = [
    (Average, [5, 2.5, 11, 7, 13, 4, 10, 9, 12, 1, 7, 2.5, 7]),
    (Min, [5, 2, 11, 6, 13, 4, 10, 9, 12, 1, 6, 2, 6]),
    (Max, [5, 3, 11, 8, 13, 4, 10, 9, 12, 1, 8, 3, 8]),
    (Dense, [4, 2, 8, 5, 10, 3, 7, 6, 9, 1, 5, 2, 5]),
    (Ordinal, [5, 2, 11, 6, 13, 4, 10, 9, 12, 1, 7, 3, 8]),
  ]
  for c in cases {
    let (ties, want) = c
    assert_eq(@describe.ranks(xs[:], ties), FixedArray::from_array(want))
  }
}

///|
test "comoment merges and gives pearson correlation" {
  let x = [1.0e9 + 1.0, 1.0e9 + 2.0, 1.0e9 + 4.0, 1.0e9 + 7.0, 1.0e9 + 8.0]
  let y = [2.0, 4.1, 7.9, 14.2, 15.8]
  let mut one = @describe.Comoment::empty()
  for i in 0..<5 {
    one = one.push(x[i], y[i])
  }
  let mut a = @describe.Comoment::empty()
  let mut b = @describe.Comoment::empty()
  for i in 0..<2 {
    a = a.push(x[i], y[i])
  }
  for i in 2..<5 {
    b = b.push(x[i], y[i])
  }
  let merged = a.merge(b)
  let r1 = one.correlation().unwrap()
  let r2 = merged.correlation().unwrap()
  assert_true((r1 - r2).abs() < 1.0e-14)
  // the same correlation as for the data without the offset
  let mut c = @describe.Comoment::empty()
  for i in 0..<5 {
    c = c.push(x[i] - 1.0e9, y[i])
  }
  assert_true((r1 - c.correlation().unwrap()).abs() < 1.0e-12)
}

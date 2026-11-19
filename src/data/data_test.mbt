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

///|
test "a table checks its shape and each index" {
  let shape = try @data.Table::new(-2, -3, [1L, 2L, 3L, 4L, 5L, 6L]) catch {
    InvalidShape(rows~, cols~) => rows == -2 && cols == -3
    _ => false
  } noraise {
    _ => false
  }
  assert_true(shape)
  let t = @data.Table::new(2, 3, [1L, 2L, 3L, 4L, 5L, 6L])
  assert_eq(t.row_totals(), [6L, 15L])
  assert_eq(t.col_totals(), [5L, 7L, 9L])
}

///|
test "panic a column past the end does not read the next row" {
  let t = @data.Table::new(2, 3, [1L, 2L, 3L, 4L, 5L, 6L])
  ignore(t.at(0, 4))
}

///|
test "cross-tabulation from two label columns" {
  let t = @data.Table::from_labels(
    [0, 0, 1, 1, 1],
    [0, 1, 1, 1, 0],
    rows=2,
    cols=2,
  )
  assert_eq(t.counts, [1L, 1L, 1L, 2L])
  let bad = try @data.Table::from_labels([0, 2], [0, 0], rows=2, cols=2) catch {
    InvalidLabel(index~, label~) => index == 1 && label == 2
    _ => false
  } noraise {
    _ => false
  }
  assert_true(bad)
}

///|
test "weights carry their kind and reject negatives" {
  let w = @data.Weights::new([1.0, 2.0, 0.0], kind=Frequency)
  let xs = @data.Weighted::new([5.0, 6.0, 7.0], w)
  assert_eq(xs.weights.kind, Frequency)
  assert_eq(xs.values, [5.0, 6.0, 7.0])
  let neg = try @data.Weights::new([1.0, -1.0], kind=Reliability) catch {
    NegativeWeight(index~) => index == 1
    _ => false
  } noraise {
    _ => false
  }
  assert_true(neg)
  let nan = try @data.Weights::new([0.0 / 0.0], kind=Reliability) catch {
    NegativeWeight(index~) => index == 0
    _ => false
  } noraise {
    _ => false
  }
  assert_true(nan)
  let short = try @data.Weighted::new([1.0], w) catch {
    LengthMismatch(expected~, got~) => expected == 1 && got == 3
    _ => false
  } noraise {
    _ => false
  }
  assert_true(short)
}

///|
test "pairs of unequal length are refused" {
  let bad = try @data.Bivariate::new([1.0, 2.0, 3.0], [1.0]) catch {
    LengthMismatch(expected~, got~) => expected == 3 && got == 1
    _ => false
  } noraise {
    _ => false
  }
  assert_true(bad)
  assert_eq(@data.Bivariate::new([1.0, 2.0], [3.0, 4.0]).length(), 2)
}

///|
test "grouped values come back by label" {
  let g = @data.Grouped::new(
    [1.0, 2.0, 3.0, 4.0, 5.0],
    [0, 2, 0, 2, 2],
    groups=3,
  )
  assert_eq(g.group(0), [1.0, 3.0])
  assert_eq(g.group(1), [])
  assert_eq(g.group(2), [2.0, 4.0, 5.0])
  let bad = try @data.Grouped::new([1.0], [3], groups=3) catch {
    InvalidLabel(index~, label~) => index == 0 && label == 3
    _ => false
  } noraise {
    _ => false
  }
  assert_true(bad)
}

///|
test "blocks and matrices are row-major and checked" {
  let b = @data.Blocks::new(2, 3, [1.0, 2.0, 3.0, 4.0, 5.0, 6.0])
  assert_eq(b.at(1, 0), 4.0)
  let m = @data.DataMatrix::from_rows([[1.0, 2.0], [3.0, 4.0], [5.0, 6.0]])
  assert_eq(m.rows, 3)
  assert_eq(m.at(2, 1), 6.0)
  assert_eq(m.row(1).to_owned(), [3.0, 4.0])
  assert_eq(m.column(1), [2.0, 4.0, 6.0])
  let ragged = try @data.DataMatrix::from_rows([[1.0, 2.0], [3.0]]) catch {
    LengthMismatch(expected~, got~) => expected == 2 && got == 1
    _ => false
  } noraise {
    _ => false
  }
  assert_true(ragged)
  let empty = try @data.DataMatrix::new(0, 2, []) catch {
    InvalidShape(..) => true
    _ => false
  } noraise {
    _ => false
  }
  assert_true(empty)
}

///|
test "a survival record starts at time 0 unless told otherwise" {
  let r = @data.SurvivalRecord::new(3.5, event=false)
  assert_eq(r.entry, 0.0)
  assert_eq(r.to_observation(), RightCensored(3.5))
  assert_eq(@data.SurvivalRecord::new(3.5, event=true, entry=1.0).entry, 1.0)
}

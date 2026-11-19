///|
/// Multivariate methods on a data matrix (observations × variables).
/// Every method here is an eigenvalue or SVD problem, so the bodies wait
/// for dense linear algebra. Failures are `@regress.ModelError`.

///|
#non_exhaustive
#derive(Eq, Debug)
pub(all) enum Scaling {
  /// Center only (covariance-based).
  Center
  /// Center and scale to unit variance (correlation-based).
  Standardize
}

///|
/// Principal components by SVD of the centred data, not by the
/// eigen-decomposition of the covariance matrix.
#derive(Eq, Debug)
pub struct Pca {
  loadings : @data.DataMatrix
  scores : @data.DataMatrix
  variances : FixedArray[Double]
  explained : FixedArray[Double]
}

///|
pub fn pca(
  x : @data.DataMatrix,
  scaling : Scaling,
) -> Pca raise @regress.ModelError {
  ignore((x, scaling))
  ...
}

///|
#non_exhaustive
#derive(Eq, Debug)
pub(all) enum Rotation {
  Unrotated
  Varimax
  Promax(power~ : Double)
}

///|
#derive(Eq, Debug)
pub struct FactorFit {
  loadings : @data.DataMatrix
  uniqueness : FixedArray[Double]
  loglik : Double
  iterations : Int32
  converged : Bool
}

///|
/// Maximum-likelihood factor analysis.
pub fn factor_analysis(
  x : @data.DataMatrix,
  factors~ : Int32,
  rotation : Rotation,
) -> FactorFit raise @regress.ModelError {
  ignore((x, factors, rotation))
  ...
}

///|
/// Linear (shared covariance) or quadratic discriminant analysis.
#non_exhaustive
#derive(Eq, Debug)
pub(all) enum Discriminant {
  Linear
  Quadratic
}

///|
/// A fitted classifier.
pub struct Classifier {
  classes : Int32
  /// Class probabilities for new rows, one row per input row.
  predict : (@data.DataMatrix) -> @data.DataMatrix raise @regress.ModelError
}

///|
/// `labels` has one class in `0 ..< classes` per row of `x`.
pub fn discriminant(
  x : @data.DataMatrix,
  labels : ArrayView[Int32],
  kind : Discriminant,
) -> Classifier raise @regress.ModelError {
  if labels.length() != x.rows {
    raise LengthMismatch(expected=x.rows, got=labels.length())
  }
  ignore(kind)
  ...
}

///|
/// ROC curve and the area under it for scores against binary labels.
#derive(Eq, Debug)
pub struct Roc {
  false_positive : FixedArray[Double]
  true_positive : FixedArray[Double]
  thresholds : FixedArray[Double]
  auc : Double
}

///|
pub fn roc(
  scores : ArrayView[Double],
  positive : ArrayView[Bool],
) -> Roc raise @regress.ModelError {
  if scores.length() != positive.length() {
    raise LengthMismatch(expected=scores.length(), got=positive.length())
  }
  ...
}

///|
#derive(Eq, Debug)
pub struct KMeans {
  centers : @data.DataMatrix
  labels : FixedArray[Int32]
  within_ss : FixedArray[Double]
  iterations : Int32
  converged : Bool
}

///|
/// k-means with k-means++ seeding.
pub fn kmeans(
  x : @data.DataMatrix,
  k~ : Int32,
  rng~ : &@num.RandomSource,
) -> KMeans raise @regress.ModelError {
  ignore((x, k, rng))
  ...
}

///|
#non_exhaustive
#derive(Eq, Debug)
pub(all) enum Linkage {
  Single
  Complete
  Average
  Ward
}

///|
/// One step of agglomerative clustering: the two clusters joined and the
/// height at which they were. Clusters 0 ..< n are the observations; the
/// one made by merge i is n + i.
#derive(Eq, Debug)
pub struct Merge {
  first : Int32
  second : Int32
  height : Double
}

///|
/// Agglomerative clustering: the merges in order.
#derive(Eq, Debug)
pub struct Dendrogram {
  merges : Array[Merge]
}

///|
pub fn hierarchical(
  distances : @data.DataMatrix,
  linkage : Linkage,
) -> Dendrogram raise @regress.ModelError {
  ignore((distances, linkage))
  ...
}

///|
/// Classical multidimensional scaling: one row of coordinates per object.
pub fn mds(
  distances : @data.DataMatrix,
  dims~ : Int32,
) -> @data.DataMatrix raise @regress.ModelError {
  ignore((distances, dims))
  ...
}

///|
/// Canonical correlations between two sets of variables measured on the
/// same rows, largest first.
pub fn cca(
  x : @data.DataMatrix,
  y : @data.DataMatrix,
) -> FixedArray[Double] raise @regress.ModelError {
  if x.rows != y.rows {
    raise LengthMismatch(expected=x.rows, got=y.rows)
  }
  ...
}

///|
/// Correspondence analysis: row and column coordinates on the principal
/// axes, and the inertia each axis carries.
#derive(Eq, Debug)
pub struct Correspondence {
  row_coordinates : @data.DataMatrix
  col_coordinates : @data.DataMatrix
  inertia : FixedArray[Double]
}

///|
pub fn correspondence(
  table : @data.Table,
  dims~ : Int32,
) -> Correspondence raise @regress.ModelError {
  ignore((table, dims))
  ...
}

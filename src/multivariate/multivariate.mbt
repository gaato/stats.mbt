///|
/// Multivariate methods on a data matrix (observations × variables).
/// Every method here is an eigenvalue or SVD problem, so the bodies wait
/// for dense linear algebra.

///|
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
pub struct Pca {
  loadings : @data.DataMatrix
  scores : @data.DataMatrix
  variances : FixedArray[Double]
  explained : FixedArray[Double]
}

///|
pub fn pca(x : @data.DataMatrix, scaling : Scaling) -> Pca {
  ignore((x, scaling))
  ...
}

///|
#derive(Eq, Debug)
pub(all) enum Rotation {
  None
  Varimax
  Promax(power~ : Double)
}

///|
pub struct FactorFit {
  loadings : @data.DataMatrix
  uniqueness : FixedArray[Double]
  loglik : Double
  converged : Bool
}

///|
/// Maximum-likelihood factor analysis.
pub fn factor_analysis(
  x : @data.DataMatrix,
  factors : Int,
  rotation : Rotation,
) -> FactorFit {
  ignore((x, factors, rotation))
  ...
}

///|
/// Linear (shared covariance) or quadratic discriminant analysis.
#derive(Eq, Debug)
pub(all) enum Discriminant {
  Linear
  Quadratic
}

///|
pub struct Classifier {
  classes : Int
  /// Class probabilities for new rows.
  predict : (@data.DataMatrix) -> @data.DataMatrix
}

///|
pub fn discriminant(
  x : @data.DataMatrix,
  labels : ArrayView[Int],
  kind : Discriminant,
) -> Classifier {
  ignore((x, labels, kind))
  ...
}

///|
/// ROC curve and the area under it for scores against binary labels.
pub struct Roc {
  false_positive : FixedArray[Double]
  true_positive : FixedArray[Double]
  thresholds : FixedArray[Double]
  auc : Double
}

///|
pub fn roc(scores : ArrayView[Double], positive : ArrayView[Bool]) -> Roc {
  ignore((scores, positive))
  ...
}

///|
pub struct KMeans {
  centers : @data.DataMatrix
  labels : FixedArray[Int]
  within_ss : FixedArray[Double]
  iterations : Int
}

///|
/// k-means with k-means++ seeding.
pub fn kmeans(
  x : @data.DataMatrix,
  k : Int,
  rng~ : &@num.RandomSource,
) -> KMeans {
  ignore((x, k, rng))
  ...
}

///|
#derive(Eq, Debug)
pub(all) enum Linkage {
  Single
  Complete
  Average
  Ward
}

///|
/// Agglomerative clustering: the merge sequence and heights.
pub struct Dendrogram {
  merges : FixedArray[(Int, Int)]
  heights : FixedArray[Double]
}

///|
pub fn hierarchical(distances : @data.DataMatrix, linkage : Linkage) -> Dendrogram {
  ignore((distances, linkage))
  ...
}

///|
/// Classical multidimensional scaling to `dims` dimensions.
pub fn mds(distances : @data.DataMatrix, dims : Int) -> @data.DataMatrix {
  ignore((distances, dims))
  ...
}

///|
/// Canonical correlations between two sets of variables.
pub fn cca(x : @data.DataMatrix, y : @data.DataMatrix) -> FixedArray[Double] {
  ignore((x, y))
  ...
}

///|
/// Correspondence analysis of a contingency table.
pub fn correspondence(table : @data.Table, dims : Int) -> (@data.DataMatrix, @data.DataMatrix) {
  ignore((table, dims))
  ...
}

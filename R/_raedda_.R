#' raedda: Robust Adaptive Eigenvalue Decomposition Discriminant Analysis
#'
#' RAEDDA is a robust generalization of the AMDA methodology (Bouveyron, 2014) that
#' accounts for outliers and label noise detecting observations with the lowest contributions to the overall likelihood employing impartial trimming.
#' It performs Supervised Learning in Presence of Outliers, Label Noise and Unobserved Classes employing MVN mixture model with Parsimonious covariance structure.
#' Parameters estimation is carried out by considering either a transductive or an inductive approach
#'
#' @section raedda functions:
#' RAEDDA_transductive estimates parameters employing a transductive (semi-supervised) approach
#'
#' @docType package
#' @name raedda
NULL

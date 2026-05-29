#' @title ConfMatrix: Statistical Tools for Thematic-Accuracy Quality Control
#'
#' @description
#' Statistical tools for validating thematic cartography using confusion 
#' matrices and column-wise quality control.
#'
#' @details
#' The package is organized into two main classes. Use the following links 
#' to access the reference documentation for each method:
#' 
#' \bold{\code{\link{ConfMatrix}} Class (Global and Class-level Evaluation):}
#' \itemize{
#'   \item \bold{Global Indices:} 
#'     \code{\link[=ConfMatrix]{Overall Accuracy (OA)}}, 
#'     \code{\link[=ConfMatrix]{Kappa}}, 
#'     \code{\link[=ConfMatrix]{Tau}}, 
#'     \code{\link[=ConfMatrix]{Ji.test}}, 
#'     \code{\link[=ConfMatrix]{Relative_Agreement_Index}}.
#'   \item \bold{Class-specific (\eqn{i}):} 
#'     \code{\link[=ConfMatrix]{UserAcc_i}}, 
#'     \code{\link[=ConfMatrix]{ProdAcc_i}}, 
#'     \code{\link[=ConfMatrix]{Error_i}}, 
#'     \code{\link[=ConfMatrix]{f_i}}, 
#'     \code{\link[=ConfMatrix]{g_i}}, 
#'     \code{\link[=ConfMatrix]{h_i}}.
#'   \item \bold{Averages and Combined:} 
#'     \code{\link[=ConfMatrix]{AvUserAcc}}, 
#'     \code{\link[=ConfMatrix]{AvProdAcc}}, 
#'     \code{\link[=ConfMatrix]{AvError}}, 
#'     \code{\link[=ConfMatrix]{CombUserAcc_i}}, 
#'     \code{\link[=ConfMatrix]{CombProdAcc_i}}.
#'   \item \bold{Advanced Tests:} 
#'     \code{\link[=ConfMatrix]{Marghit}} (Marginal Homogeneity), 
#'     \code{\link[=ConfMatrix]{Quasi_Independence}}.
#' }
#' 
#' \bold{\code{\link{QCCS}} Class (Quality Control):}
#' \itemize{
#'   \item \bold{Tests:} 
#'     \code{\link[=QCCS]{Exact.test}}, 
#'     \code{\link[=QCCS]{Ji.test}}, 
#'     \code{\link[=QCCS]{JiGlobal.test}}.
#' }
#'
#' @name ConfMatrix-package
#' @aliases ConfMatrix-package
#' @docType package
"_PACKAGE"
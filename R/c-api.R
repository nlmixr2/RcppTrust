#' Thread-safe C interfaces
#'
#' Every optimizer in this package is implemented as a C-callable core with
#' no R API calls and no global state, so each can run concurrently from
#' OpenMP/parallel C++ code (e.g. as an inner or outer optimizer in
#' \code{nlmixr2est}). The R functions [trust()], [steihaug()], [bobyqa()]
#' and [newuoa()] are thin wrappers around these same cores.
#'
#' \tabular{lll}{
#'   \strong{Solver} \tab \strong{Header} \tab \strong{Entry point / free} \cr
#'   trust \tab \code{trust_types.h} \tab \code{trust_solve_c} / \code{trust_result_free} \cr
#'   BOBYQA \tab \code{minqa_types.h} \tab \code{bobyqa_solve_c} / \code{minqa_result_free} \cr
#'   NEWUOA \tab \code{minqa_types.h} \tab \code{newuoa_solve_c} / \code{minqa_result_free} \cr
#'   Steihaug \tab \code{steihaug_types.h} \tab \code{steihaug_solve_c} / \code{steihaug_result_free} \cr
#' }
#'
#' A consuming package adds \code{RcppTrust} to \code{LinkingTo} and
#' \code{Imports}, includes \code{RcppTrust.h}, and resolves the pointers at
#' load time from [.RcppTrustPtr()] (slots, in order: \code{trust_solve_c},
#' \code{trust_result_free}, \code{bobyqa_solve_c}, \code{newuoa_solve_c},
#' \code{minqa_result_free}, \code{steihaug_solve_c},
#' \code{steihaug_result_free}). See \code{vignette("RcppTrust")} for a
#' worked example.
#'
#' @name minqa_c_api
#' @aliases steihaug_c_api trust_c_api
#' @keywords internal
NULL

#' Textile data on breaking strength by starch type and chip thickness
#'
#' Furry (1939) analyzes the breaking strength of a starch chip
#' as a function of the chip’s thickness (measured in 10^-4 inches) and the type
#' of plant from which the starch was derived (corn, canna, or potato).
#' @docType data
#'
#' @usage data(textile)
#'
#' @format A data frame with 49 rows and 3 variables:
#' \describe{
#'   \item{strength}{the response, the breaking strength.}
#'   \item{film}{the chip's film thickness, measured in 10^-4 inches.}
#'   \item{starch}{the chip's starch component: canna, corn, or potato}
#' }
#'
#' @references
#' \insertRef{Furry}{slgf}
#'
#' @keywords datasets
#'
#'
"textile"

#' Locknut data on torque required to tighten a fixture by plating method
#'
#' Meek and Ozgur (1991) analyzes the torque required to strengthen a
#' fixture (bolt or mandrel) as a function of the fixture's plating method
#' (cadmium and wax, heat treating, and phosphate and oil, denoted CW, HT,
#' and PO, respectively).
#'
#' @docType data
#'
#' @usage data(torque)
#'
#' @format A data frame with 60 rows and 3 variables:
#' \describe{
#'   \item{Torque}{the response, the torque required to tighten the fixture.}
#'   \item{Fixture}{the type of fixture, bolt or mandrel.}
#'   \item{Plating}{the plating treatment, CW, HT, or PO.}
#' }
#'
#' @references
#' \insertRef{torque}{slgf}
#'
#' @keywords datasets
#'
#'
"torque"

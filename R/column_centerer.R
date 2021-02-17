#' Column centerer for a design matrix.
#'
#' @description \code{column_centerer} Centers the columns of a matrix by the column mean.
#'
#' @param mm a model matrix.
#'
#' @return \code{column_centerer} centers the columns of a design matrix by each column's mean.
#'
#' @export
#' @examples
#' set.seed(314159)
#' test.data <- data.frame("y"=c(rnorm(10,0,1), rnorm(10,3,1), rnorm(10,5,3)),
#'                         "x1"=c(rep("A",10), rep("B",10), rep("C",10)),
#'                         "x2"=rnorm(30,0,1))
#' m <- lm(y~x1+x2, data=test.data)
#' mm <- model.matrix(m)
#' column_centerer(mm)
#'
column_centerer <- function(mm){
  mm <- mm-matrix(rep(as.vector(apply(mm,2,mean)),nrow(mm)),nrow=nrow(mm),byrow=TRUE)
  return(mm)
}

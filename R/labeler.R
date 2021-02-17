#' Labels observations according to group membership.
#'
#' @description \code{labeler} Returns a Boolean indicator for each observation's group membership for a two-way layout.
#'
#' @param nrows the number of rows in the data matrix.
#'
#' @param ncols the number of columns in the data matrix.
#'
#' @param combo.iteration the index of the grouping scheme under consideration.
#'
#' @return \code{labeler} returns a vector of 1s and 0s corresponding to the input vector's group membership by index.
#'
#' @export
#'
#'
labeler=function(nrows,ncols,combo.iteration){
  Gid=rep(0,nrows)
  Gid[combo.iteration]=1
  labels=rep(Gid,each=ncols)
  return(labels)
}

#' Groupings finder for two-way layouts.
#'
#' @importFrom utils combn
#'
#' @description \code{groupings} Computes the possible grouping schemes for a two-way layout with r rows and c columns.
#'
#' @param data_matrix an r by c data matrix.
#'
#' @return \code{groupings} returns the unique possible row-wise groupings of the input two-way layout.
#'
#' @export
#' @examples
#' # Determine the possible row-wise groupings for an 8 by 5 matrix.

#' groupings(matrix(NA, nrow=8, ncol=4))
#'
groupings=function(data_matrix){
r=nrow(data_matrix)
c=ncol(data_matrix)
if(r>20){
  cat("Due to computation time, \n groupings is unavailable for a > 20 rows\n")
  stop
}
else{
  possible=as.list(rep(NA,trunc(r/2)))
  if(r>3){
    if(r%%2==1){
      niter=round((r-1)/2)
      for(j in 2:niter){
        combs=utils::combn(r,j)
        possible[[j]]=combs
      }
    }
    else {
      niter=round(r/2-1)
      if(niter>1){
        for(j in 2:niter){
          combs=utils::combn(r,j)
          possible[[j]]=combs
        }
      }
      niter=niter+1
      combs=utils::combn(r,niter)
      if(ncol(combs)%%2==0){
        combs=combs[,-seq(.5*ncol(combs)+1,ncol(combs))]
      }else{
        combs=combs[,-seq(.5*ncol(combs)+2,ncol(combs))]
      }
      possible[[trunc(r/2)]]=combs
    }
  }
}
return(possible[-1])
}


#' Converts a two-way layout into tall format with row and column index labels.
#'
#' @description \code{maketall} Converts a two-way layout into tall format with row and column index labels.
#'
#' @param data_matrix an r by c data matrix.
#'
#' @return \code{maketall} returns a data frame containing the original observations, row labels, and column labels.
#'
#' @export
#' @examples
#' library(slgf)
#' data(lymphoma)
#' maketall(lymphoma)
#'

maketall=function(data_matrix){
  n=prod(dim(data_matrix)) #compute the number of elements in the matrix
  y=rep(NA,n) #create empty vector for observations
  rows=rep(NA,n) #create empty vector for row indicator
  cols=rep(NA,n) #create empty vector for column indicator
  r=nrow(data_matrix) #count number of rows
  c=ncol(data_matrix) #count number of columns
  counter=0
  for(i in 1:r)
    for(j in 1:c){
      counter=counter+1
      rows[counter]=i #set row label
      cols[counter]=j #set column label
      y[counter]=data_matrix[i,j] #fill vector with observations
    }
  data_and_labels=matrix=cbind(y,rows,cols)
  return(data_and_labels)
}

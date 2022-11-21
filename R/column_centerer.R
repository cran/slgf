.column_centerer <- function(mm){
  mm <- mm-matrix(rep(as.vector(apply(mm,2,mean)),nrow(mm)),nrow=nrow(mm),byrow=TRUE)
  return(mm)
}

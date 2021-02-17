

#' @name labeler

#' @rdname labeler
#' @export
labeler=function(nrows,ncols,combo.iteration){
  Gid=rep(0,nrows)
  Gid[combo.iteration]=1
  labels=rep(Gid,each=ncols)
  return(labels)
}

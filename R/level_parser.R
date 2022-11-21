.level_parser <- function(scheme_string, lgf_vec){
  if(scheme_string=="None"){
    return(NA)
  }else{
    split_scheme <- lapply(lapply(strsplit(scheme_string,"}"), function(x) sub(".","",x)), function(x) strsplit(x,","))[[1]]
    group1 <- levels(lgf_vec)%in%split_scheme[[1]]
    group2 <- levels(lgf_vec)%in%split_scheme[[2]]
    return(list(group1=group1, group2=group2))
  }
}

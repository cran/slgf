.grouping_info <- function(dataf, lgf_beta, lgf_Sigma, min_levels_beta, min_levels_Sigma, response, same_scheme){
  lgf_beta_levels <- levels(dataf[,which(colnames(dataf)==lgf_beta)])
  lgf_Sigma_levels <- levels(dataf[,which(colnames(dataf)==lgf_Sigma)])

  K_beta <- length(lgf_beta_levels)
  K_Sigma <- length(lgf_Sigma_levels)

  # create grouping schemes
  if(!is.na(lgf_beta)){
    schemes_beta <- matrix(NA, ncol=1, nrow=K_beta)
    for(m in min_levels_beta:trunc(K_beta/2)){
      schemes_beta <- cbind(schemes_beta, apply(combn(1:K_beta, m=m), 2, function(x) as.numeric(as.factor(lgf_beta_levels))%in%x))
    }
    rownames(schemes_beta) <- lgf_beta_levels
    schemes_beta <- schemes_beta[,-1]
    nschemes_beta <- ncol(schemes_beta)
  }

  if(!is.na(lgf_Sigma)){
    schemes_Sigma <- matrix(NA, ncol=1, nrow=K_Sigma)
    for(m in min_levels_Sigma:trunc(K_Sigma/2)){
      schemes_Sigma <- cbind(schemes_Sigma, apply(combn(1:K_Sigma, m=m), 2, function(x) as.numeric(as.factor(lgf_Sigma_levels))%in%x))
    }
    rownames(schemes_Sigma) <- lgf_Sigma_levels
    schemes_Sigma <- schemes_Sigma[,-1]
    nschemes_Sigma <- ncol(schemes_Sigma)
  }

  if(is.na(lgf_beta)){nschemes_beta <- 0;K_beta <- 0}
  if(is.na(lgf_Sigma)){nschemes_Sigma <- 0;K_Sigma <- 0}
  if(K_Sigma>0 & K_Sigma%%2==0){
    duplicated <- which(apply(schemes_Sigma, 2, function(x) sum(x==TRUE))==apply(schemes_Sigma, 2, function(x) sum(x==FALSE)))
    schemes_Sigma <- schemes_Sigma[,-duplicated[(length(duplicated)/2+1):(length(duplicated))]]
    nschemes_Sigma <- ncol(schemes_Sigma)
  }
  if(K_beta>0 & K_beta%%2==0){
    duplicated <- which(apply(schemes_beta, 2, function(x) sum(x==TRUE))==apply(schemes_beta, 2, function(x) sum(x==FALSE)))
    schemes_beta <- schemes_beta[,-duplicated[(length(duplicated)/2+1):(length(duplicated))]]
    nschemes_beta <- ncol(schemes_beta)
  }
  # if((nschemes_beta*nschemes_Sigma>2000 & same_scheme==FALSE) | nschemes_Sigma>2000 | nschemes_beta>2000){
  #   stop(paste0("Error 4\nThere are too many grouping schemes (", nschemes_beta*nschemes_Sigma, ") to consider.\nSpecify a different latent grouping factor, or, make min_levels\ncloser to ",trunc(max(K_beta,K_Sigma)/2),", or, set same_scheme=TRUE."))
  # }
  if(K_Sigma>0 & K_Sigma==2){
    schemes_Sigma <- matrix(schemes_Sigma, ncol=1)
  }

  if(is.na(lgf_beta)){
    return(list(schemes_Sigma=as.matrix(schemes_Sigma), K_Sigma=K_Sigma))
  }
  if(is.na(lgf_Sigma)){
    return(list(schemes_beta=as.matrix(schemes_beta), K_beta=K_beta))
  }
  if(same_scheme==FALSE){
    return(list(schemes_beta=as.matrix(schemes_beta), schemes_Sigma=as.matrix(schemes_Sigma), K_beta=K_beta, K_Sigma=K_Sigma))
  }
  if(same_scheme==TRUE){
    return(list(schemes=as.matrix(schemes_beta), K=K_beta))
  }
}

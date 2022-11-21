.prior_maker <- function(models_out){
  beta_indicator <- models_out$models$Scheme.beta!="None"
  beta_indicator[as.logical(beta_indicator)] <- "; group-based_effects=TRUE"
  beta_indicator[!as.logical(beta_indicator)] <- "; group-based_effects=FALSE"

  Sigma_indicator <- models_out$models$Scheme.Sigma!="None"
  Sigma_indicator[as.logical(Sigma_indicator)] <- "; group-based_variances=TRUE"
  Sigma_indicator[!as.logical(Sigma_indicator)] <- "; group-based_variances=FALSE"

  classes <- as.factor(paste0(models_out$models$Model, beta_indicator, Sigma_indicator))
  nclasses <- length(levels(classes))
  modpriors <- rep(NA,nrow(models_out$models))

  for(c in levels(classes)){
    model_index <- which(classes==c)
    modpriors[model_index] <- (1/nclasses)/(length(model_index))
  }
  return(list(modpriors=modpriors, classes=classes))
}

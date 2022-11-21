.model_builder <- function(usermodels, het, grouping_info, same_scheme, dataf,
                           lgf_beta=lgf_beta, lgf_beta_vec, lgf_Sigma_vec){
  models <- matrix(NA, ncol=4, nrow=2^16)
  model_fits <- as.list(rep(NA,2^16))
  coefficients <- as.list(rep(NA,2^16))
  variances <- as.list(rep(NA,2^16))
  gs <- as.list(rep(NA,2^16))
  complexities <- as.list(rep(NA,2^16))
  colnames(models) <- c("Model", "Scheme.beta", "Scheme.Sigma", "Log-Marginal")
  if(same_scheme==TRUE){
    schemes_beta <- grouping_info$schemes
    schemes_Sigma <- grouping_info$schemes
  }else{
    schemes_beta <- grouping_info$schemes_beta
    schemes_Sigma <- grouping_info$schemes_Sigma
  }
  mm_index <- 1
  for(m in 1:length(usermodels)){
    tempmod <- unlist(usermodels)[[m]]
    if(!grepl("group",tempmod)){
      out <- lm(as.formula(tempmod), data=dataf)
      if(out$df.residual==0){
        stop("Error 11\nAt least one model is overparametrized and has 0 residual degrees of freedom. Respecify your models.")
      }
      x1 <- names(out$model)
      # if("group"%in%x1 & lgf_beta%in%x1){
      #   errmess <- paste0("Error 16\nYou specified ", lgf_beta, " as lgf_beta, but some models contain both ", lgf_beta, " and a group effect.\nSpecify a different lgf_beta, or, remove models containing both ", lgf_beta, " and a group effect.")
      #   stop(errmess)
      # }
      model_fits[[mm_index]] <- out
      coefficients[[mm_index]] <- out$coefficients
      variances[[mm_index]] <- summary(out)$sigma^2
      complexities[[mm_index]] <- length(out$coefficients) + 1
      models[mm_index,] <- c(tempmod, "None", "None", NA)
      mm_index <- mm_index + 1
    }
    if(!grepl("group",tempmod) & het[m]==1){
      for(i in 1:ncol(schemes_Sigma)){
        group1 <- paste0("{",paste0(levels(lgf_Sigma_vec)[schemes_Sigma[,i]], collapse=","),"}")
        group2 <- paste0("{",paste0(levels(lgf_Sigma_vec)[!schemes_Sigma[,i]], collapse=","),"}")
        tempscheme <- paste0(group1,group2)
        models[mm_index,] <- c(tempmod, "None", tempscheme, NA)
        model_fits[[mm_index]] <- out
        coefficients[[mm_index]] <- out$coefficients
        complexities[[mm_index]] <- length(out$coefficients) + 2
        variances[[mm_index]] <- NA
        mm_index <- mm_index + 1
      }
    }
    if(grepl("group",tempmod)){
      for(i in 1:ncol(schemes_beta)){
        group1 <- paste0("{",paste0(levels(lgf_beta_vec)[schemes_beta[,i]], collapse=","),"}")
        group2 <- paste0("{",paste0(levels(lgf_beta_vec)[!schemes_beta[,i]], collapse=","),"}")
        tempscheme <- paste0(group1,group2)
        tempdf <- data.frame(dataf, "group"=lgf_beta_vec%in%levels(lgf_beta_vec)[schemes_beta[,i]])
        tempdf$group[tempdf$group==TRUE] <- group1
        tempdf$group[tempdf$group==FALSE] <- group2
        out <- lm(as.formula(tempmod), data=tempdf)
        if(out$df.residual==0){
          stop("Error 11\nAt least one model is overparametrized and has 0 residual degrees of freedom. Respecify your models.")
        }
        x1 <- names(out$model)
        # if("group"%in%x1 & lgf_beta%in%x1){
        #   errmess <- paste0("Error 16\nYou specified ", lgf_beta, " as lgf_beta, but some models contain both ", lgf_beta, " and a group effect.\nSpecify a different lgf_beta, or, remove models containing both ", lgf_beta, " and a group effect.")
        #   stop(errmess)
        # }
        models[mm_index,] <- c(tempmod, tempscheme, "None", NA)
        model_fits[[mm_index]] <- out
        coefficients[[mm_index]] <- out$coefficients
        complexities[[mm_index]] <- length(out$coefficients) + 1
        variances[[mm_index]] <- summary(out)$sigma^2
        mm_index <- mm_index + 1
      }
    }
    if(same_scheme==TRUE){
      if(grepl("group",tempmod) & het[m]==1){
        for(i in 1:ncol(schemes_Sigma)){
          group1 <- paste0("{",paste0(levels(lgf_beta_vec)[schemes_Sigma[,i]], collapse=","),"}")
          group2 <- paste0("{",paste0(levels(lgf_beta_vec)[!schemes_Sigma[,i]], collapse=","),"}")
          tempscheme <- paste0(group1,group2)
          tempdf <- data.frame(dataf, "group"=as.character(lgf_beta_vec%in%levels(lgf_beta_vec)[schemes_beta[,i]]))
          tempdf$group[tempdf$group==TRUE] <- group1
          tempdf$group[tempdf$group==FALSE] <- group2
          out <- lm(as.formula(tempmod), data=tempdf)
          if(out$df.residual==0){
            stop("Error 11\nAt least one model is overparametrized and has 0 residual degrees of freedom. Respecify your models.")
          }
          x1 <- names(out$model)
          # if("group"%in%x1 & lgf_beta%in%x1){
          #   errmess <- paste0("Error 16\nYou specified ", lgf_beta, " as lgf_beta, but some models contain both ", lgf_beta, " and a group effect.\nSpecify a different lgf_beta, or, remove models containing both ", lgf_beta, " and a group effect.")
          #   stop(errmess)
          # }
          models[mm_index,] <- c(tempmod, tempscheme, tempscheme, NA)
          model_fits[[mm_index]] <- out
          coefficients[[mm_index]] <- out$coefficients
          complexities[[mm_index]] <- length(out$coefficients) + 2
          variances[[mm_index]] <- NA
          mm_index <- mm_index + 1
        }
      }
    }
    if(same_scheme==FALSE){
      if(grepl("group",tempmod) & het[m]==1){
        for(i in 1:ncol(schemes_Sigma)){
          for(j in 1:ncol(schemes_beta)){
            group1 <- paste0("{",paste0(levels(lgf_beta_vec)[schemes_beta[,j]], collapse=","),"}")
            group2 <- paste0("{",paste0(levels(lgf_beta_vec)[!schemes_beta[,j]], collapse=","),"}")
            tempdf <- data.frame(dataf, "group"=lgf_beta_vec%in%levels(lgf_beta_vec)[schemes_beta[,j]])
            tempdf$group[tempdf$group==TRUE] <- group1
            tempdf$group[tempdf$group==FALSE] <- group2
            tempschemebeta <- paste0(group1,group2)
            group1 <- paste0("{",paste0(levels(lgf_Sigma_vec)[schemes_Sigma[,i]], collapse=","),"}")
            group2 <- paste0("{",paste0(levels(lgf_Sigma_vec)[!schemes_Sigma[,i]], collapse=","),"}")
            tempschemeSigma <- paste0(group1,group2)
            out <- lm(as.formula(tempmod), data=tempdf)
            if(out$df.residual==0){
              stop("Error 11\nAt least one model is overparametrized and has 0 residual degrees of freedom. Respecify your models.")
            }
            x1 <- names(out$model)
            # if("group"%in%x1 & lgf_beta%in%x1){
            #   errmess <- paste0("Error 16\nYou specified ", lgf_beta, " as lgf_beta, but some models contain both ", lgf_beta, " and a group effect.\nSpecify a different lgf_beta, or, remove models containing both ", lgf_beta, " and a group effect.")
            #   stop(errmess)
            # }
            models[mm_index,] <- c(tempmod, tempschemebeta, tempschemeSigma, NA)
            model_fits[[mm_index]] <- out
            coefficients[[mm_index]] <- out$coefficients
            complexities[[mm_index]] <- length(out$coefficients) + 2
            variances[[mm_index]] <- NA
            mm_index <- mm_index + 1
          }
        }
      }
    }

    nmodels <- which(is.na(models[,1]))[1]-1
    model_fits <- model_fits[1:nmodels]
    coefficients <- coefficients[1:nmodels]
    variances <- variances[1:nmodels]
    complexities <- complexities[1:nmodels]
    gs <- gs[1:nmodels]
  }
  models <- as.data.frame(models[1:nmodels,])
  if(nmodels==1){
    stop("Error 12\nThere is only one model under consideration. Model selection requires at least two models to be considered.")
  }
  rownames(models) <- 1:nmodels
  return(list(models=models,
              model_fits=model_fits,
              coefficients=coefficients,
              variances=variances, gs=gs,
              complexities=complexities))
}




#' Bayesian Model Selection with Latent Group-Based Regression Effects and Heteroscedasticity
#'
#' @importFrom Rdpack reprompt
#' @importFrom stats lm
#' @importFrom stats as.formula
#' @importFrom stats model.frame
#' @importFrom stats model.matrix
#' @importFrom stats optim
#' @importFrom stats uniroot
#' @importFrom stats var
#' @importFrom utils combn
#'
#' @description \code{ms_slgf} Implements the model selection method proposed by \insertCite{metzger2019}{slgf}.
#'
#' @references
#' \insertRef{metzger2019}{slgf}
#'
#' @author Thomas A. Metzger and Christopher T. Franck
#'
#' @param dataf A data frame containing a continuous response, at least one categorical predictor, and any other covariates of interest. This data frame should not contain column names with the character string \code{group}.
#' @param response A character string indicating the column of \code{dataf} that contains the response.
#' @param lgf_beta An optional character string indicating the column of `dataf` that contains the suspected latent grouping factor (SLGF) for the regression effects.
#' @param min_levels_beta A numeric value indicating the minimum number of levels of `lgf_beta` that can comprise a group. Defaults to 1.
#' @param lgf_Sigma An optional character string indicating the column of `dataf` that contains the suspected latent grouping factor (SLGF) for the residual variances.
#' @param min_levels_Sigma A numeric value indicating the minimum number of levels of `lgf_Sigma` that can comprise a group. Defaults to 1.
#' @param same_scheme A Boolean operator indicating whether the schemes for `lgf_beta` and `lgf_Sigma` must be the same.
#' @param usermodels A list of length \code{M} where each element contains a string of *R* class \code{formula} or \code{character} indicating the models to consider. The term \code{group} should be used to replace the name of the SLGF in models with group-based regression effects.
#' @param het A vector of 0s and 1s of length \code{M}. If the mth element of \code{het} is 0, then the mth model of \code{usermodels} is considered in a homoscedastic context only; if the mth element of \code{het} is 1, the mth model of \code{usermodels} is considered in both homoscedastic and heteroscedastic contexts.
#' @param prior A character string \code{"flat"} or \code{"zs"} indicating whether to implement the flat or Zellner-Siow mixture g-prior on regression effects, respectively. Defaults to \code{"flat"}.
#' @param m0 An integer value indicating the minimum training sample size. Defaults to NULL. If no value is provided, the lowest value that leads to convergence for all considered posterior model probabilities will be used. If the value provided is too low for convergence, it will be increased automatically.
#'
#' @import numDeriv
#' @return \code{ms_slgf} returns a list of five elements if the flat prior is used, and six elements if the Zellner-Siow mixture g-prior is used:\cr
#'     1) \code{models}, an \code{M} by 7 matrix where columns contain the model selection results and information for each model, including: \cr
#'     - \code{Model}, the formula associated with each model; \cr
#'     - \code{Scheme.beta}, the grouping scheme associated with the fixed effects; \cr
#'     - \code{Scheme.Sigma}, the grouping scheme associated with the variances; \cr
#'     - \code{Log-Marginal}, the fractional log-marginal likelihood associated with each model; \cr
#'     - \code{FmodProb}, the fractional posterior probability associated with each model; \cr
#'     - \code{ModPrior}, the prior assigned to each model; \cr
#'     - \code{Cumulative}, the cumulative fractional posterior probability associated with a given model and the previous models; \cr
#'     2) \code{class_probabilities}, a vector containing cumulative posterior probabilities associated with each model class; \cr
#'     3) \code{coefficients}, MLEs for each model's regression effects; \cr
#'     4) \code{variances}, MLEs based on concentrated likelihood for each model's variance(s); \cr
#'     5) \code{gs}, MLEs based on concentrated likelihood for each model's \code{g}; only included if \code{prior="zs"}.

#' @export
#' @examples
#' # Analyze the smell and textile data sets.
#'
#' library(numDeriv)
#'
#'
# smell data set
#' data(smell)
#' out_smell <- ms_slgf(dataf = smell, response = "olf", het=c(1,1),
#'                      lgf_beta = "agecat", lgf_Sigma = "agecat",
#'                      same_scheme=TRUE, min_levels_beta=1, min_levels_Sigma=1,
#'                      usermodels = list("olf~agecat", "olf~group"), m0=4)
#' out_smell$models[1:5,]
#' out_smell$coefficients[[46]]
#' out_smell$variances[[46]]
#'
#' # textile data set
#' data(textile)
#' out_textile <- ms_slgf(dataf = textile, response = "strength",
#'                      lgf_beta = "starch", lgf_Sigma = "starch",
#'                      same_scheme=FALSE, min_levels_beta=1, min_levels_Sigma=1,
#'                      usermodels = list("strength~film+starch", "strength~film*starch",
#'                                        "strength~film+group", "strength~film*group"),
#'                      het=c(1,1,1,1), prior="flat", m0=8)
#' out_textile$models[1:5,c(1,2,3,5)]
#' out_textile$class_probabilities
#' out_textile$coefficients[31]
#' out_textile$variances[31]


ms_slgf <- function(dataf,
                    response,
                    lgf_beta,
                    min_levels_beta=1,
                    lgf_Sigma,
                    min_levels_Sigma=1,
                    same_scheme=TRUE,
                    usermodels,
                    het=rep(0,length(usermodels)),
                    prior="flat", m0=NULL){

  # initial error checking
  {
  if(length(het)!=length(usermodels)){
    stop("Error 0\nlength(het) must equal length(usermodels).")
  }
  if(same_scheme==TRUE & (is.na(lgf_beta)|is.na(lgf_Sigma))){
    stop("Error 1\nsame_scheme cannot be TRUE when lgf_beta=NA or lgf_Sigma=NA")
  }
  if(is.na(lgf_Sigma) & sum(het)>=1 & same_scheme==FALSE){
    stop("Error 13\nYou must specify a latent grouping factor for the variance.")
  }
  if(same_scheme==TRUE & (lgf_beta!=lgf_Sigma)){
    stop("Error 2\nsame_scheme cannot be TRUE when lgf_beta!=lgf_Sigma.")
  }
  if(is.na(lgf_beta)){min_levels_beta <- 0}
  if(is.na(lgf_Sigma)){min_levels_Sigma <-0}
  if(min_levels_beta!=min_levels_Sigma & same_scheme==TRUE){
    stop("Error 3\nsame_scheme cannot be TRUE when min_levels_beta!=min_levels_Sigma.")
  }
  if(!is.null(m0)){
    if(m0%%1!=0){
      stop("Error 17\nm0 must be an integer.")
    }
  }
  if(sum(het)==0 & !is.na(lgf_Sigma)){
    stop("Error 18\nYou specified lgf_Sigma, but there are no heteroscedastic models specified.\nLet lgf_Sigma=NA, or, change het() to contain at least one 1.")
  }
  }

  # ensure variable types are correct
  {
  response <- as.character(response)
  lgf_beta <- as.character(lgf_beta)
  lgf_Sigma <- as.character(lgf_Sigma)
  dataf <- as.data.frame(dataf)
  dataf[,which(colnames(dataf)==response)] <- as.numeric(as.character(dataf[,which(colnames(dataf)==response)]))
  dataf[,which(colnames(dataf)==lgf_beta)] <- as.factor(dataf[,which(colnames(dataf)==lgf_beta)])
  dataf[,which(colnames(dataf)==lgf_Sigma)] <- as.factor(dataf[,which(colnames(dataf)==lgf_Sigma)])
  y <- dataf[,which(colnames(dataf)==response)]
  N <- length(y)
  }

  grouping_info <- .grouping_info(dataf, lgf_beta, lgf_Sigma, min_levels_beta, min_levels_Sigma, response, same_scheme)

  # second round of error checking
  {
  if(sum(grepl("group",usermodels))>0 & is.null(grouping_info$K_beta)){
    if(grouping_info$K==2){
      stop("Error 10\nGroup-based regression effects are equivalent to categorical effects when the categorical predictor has only 2 levels.\nSpecify a latent grouping factor with 3 or more levels, or, remove models with group-based regression effects.")
    }
  }
  if(sum(grepl("group",usermodels))>0 & is.null(grouping_info$K)){
    if(grouping_info$K_beta==2){
      stop("Error 10\nGroup-based regression effects are equivalent to categorical effects when the categorical predictor has only 2 levels.\nSpecify a latent grouping factor with 3 or more levels, or, remove models with group-based regression effects.")
    }
  }
  if(same_scheme==TRUE){
    if(grouping_info$K/2<min_levels_beta | grouping_info$K/2<min_levels_Sigma){
      stop("Error 5\nmin_levels_beta and/or min_levels_Sigma cannot be more than half the number of levels of the corresponding LGF.\nSpecify a different latent grouping factor, or, reduce min_levels_beta and/or min_levels_Sigma.")
    }
  }
  if(is.na(lgf_beta)){
    if(grouping_info$K_Sigma/2<min_levels_Sigma){
      stop("Error 6\nmin_levels_Sigma cannot be more than half the number of levels of the corresponding LGF.\nSpecify a different latent grouping factor, or, reduce min_levels_Sigma.")
    }
  }
  if(is.na(lgf_Sigma)){
    if(grouping_info$K_beta/2<min_levels_beta){
      stop("Error 7\nmin_levels_beta cannot be more than half the number of levels of the corresponding LGF.\nSpecify a different latent grouping factor, or, reduce min_levels_beta.")
    }
  }
  if((!is.na(lgf_beta) & !is.na(lgf_Sigma)) & same_scheme==TRUE){
    if(grouping_info$K/2<min_levels_beta | grouping_info$K/2<min_levels_Sigma){
      stop("Error 8\nmin_levels_beta and/or min_levels_Sigma cannot be more than half the number of levels of the corresponding LGF.\nSpecify a different latent grouping factor, or, reduce min_levels_beta and/or min_levels_Sigma.")
    }
  }
  if((!is.na(lgf_beta) & !is.na(lgf_Sigma)) & same_scheme==FALSE){
    if(grouping_info$K_beta/2<min_levels_beta | grouping_info$K_Sigma/2<min_levels_Sigma){
      stop("Error 9\nmin_levels_beta and/or min_levels_Sigma cannot be more than half the number of levels of the corresponding LGF.\nSpecify a different latent grouping factor, or, reduce min_levels_beta and/or min_levels_Sigma.")
    }
  }
}

  lgf_Sigma_vec <- dataf[,which(colnames(dataf)==lgf_Sigma)]
  lgf_beta_vec <- dataf[,which(colnames(dataf)==lgf_beta)]

  models_out <- .model_builder(usermodels=usermodels, het=het, grouping_info=grouping_info,
                           same_scheme=same_scheme, dataf=dataf, lgf_beta=lgf_beta,
                           lgf_beta_vec=lgf_beta_vec, lgf_Sigma_vec=lgf_Sigma_vec)
  nmodels <- nrow(models_out$models)
  if(nmodels>=(2^16)){
    stop("Error 15\nToo many models are under consideration. Lower min_levels_beta or min_levels_Sigma, or,\nconsider a different set of models.")
  }
  model_order <- order(unlist(models_out$complexities), decreasing=TRUE)
  if(is.null(m0)){m0 <- max(unlist(models_out$complexities))}
  b <- m0/N

  ####################
  return_to_start <- FALSE

  if(prior=="flat"){
    fit <- TRUE
    changed <- FALSE

    while(fit==TRUE){
      return_to_start <- FALSE

      for(i in model_order){

      if(abs(b-1)<1e-10){
        stop("Error 14\nThe minimal training fraction has reached 1. There is not enough data to consider the specified models.")
      }
      if(return_to_start){break}

      tempfit <- models_out$model_fits[i][[1]]
      tempschemeSigma <- models_out$models$Scheme.Sigma[i]
      if(tempschemeSigma=="None"){
        models_out$models$`Log-Marginal`[i] <- .marginal_calculator_hom(tempfit=tempfit, N=N, b=b)
        if(is.infinite(as.numeric(models_out$models$`Log-Marginal`[i])) |
           models_out$models$`Log-Marginal`[i]==0){b <- b+(1/N); changed <- TRUE; break}
      }
      if(tempschemeSigma!="None"){

        scheme_string <- tempschemeSigma
        split_scheme <- lapply(lapply(strsplit(scheme_string,"}"), function(x) sub(".","",x)), function(x) strsplit(x,","))[[1]]
        group1 <- levels(lgf_Sigma_vec)%in%split_scheme[[1]]
        group2 <- levels(lgf_Sigma_vec)%in%split_scheme[[2]]
        parsed <- list(group1=group1, group2=group2)

        MM <- model.matrix(tempfit)
        tempP <- ncol(MM)
        n1 <- length(tempfit$residuals[lgf_Sigma_vec%in%levels(lgf_Sigma_vec)[parsed[[1]]]])
        n2 <- length(tempfit$residuals[lgf_Sigma_vec%in%levels(lgf_Sigma_vec)[parsed[[2]]]])
        lv1hat <- log(sum(tempfit$residuals[lgf_Sigma_vec%in%levels(lgf_Sigma_vec)[parsed[[1]]]]^2)/(n1-tempP))
        lv2hat <- log(sum(tempfit$residuals[lgf_Sigma_vec%in%levels(lgf_Sigma_vec)[parsed[[2]]]]^2)/(n2-tempP))

        lvloghet <- function(arg){
        phi1 <- arg[1]
        phi2 <- arg[2]
        gam1 <- exp(-phi1)
        gam2 <- exp(-phi2)
        # J <- exp(-(phi1+phi2))
        logJ <- -(phi1+phi2)
        bee <- 1
        Sigmad <- rep(gam1, N)
        Sigmad[lgf_Sigma_vec%in%levels(lgf_Sigma_vec)[parsed[[2]]]] <- gam2
        Sigma <- diag(Sigmad)
        d <- det(bee*t(MM)%*%Sigma%*%MM)^-.5
        s <- (gam1^((bee*n1/2)-1))*(gam2^((bee*n2/2)-1))
        logs <- ((bee*n1/2)-1)*(-phi1) + ((bee*n2/2)-1)*(-phi2)
        e1 <- bee*t(y)%*%Sigma%*%y
        e2 <- bee*t(y)%*%Sigma%*%MM%*%solve(t(MM)%*%Sigma%*%MM)%*%t(MM)%*%Sigma%*%y
        e <- exp(-.5*as.numeric(e1-e2))
        return(log(d)+logs+(-.5*as.numeric(e1-e2))+logJ)
      }

        mode3 <- tryCatch(optim(c(lv1hat,lv2hat), lvloghet,
                        control=list(fnscale=-1),
                        method="Nelder-Mead")$par,
                        error = function(e){return_to_start <<- TRUE})
        if(return_to_start){b <- b+(1/N); changed <- TRUE; break}

        H3 <- tryCatch(det(-1*optim(mode3, lvloghet,
                           control=list(fnscale=-1),hessian=TRUE,
                           method="Nelder-Mead")$hessian)^-.5,
                           error = function(e){return_to_start <<- TRUE})
        if(return_to_start){b <- b+(1/N); changed <- TRUE; break}

        lvloghetb <- function(arg){
        phi1 <- arg[1]; phi2 <- arg[2]
        gam1 <- exp(-phi1); gam2 <- exp(-phi2)
        J <- exp(-(phi1+phi2))
        logJ <- -(phi1+phi2)
        bee <- b
        Sigmad <- rep(gam1, N)
        Sigmad[lgf_Sigma_vec%in%levels(lgf_Sigma_vec)[parsed[[2]]]] <- gam2
        Sigma <- diag(Sigmad)
        d <- det(bee*t(MM)%*%Sigma%*%MM)^-.5
        s <- (gam1^((bee*n1/2)-1))*(gam2^((bee*n2/2)-1))
        logs <- ((bee*n1/2)-1)*(-phi1) + ((bee*n2/2)-1)*(-phi2)
        e1 <- bee*t(y)%*%Sigma%*%y
        e2 <- bee*t(y)%*%Sigma%*%MM%*%solve(t(MM)%*%Sigma%*%MM)%*%t(MM)%*%Sigma%*%y
        e <- exp(-.5*as.numeric(e1-e2))
        return((-N*bee/2)*log(2*pi)+log(d)+logs+(-.5*as.numeric(e1-e2))+logJ)
      }

        mode3b <- tryCatch(optim(c(lv1hat,lv2hat), lvloghetb,
                        control=list(fnscale=-1),
                        method="Nelder-Mead")$par,
                        error = function(e){return_to_start <<- TRUE})
        if(return_to_start){b <- b+(1/N); changed <- TRUE; break}

        H3b <- tryCatch(det(-1*optim(mode3b, lvloghetb,
                          control=list(fnscale=-1),hessian=TRUE,
                          method="Nelder-Mead")$hessian)^-.5,
                          error = function(e){return_to_start <<- TRUE})
        if(return_to_start){b <- b+(1/N); changed <- TRUE; break}

        templogPY <- (-(N-tempP)/2)*log(2*pi) - (-(N*b-tempP)/2)*log(((2*pi))) +
          (log(H3) + lvloghet(mode3)) - log(H3b) - lvloghetb(mode3)
        models_out$models$`Log-Marginal`[i] <- templogPY
        vars <- c(mode3[1], mode3[2])
        names(vars) <- c(paste0(strsplit(tempschemeSigma, "}")[[1]][1],"}"),
                         paste0(strsplit(tempschemeSigma, "}")[[1]][2],"}"))
        models_out$variances[[i]] <- exp(vars)
      }
      }
      if(i==model_order[nmodels]){fit=FALSE}
    }
  }
  if(prior=="zs"){
    fit <- TRUE
    changed <- FALSE

    SST <- sum((y-mean(y))^2)
    Z <- rep(1,N)

    while(fit==TRUE){
      return_to_start <- FALSE

      for(i in model_order){

        if(return_to_start){break}

        tempfit <- models_out$model_fits[i][[1]]
        tempschemeSigma <- models_out$models$Scheme.Sigma[i]

        if(tempschemeSigma=="None"){

          tempP <- sum(!is.na(tempfit$coefficients))-1
          R2A <- summary(tempfit)$r.squared

          faddmode <- function(gee){
          Q <- 1-R2A
          value <- -Q*(tempP+3)*(gee^3)+
            (N-tempP-4-2*(1-R2A))*(gee^2)+
            ((N*(2-R2A)-3)*gee)+N
          return(value)
        }
          faddbmode <- function(gee){
          Q <- 1-R2A
          value=-Q*(b^2)*(tempP+3)*(gee^3)+
            (b*(N*b-tempP-4)-2*Q)*(gee^2)+
            (N*b*(2-R2A)-3)*gee+N
          return(value)
        }

          addmode <- tryCatch(uniroot(faddmode,c(0,1e9),check.conv=TRUE,tol=1e-10)$root,
                              error = function(e){return_to_start <<- TRUE})
          if(return_to_start){b <- b+(1/N); changed <- TRUE; break}
          addbmode <- tryCatch(uniroot(faddbmode,c(1e-9,1e9),check.conv=TRUE,tol=1e-10)$root,
                               error = function(e){return_to_start <<- TRUE})
          if(return_to_start){b <- b+(1/N); changed <- TRUE; break}

          logfadd <- function(gee){
          value <- ((N-tempP-1)/2)*log(1+gee) +
            (-(N-1)/2)*log(1+(1-R2A)*gee) +
            -1.5*log(gee) +
            (-N/(2*gee))
          return(value)
        }
          logfaddb <- function(gee){  #Put on log-scale
          value <- ((N*b-1-tempP)/2)*log(1+b*gee) +
            (-(N*b-1)/2)*log(1+b*gee*(1-R2A)) +
            -1.5*log(gee) + (-N/(2*gee))
          return(value)
        }

          faddH <- function(gee){
          value <- .5*((((N-1)*(1-R2A)^2)/((1+gee*(1-R2A))^2))-
                         ((N-tempP-1)/((1+gee)^2))+
                         (3/(gee^2))-
                         ((2*N)/(gee^3)))
          return(value)
        }
          faddHb <- function(gee){
          value <- .5*((((N*b-1)*(b^2)*(1-R2A)^2)/((1+gee*b*(1-R2A))^2))-
                         ((N*b-tempP-1)*(b^2)/((1+b*gee)^2))+
                         (3/(gee^2))-
                         ((2*N)/(gee^3)))
          return(value)
        }

          addH <- (-1*faddH(addmode))^-.5
          addHb <- (-1*faddHb(addbmode))^-.5

          qaddLA <- lgamma((N-1)/2) + .5*log(N/2) + (-(N-1)/2)*log(SST) +
            .5*log(N) - ((N-1)/2)*log(pi) - lgamma(.5) + log(sqrt(2*pi)) +
            log(addH) + logfadd(addmode)

          qaddbLA <- lgamma((N*b-1)/2) + .5*log(N/2) + (-(N*b-1)/2)*log(SST) +
            .5*log(N) - ((N*b-1)/2)*log(pi) - lgamma(.5) + log(sqrt(2*pi)) +
            log(addHb) + logfaddb(addbmode) + (-N*b/2)*log(b)

          templogPY <- qaddLA - qaddbLA

          models_out$models$`Log-Marginal`[i] <- templogPY
          models_out$gs[[i]] <- addmode
      }
        if(tempschemeSigma!="None"){

        scheme_string <- tempschemeSigma
        split_scheme <- lapply(lapply(strsplit(scheme_string,"}"), function(x) sub(".","",x)), function(x) strsplit(x,","))[[1]]
        group1 <- levels(lgf_Sigma_vec)%in%split_scheme[[1]]
        group2 <- levels(lgf_Sigma_vec)%in%split_scheme[[2]]
        parsed <- list(group1=group1, group2=group2)

        MM <- .column_centerer(model.matrix(tempfit))[,!is.na(tempfit$coefficients)]
        MM <- MM[,-1]
        tempP <- ncol(MM); if(is.null(tempP)){tempP <- 1}
        n1 <- length(tempfit$residuals[lgf_Sigma_vec%in%levels(lgf_Sigma_vec)[parsed[[1]]]])
        n2 <- length(tempfit$residuals[lgf_Sigma_vec%in%levels(lgf_Sigma_vec)[parsed[[2]]]])
        R2A <- summary(tempfit)$r.squared
        gamhat <- 1/var(y)
        lvhat <- log(1/gamhat)

        lvlogf3 <- function(arg){
          gee <- arg[1]; u <- arg[2]; v <- arg[3]
          gam1 <- exp(-u); gam2 <- exp(-v)
          J <- exp(-(u+v)); logJ <- -(u+v)
          Sigmad <- rep(gam1, N)
          Sigmad[lgf_Sigma_vec%in%levels(lgf_Sigma_vec)[parsed[[2]]]] <- gam2
          S <- diag(Sigmad)
          Z_S <- S%*%Z%*%solve(t(Z)%*%S%*%Z)%*%t(Z)%*%S

          value <- (n1/2-1)*log(gam1) + (n2/2-1)*log(gam2) - ((tempP/2)*log(gee)) +
            .5*log(det(t(MM)%*%S%*%MM)) -
            .5*log(det(t(Z)%*%S%*%Z)) -
            .5*log(det(((gee+1)/gee)*t(MM)%*%S%*%MM-t(MM)%*%Z_S%*%MM)) -
            .5*(t(y)%*%S%*%y - t(y)%*%Z_S%*%y -
                t(y)%*%(S-Z_S)%*%MM%*%solve(((gee+1)/gee)*t(MM)%*%S%*%MM - t(MM)%*%Z_S%*%MM)%*%t(MM)%*%(S-Z_S)%*%y) -
            1.5*log(gee) - (N/(2*gee)) + logJ
          return(value)
        }
        lvlogf3b <- function(arg){
          gee <- arg[1]; u <- arg[2]; v <- arg[3]
          gam1 <- exp(-u); gam2=exp(-v)
          J <- exp(-(u+v)); logJ <- -(u+v)
          Sigmad <- rep(gam1, N)
          Sigmad[lgf_Sigma_vec%in%levels(lgf_Sigma_vec)[parsed[[2]]]] <- gam2
          S <- diag(Sigmad)
          Z_S <- S%*%Z%*%solve(t(Z)%*%S%*%Z)%*%t(Z)%*%S
          bg <- b*gee

          value <- (n1*b/2-1)*log(gam1) + (n2*b/2-1)*log(gam2) - ((tempP/2)*log(gee)) +
            .5*log(det(t(MM)%*%S%*%MM)) -
            .5*log(det(t(Z)%*%S%*%Z)) -
            .5*log(det(((bg+1)/bg)*t(MM)%*%S%*%MM-t(MM)%*%Z_S%*%MM))-
            .5*b*(t(y)%*%S%*%y-t(y)%*%Z_S%*%y-
                  t(y)%*%(S-Z_S)%*%MM%*%solve(((bg+1)/bg)*t(MM)%*%S%*%MM - t(MM)%*%Z_S%*%MM)%*%t(MM)%*%(S-Z_S)%*%y) -
            1.5*log(gee) - (N/(2*gee)) + logJ
          return(value)
        }

        suppressWarnings({
        lvmode3 <- tryCatch(optim(c(N,lvhat,lvhat),
                         lvlogf3,
                         control=list(fnscale=-1),
                         method="Nelder-Mead")$par,
                         error = function(e){return_to_start <<- TRUE})
        })
        if(return_to_start){b <- b+(1/N); changed <- TRUE; break}
        suppressWarnings({
        lvmode3b <- tryCatch(optim(c(N,lvhat,lvhat),
                          lvlogf3b,
                          control=list(fnscale=-1),
                          method="Nelder-Mead")$par,
                          error = function(e){return_to_start <<- TRUE})
        })
        if(return_to_start){b <- b+(1/N); changed <- TRUE; break}

        # H3 <- (-1*det(numDeriv::hessian(lvlogf3,lvmode3)))^-.5
        logH3 <- tryCatch(-.5*determinant(-1*numDeriv::hessian(lvlogf3,lvmode3), logarithm=TRUE)$modulus,
                          error = function(e){return_to_start <<- TRUE})
        if(return_to_start){b <- b+(1/N); changed <- TRUE; break}
        # H3b <- (-1*det(numDeriv::hessian(lvlogf3b,lvmode3b)))^-.5
        logH3b <- tryCatch(-.5*determinant(-1*numDeriv::hessian(lvlogf3b,lvmode3b), logarithm=TRUE)$modulus,
                           error = function(e){return_to_start <<- TRUE})
        if(return_to_start){b <- b+(1/N); changed <- TRUE; break}

        q3LA <- (-(N-1)/2)*log(2*pi) + .5*log(N/2) - lgamma(.5) +
          1.5*log(2*pi) + logH3 + lvlogf3(lvmode3)
        q3bLA <- (-(N*b-1)/2)*log(2*pi) + (-(tempP+1)/2)*log(b) +
          .5*log(N/2) - lgamma(.5) + 1.5*log(2*pi) +
          logH3b + lvlogf3b(lvmode3b)

        templogPY <- q3LA - q3bLA

        models_out$models$`Log-Marginal`[i] <- templogPY
        vars <- c(lvmode3[2], lvmode3[3])
        names(vars) <- c(paste0(strsplit(tempschemeSigma, "}")[[1]][1],"}"),
                         paste0(strsplit(tempschemeSigma, "}")[[1]][2],"}"))
        models_out$variances[[i]] <- exp(vars)
        models_out$gs[[i]] <- lvmode3[1]
        }
        if(is.infinite(as.numeric(templogPY)) | templogPY==0){b <- b+(1/N); changed <- TRUE; break}
      }
      if(i==model_order[nmodels]){fit=FALSE}
    }
  }

  models_out$models$`Log-Marginal` <- as.numeric(models_out$models$`Log-Marginal`)
  modevs <- exp(models_out$models$`Log-Marginal` + abs(max(models_out$models$`Log-Marginal`)))
  models_out$models$FModProb <- rep(NA, nmodels)

  #Create uniform prior by class, create classes
  prior_out <- .prior_maker(models_out)
  models_out$models$ModPrior <- prior_out[[1]]

  # compute posterior model probabilities
  for(i in 1:length(modevs)){
    models_out$models$FModProb[i]=(modevs[i]*models_out$models$ModPrior[i])/(modevs%*%models_out$models$ModPrior)
  }

  # compute posterior class probabilities
  classes <- prior_out[[2]]
  classprobs <- rep(NA,length(unique(classes)))
  names(classprobs) <- levels(classes)
  for(c in levels(classes)){
    classprobs[c] <- sum(models_out$models$FModProb[classes==c])
  }

  models_out$models <- models_out$models[order(-models_out$models$FModProb),]
  models_out$models$Cumulative <- cumsum(models_out$models$FModProb)

  # compute posterior grouping scheme probabilities
  scheme_probs_beta <- "No lgf_beta specified."
  scheme_probs_Sigma <- "No lgf_Sigma specified."

  if(!is.na(lgf_beta)){
    scheme_probs_beta <- data.frame("Scheme.beta"=unique(levels(as.factor(models_out$models$Scheme.beta))),
                                    "Cumulative"=rep(NA, length(unique(levels(as.factor(models_out$models$Scheme.beta))))))
    for(s in scheme_probs_beta$Scheme.beta){
      scheme_probs_beta$Cumulative[which(scheme_probs_beta$Scheme.beta==s)] <- sum(models_out$models$FModProb[models_out$models$Scheme.beta==s])
    }
    scheme_probs_beta <- scheme_probs_beta[order(scheme_probs_beta$Cumulative, decreasing=TRUE),]
  }

  if(!is.na(lgf_Sigma)){
    scheme_probs_Sigma <- data.frame("Scheme.Sigma"=unique(levels(as.factor(models_out$models$Scheme.Sigma))),
                                    "Cumulative"=rep(NA, length(unique(levels(as.factor(models_out$models$Scheme.Sigma))))))
    for(s in scheme_probs_Sigma$Scheme.Sigma){
      scheme_probs_Sigma$Cumulative[which(scheme_probs_Sigma$Scheme.Sigma==s)] <- sum(models_out$models$FModProb[models_out$models$Scheme.Sigma==s])
    }
    scheme_probs_Sigma <- scheme_probs_Sigma[order(scheme_probs_Sigma$Cumulative, decreasing=TRUE),]
  }


  if(changed){
    warning(paste0("Due to convergence issues, m0 was increased to ", b*N, "."))
  }

  if(prior=="flat"){
    return(list(models=models_out$models,
                class_probabilities=classprobs,
                scheme_probabilities_beta=scheme_probs_beta,
                scheme_probabilities_Sigma=scheme_probs_Sigma,
                coefficients=models_out$coefficients,
                variances=models_out$variances,
                model_fits=models_out$model_fits,m0=b*N))}
  if(prior=="zs"){
    return(list(models=models_out$models,
                class_probabilities=classprobs,
                scheme_probabilities_beta=scheme_probs_beta,
                scheme_probabilities_Sigma=scheme_probs_Sigma,
                coefficients=models_out$coefficients,
                variances=models_out$variances,
                model_fits=models_out$model_fits,
                gs=models_out$gs, m0=b*N))}
}

# add message about considering ZS if m0/N exceeds 10% and prior=flat.



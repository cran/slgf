#' Bayesian Model Selection with Latent Group-Based Regression Effects and Heteroscedasticity
#'
#' @importFrom Rdpack reprompt
#' @importFrom stats lm
#' @importFrom stats as.formula
#' @importFrom stats model.frame
#' @importFrom stats model.matrix
#' @importFrom stats optim
#' @importFrom stats uniroot
#'
#' @description \code{ms.slgf} Implements the model selection method proposed by \insertCite{metzger2019}{slgf}.
#'
#' @references
#' \insertRef{metzger2019}{slgf}
#'
#'
#' @author Thomas A. Metzger and Christopher T. Franck
#'
#' @param dataf A data frame containing a continuous response, at least one categorical predictor, and any other covariates of interest. This data frame should not contain column names with the character string \code{group}.
#' @param response A character string indicating the column of \code{dataf} that contains the response.
#' @param lgf A character string indicating the column of `dataf` that contains the suspected latent grouping factor (SLGF).
#' @param usermodels A list of length \code{M} where each element contains a string of *R* class \code{formula} or \code{character} indicating the models to consider. The term \code{group} should be used to replace the name of the slgf in models with group-based regression effects. This list must contain at least one model with group-based regression effects.
#' @param prior A character string \code{"flat"} or \code{"zs"} indicating whether to implement the flat or Zellner-Siow mixture g-prior on regression effects, respectively. Defaults to \code{"flat"}.
#' @param het A vector of 0s and 1s of length \code{M}. If the mth element of \code{het} is 0, then the mth model of \code{usermodels} is considered in a homoscedastic context only; if the mth element of \code{het} is 1, the mth model of \code{usermodels} is considered in both homoscedastic and heteroscedastic contexts.
#' @param min.levels A numeric value indicating the minimum number of levels of the SLGF that can comprise a group. Defaults to 1.
#'
#' @import numDeriv
#' @return \code{ms.slgf} returns a list of six elements:\cr
#'     1) \code{results}, an \code{M} by 11 matrix where columns contain the model selection results and information for each model, including: \cr
#'     - \code{Model}, the formula associated with each model; \cr
#'     - \code{Scheme}, the grouping scheme associated with each model; \cr
#'     - \code{Variance}, a label of whether each model is homoscedastic or heteroscedastic; \cr
#'     - \code{logFlik}, the fractional log-likelihood associated with each model; \cr
#'     - \code{Mod.Prior}, the prior assigned to each model; \cr
#'     - \code{Fmodprob}, the fractional posterior probability associated with each model; \cr
#'     - \code{Cumulative}, the cumulative fractional posterior probability associated with a given model and the previous models; \cr
#'     - \code{dataf.Index}, an index indicating which element of \code{group.datafs} contains the corresponding group \code{dataframe}; \cr
#'     - \code{mle.index}, an index indicating which element of \code{coefficients}, \code{variances}, and \code{gs} contains the corresponding estimates; \cr
#'     - \code{Model.Index}, an index indicating where the model ranks in its posterior model probability; \cr
#'     - \code{Class}, a label of the model with its group variance specification; \cr
#'     2) \code{group.datafs}, a list containing dataframes associated with each model class containing the appropriate effects, including group effects; \cr
#'     3) \code{scheme.Probs}, a \code{data.frame} containing the total posterior probability for each grouping scheme considered; \cr
#'     4) \code{class.Probs}, a \code{data.frame} containing the total posterior probability for each model class considered; \cr
#'     5) \code{coefficients}, MLEs for each model's regression effects; \cr
#'     6) \code{variances}, MLEs based on concentrated likelihood for each model's variance(s); \cr
#'     7) \code{gs}, MLEs based on concentrated likelihood for each model's \code{g}; only included if \code{prior="zs"}.

#' @export
#' @examples
#' # Fit a a heteroscedastic ANOVA example with distinct means by level of the LGF.
#'
#' library(numDeriv)
#'
#' set.seed(314159)
#' test.data <- data.frame("y"=c(rnorm(10,0,1), rnorm(10,3,1), rnorm(10,5,3)),
#'                         "x"=c(rep("A",10), rep("B",10), rep("C",10)))
#' test.models <- list("y~1", "y~x", "y~group")
#' test.out <- ms.slgf(dataf=test.data, response="y", lgf="x",
#'                     usermodels=test.models,
#'                     prior="flat", het=c(1,1,1), min.levels=1)
#' test.out$results[1:3,c(1:4,6,7)]
ms.slgf <- function(dataf,response,lgf=NA,
                    usermodels,prior="flat",
                    het=rep(0,length(usermodels)),
                    min.levels=1){

  df <- dataf
  df <- as.data.frame(df)
  #Ensure variable types are correct
  df[,which(colnames(df)==response)] <- as.numeric(as.character(df[,which(colnames(df)==response)]))
  df[,which(colnames(df)==lgf)] <- as.factor(df[,which(colnames(df)==lgf)])

  #Create a vector for the LGF
  lgf <- as.factor(df[,which(colnames(df)==as.character(lgf))])
  #K = number of levels of LGF
  K <- length(levels(lgf))

  #Create a centered vector with the response
  y <- df[,which(colnames(df)==as.character(response))]
  ybar <- mean(y)
  # y <- y-ybar

  #Sample size
  N <- length(y)

  if(min.levels>=3){
    d_groups <- groupings(matrix(NA,nrow=K,ncol=1))
    while(length(d_groups[[1]][,1])>min.levels){
      d_groups <- (d_groups[-length(d_groups)])
    }
  }
  if(min.levels==2){
    d_groups <- groupings(matrix(NA,nrow=K,ncol=1))
    if(K==3){d_groups[[1]][,3] <- c(2,3)}
  }
  if(min.levels==1){
    d_groups <- as.list(c(NA,groupings(matrix(NA,nrow=K,ncol=1))))
    d_groups[[1]] <- as.matrix(combn(K,1))
  }
  # if(min.levels<1){print("Enter a valid value for min.levels.")}
  if(K==2){
    d_groups <- list(matrix(1,nrow=1,ncol=1))
  }

  #Determine the number of possible groupings
  if(min.levels>=3){ngroups <- sum(unlist(lapply(d_groups, ncol)))}
  if(min.levels==2){ngroups <- (2^(K-1))-K-1}

  if(min.levels==1){ngroups <- (2^(K-1))-1}

  #break if there are too many grouping schemes
  if(ngroups>1024){
    paste0("There are too many grouping schemes (", ngroups,"); try increasing min.levels to reduce the number of schemes.")
    # break
  }

  #Create a minimial training sample size of 0 and
  #fractional exponent b = m0/N that will be
  #updated by more complex models as needed
  m0 <- 2
  b <- m0/N

  #Fix for heteroscedastic models but no group-based FEs
  fix <- 0
  if(sum(het)>0 & sum(grepl("group", usermodels))==0){
    fix.model <- paste0(response,"~group")
    usermodels <- as.list(c(unlist(usermodels), fix.model))
    het.orig <- het
    het=c(het, 0)
    fix=1
  }

  #Determine which models do and do not contain a group-based fixed effect

  index.lgf <- which(grepl("group",usermodels))
  index.no.lgf <- which(!(1:length(usermodels)%in%index.lgf))

  #Determine which models should have a heteroscedastic counterpart
  if(sum(het)>0){
    index.het <- seq(1:sum(het)) + length(usermodels)
  }else{
    index.het <- NA
  }

  allmodels <- unlist(usermodels)
  allmodels <- c(allmodels,paste0(allmodels[which(het==1)]))

  #Create vector indicating how many models are in each class
  repvec <- rep(1,length(allmodels))
  repvec[c(index.lgf,index.het)] <- ngroups

  #Create vector indicating which models will be fit, including
  #duplicates for group-based models
  modelvec <- rep(allmodels,times=repvec)

  #Create indicator for heteroscedastic models
  if(sum(het)>0){
    hetvec <- rep(FALSE,length(modelvec))
    hetvec[(length(modelvec)-sum(het==1)*ngroups+1):length(modelvec)] <- TRUE
  }else{
    hetvec <- rep(FALSE, length(modelvec))
  }

  #Create vector for each model's scheme
  schemevec <- rep("None",length(modelvec))

  #Create vector to store each model's marginal probability
  marginalvec <- rep(0,length(modelvec))

  classvec <- rep(1:length(repvec),times=repvec)

  #Create result output matrix
  result.mat <- as.data.frame(cbind(modelvec,schemevec,classvec,hetvec,marginalvec))
  result.mat$marginalvec <- as.numeric(as.character(result.mat$marginalvec))
  result.mat$schemevec <- as.character(result.mat$schemevec)

  fitted.models <- as.list(rep(NA,length(unique(usermodels))))
  names(fitted.models) <- unlist(unique(usermodels))

  #How can we determine whether the model results in an
  #expansion or a contraction?
  group.dfs <- as.list(rep(NA,ngroups + length(index.no.lgf)))
  all.coefs <- as.list(rep(NA, nrow(result.mat)))
  all.vars <- as.list(rep(NA, nrow(result.mat)))

  #Create list for the models with a group effect

  #Fit the models,  compute m0 - noninformative prior case
  for(tempmod in unlist(usermodels)){

    m <- which(tempmod==names(fitted.models))

    #Determine the class of the model in question.
    #If it's a class that contains only one model,
    #there is no LGF so the model is fit accordingly.

    if(grepl("group",tempmod)==FALSE){

      templm <- lm(as.formula(tempmod),data=df)
      tempmm <- model.frame(templm)

      group.dfs[[ngroups+which(index.no.lgf==which(tempmod==usermodels))]] <- tempmm
      names(group.dfs)[[ngroups+which(index.no.lgf==which(tempmod==usermodels))]] <- tempmod

      #Update m0 if model is more complex
      if((length(templm$coefficients)+3)>m0){
        m0 <- length(templm$coefficients)+2
        b <- m0/N
      }

      #For group-based variances only
      if(het[m]==TRUE){

        index <- 1
        templist <- as.list(rep(NA,ngroups))

        for(i in 1:length(d_groups)){
          for(j in 1:ncol(d_groups[[i]])){
            {
              #Create the distinct group effect
              group <- as.numeric(as.numeric(lgf)%in%d_groups[[i]][,j])
              tempscheme1 <- paste(unique(lgf[group==1]),collapse=",")
              tempscheme0 <- paste(unique(lgf[group==0]),collapse=",")
              group[group==1] <- tempscheme1
              group[group==0] <- tempscheme0

              #Update m0 if necessary
              if((length(templm$coefficients)+4)>m0){
                m0=length(templm$coefficients)+4 #<-need to delete coefs with NA
                b=m0/N
              }

              #Name the model
              tempscheme <- paste0("{",paste(c(tempscheme1,tempscheme0),collapse="}{"),"}")

              result.mat$schemevec[which(modelvec==tempmod&hetvec==TRUE)[index]] <-
                paste0("{",paste(c(tempscheme1,tempscheme0),collapse="}{"),"}")

              names(all.coefs)[[which(modelvec==tempmod&result.mat$schemevec==tempscheme)]] <- tempscheme
              names(all.vars)[[which(modelvec==tempmod&result.mat$schemevec==tempscheme)]] <- tempscheme

              tempscheme <- paste0("scheme=",tempscheme)
              names(templist)[[index]] <- paste(tempscheme)

              index=index+1
            }
          }
        }
      }

      m0 <- m0+1
      b <- m0/N
      fitted.models[[m]] <- templm
    }

    #Group-based fixed effect models
    if(grepl("group",tempmod)==TRUE){

      index <- 1
      templist <- as.list(rep(NA,ngroups))

      for(i in 1:length(d_groups)){
        for(j in 1:ncol(d_groups[[i]])){
          {
            #Create the distinct group effect
            group <- as.numeric(as.numeric(lgf)%in%d_groups[[i]][,j])
            tempscheme1 <- paste(unique(lgf[group==1]),collapse=",")
            tempscheme0 <- paste(unique(lgf[group==0]),collapse=",")
            group[group==1] <- tempscheme1
            group[group==0] <- tempscheme0

            tempdf <- data.frame(df,"group"=as.factor(group))
            # tempdf$group <- relevel(tempdf$group, ref=tempscheme0) # added 9/15

            templm <- lm(as.formula(tempmod),data=tempdf)

            templist[[index]] <- templm
            tempdf[,which(colnames(tempdf)==response)] <- tempdf[,which(colnames(tempdf)==response)] # + ybar
            group.dfs[[index]] <- tempdf
            names(group.dfs)[index] <- paste0("{",paste(c(tempscheme1,tempscheme0),collapse="}{"),"}")

            #Update m0 if necessary
            if((length(templm$coefficients)+2)>m0){
              m0=length(templm$coefficients)+2 #<-need to delete coefs with NA
              b=m0/N
            }

            #Name the model
            tempscheme <- paste0("{",paste(c(tempscheme1,tempscheme0),collapse="}{"),"}")

            result.mat$schemevec[which(modelvec==tempmod&hetvec==FALSE)[index]] <- tempscheme

            result.mat$schemevec[which(modelvec==tempmod&hetvec==TRUE)[index]] <- tempscheme

            nameindex <- which(modelvec==tempmod&result.mat$schemevec==tempscheme)

            names(all.coefs)[nameindex] <- rep(tempscheme, length(nameindex))
            names(all.vars)[nameindex] <- rep(tempscheme, length(nameindex))

            # names(all.coefs)[[which(modelvec==tempmod&hetvec==TRUE&result.mat$schemevec==tempscheme)]] <- tempscheme
            # names(all.coefs)[[which(modelvec==tempmod&hetvec==FALSE&result.mat$schemevec==tempscheme)]] <- tempscheme
            #
            # names(all.vars)[[which(modelvec==tempmod&hetvec==TRUE&result.mat$schemevec==tempscheme)]] <- tempscheme
            # names(all.vars)[[which(modelvec==tempmod&hetvec==FALSE&result.mat$schemevec==tempscheme)]] <- tempscheme

            tempscheme <- paste0("scheme=",tempscheme)
            names(templist)[[index]] <- paste(tempscheme)

            index <- index+1
          }
        }
      }

      m0 <- m0+3
      b <- m0/N
      fitted.models[[m]] <- templist

    }

  }

  #Have smaller m0 if using ZS prior
  if(prior=="zs"){
    # m0 <- 2 + 1 + 1 #two variances + g + intercept
    # m0 <- 6
    m0 <- 4
    b <- m0/N
    all.gs <- as.list(rep(NA, nrow(result.mat)))
    names(all.gs) <- paste0(result.mat$modelvec,", ", result.mat$schemevec)
  }

  names(all.coefs) <- paste0(result.mat$modelvec,", ", result.mat$schemevec)
  names(all.vars) <- paste0(result.mat$modelvec,", ", result.mat$schemevec)

  #Compute marginal model probabilities - noninformative prior case
  #Store fixed effect and variance MLEs



  if(prior=="flat"){
    for(tempmod in unlist(usermodels)){

      m <- which(tempmod==names(fitted.models))
      tempfit <- fitted.models[[m]]

      #Compute marginal model probability for homoscedastic
      #models with global fixed effects
      if(grepl("group",tempmod)==FALSE){

        tempSSResid <- sum(tempfit$residuals^2)
        tempP <- length(tempfit$coefficients)
        templogPY <- (-N*(1-b)/2)*log(pi) + log(b^((N*b-1)/2)) +
          (lgamma((N-tempP)/2)) - (lgamma((N*b-tempP)/2)) +
          ((-N*(1-b)/2)*log(tempSSResid))

        result.mat$marginalvec[which(result.mat$modelvec==tempmod&
                                       result.mat$hetvec==FALSE)] <- templogPY
        all.coefs[[which(result.mat$modelvec==tempmod&result.mat$schemevec=="None")]] <- tempfit$coefficients
        all.vars[[which(result.mat$modelvec==tempmod&result.mat$schemevec=="None")]] <- summary(tempfit)$sigma^2
      }

      #Heteroscedastic models with global fixed effects
      if(grepl("group",tempmod)==FALSE & het[m]==TRUE){

        index <- 1

        for(i in 1:length(d_groups)){
          for(j in 1:ncol(d_groups[[i]])){
            {
              #Create the distinct group effect
              MM <- model.matrix(tempfit)
              tempdf <- group.dfs[[index]]

              group <- as.numeric(as.numeric(lgf)%in%d_groups[[i]][,j])
              tempscheme1 <- paste(unique(lgf[group==1]),collapse=",")
              tempscheme0 <- paste(unique(lgf[group==0]),collapse=",")
              tempscheme <- paste0("{",paste(c(tempscheme1,tempscheme0),collapse="}{"),"}")

              level1 <- tempscheme1
              y1 <- tempdf[tempdf$group==level1,which(names(tempdf)==response)]
              n1 <- sum(tempdf$group==level1)
              level2 <- tempscheme0
              y2 <- tempdf[tempdf$group==level2,which(names(tempdf)==response)]
              n2 <- sum(tempdf$group==level2)

              tempP <- ncol(MM)

              f3 =function(arg){
                why=y
                bee=1
                gam1=arg[1]
                gam2=arg[2]
                Sigmad=rep(gam1, N)
                Sigmad[group==0]=gam2
                Sigma=diag(Sigmad)
                d=det(bee*t(MM)%*%Sigma%*%MM)^-.5
                s=(gam1^((bee*n1/2)-1))*(gam2^((bee*n2/2)-1))
                e1=bee*t(why)%*%Sigma%*%why
                e2=bee*t(why)%*%Sigma%*%MM%*%solve(t(MM)%*%Sigma%*%MM)%*%t(MM)%*%Sigma%*%why
                e=exp(-.5*as.numeric(e1-e2))
                return(d*s*e)
              }
              f3b=function(arg){
                why=y
                bee=b
                gam1=arg[1]
                gam2=arg[2]
                Sigmad=rep(gam1, N)
                Sigmad[group==0]=gam2
                Sigma=diag(Sigmad)
                d=det(bee*t(MM)%*%Sigma%*%MM)^-.5
                s=(gam1^((bee*n1/2)-1))*(gam2^((bee*n2/2)-1))
                e1=bee*t(why)%*%Sigma%*%why
                e2=bee*t(why)%*%Sigma%*%MM%*%solve(t(MM)%*%Sigma%*%MM)%*%t(MM)%*%Sigma%*%why
                e=exp(-.5*as.numeric(e1-e2))
                return(d*s*e)
              }

              lvf3 =function(arg){
                phi1=arg[1]
                phi2=arg[2]
                gam1=exp(-phi1)
                gam2=exp(-phi2)
                J=exp(-(phi1+phi2))
                why=y
                bee=1
                Sigmad=rep(gam1, N)
                Sigmad[group==0]=gam2
                Sigma=diag(Sigmad)
                d=det(bee*t(MM)%*%Sigma%*%MM)^-.5
                s=(gam1^((bee*n1/2)-1))*(gam2^((bee*n2/2)-1))
                e1=bee*t(why)%*%Sigma%*%why
                e2=bee*t(why)%*%Sigma%*%MM%*%solve(t(MM)%*%Sigma%*%MM)%*%t(MM)%*%Sigma%*%why
                e=exp(-.5*as.numeric(e1-e2))
                return(d*s*e*J)
              }
              lvf3b=function(arg){
                phi1=arg[1]
                phi2=arg[2]
                gam1=exp(-phi1)
                gam2=exp(-phi2)
                J=exp(-(phi1+phi2))
                why=y
                bee=b
                Sigmad=rep(gam1, N)
                Sigmad[group==0]=gam2
                Sigma=diag(Sigmad)
                d=det(bee*t(MM)%*%Sigma%*%MM)^-.5
                s=(gam1^((bee*n1/2)-1))*(gam2^((bee*n2/2)-1))
                e1=bee*t(why)%*%Sigma%*%why
                e2=bee*t(why)%*%Sigma%*%MM%*%solve(t(MM)%*%Sigma%*%MM)%*%t(MM)%*%Sigma%*%why
                e=exp(-.5*as.numeric(e1-e2))
                return(d*s*e*J)
              }

              lvlogf3 =function(arg){
                phi1=arg[1]
                phi2=arg[2]
                gam1=exp(-phi1)
                gam2=exp(-phi2)
                J=exp(-(phi1+phi2))
                why=y
                bee=1
                Sigmad=rep(gam1, N)
                Sigmad[group==0]=gam2
                Sigma=diag(Sigmad)
                d=det(bee*t(MM)%*%Sigma%*%MM)^-.5
                s=(gam1^((bee*n1/2)-1))*(gam2^((bee*n2/2)-1))
                e1=bee*t(why)%*%Sigma%*%why
                e2=bee*t(why)%*%Sigma%*%MM%*%solve(t(MM)%*%Sigma%*%MM)%*%t(MM)%*%Sigma%*%why
                e=exp(-.5*as.numeric(e1-e2))
                return(log(d)+log(s)+(-.5*as.numeric(e1-e2))+log(J))
              }
              lvlogf3b=function(arg){
                phi1=arg[1]
                phi2=arg[2]
                gam1=exp(-phi1)
                gam2=exp(-phi2)
                J=exp(-(phi1+phi2))
                why=y
                bee=b
                Sigmad=rep(gam1, N)
                Sigmad[group==0]=gam2
                Sigma=diag(Sigmad)
                d=det(bee*t(MM)%*%Sigma%*%MM)^-.5
                s=(gam1^((bee*n1/2)-1))*(gam2^((bee*n2/2)-1))
                e1=bee*t(why)%*%Sigma%*%why
                e2=bee*t(why)%*%Sigma%*%MM%*%solve(t(MM)%*%Sigma%*%MM)%*%t(MM)%*%Sigma%*%why
                e=exp(-.5*as.numeric(e1-e2))
                return(log(d)+log(s)+(-.5*as.numeric(e1-e2))+log(J))
              }

              gam1hat <- sum(fitted.models[[m]]$residuals^2)/(N-tempP)
              gam2hat <- gam1hat
              lvhat <- log(c(gam1hat, gam2hat))

              mode3=optim(lvhat,lvlogf3,
                          control=list(warn.1d.NelderMead=FALSE, fnscale=-1),
                          method="Nelder-Mead")$par

              mode3b=optim(lvhat,lvlogf3b,
                           control=list(warn.1d.NelderMead=FALSE, fnscale=-1),
                           method="Nelder-Mead")$par

              H3=det(-1*optim(mode3,lvlogf3,
                              control=list(warn.1d.NelderMead=FALSE, fnscale=-1),hessian=TRUE,
                              method="Nelder-Mead")$hessian)^-.5
              H3b=det(-1*optim(mode3b,lvlogf3b,
                               control=list(warn.1d.NelderMead=FALSE, fnscale=-1),hessian=TRUE,
                               method="Nelder-Mead")$hessian)^-.5

              templogPY <- log(((2*pi)^(-(N-tempP)/2)))-
                log(((2*pi)^(-(N*b-tempP)/2)))+
                # log(H3*(lvf3(mode3))/(H3b*lvf3b(mode3b)))
                log(H3) + lvlogf3(mode3) - log(H3b) - lvlogf3b(mode3b)

              result.mat$marginalvec[which(result.mat$modelvec==tempmod&
                                             result.mat$hetvec==TRUE)[index]] <- templogPY

              all.coefs[[which(result.mat$modelvec==tempmod&
                                 result.mat$schemevec==tempscheme&
                                 result.mat$hetvec==TRUE)]] <- tempfit$coefficients
              all.vars[[which(result.mat$modelvec==tempmod&
                                result.mat$schemevec==tempscheme&
                                result.mat$hetvec==TRUE)]] <- exp(c(mode3[2], mode3[1]))
              index=index+1
            }
          }
        }
      }

      #Homoscedastic models, group-based fixed effects
      if(grepl("group",tempmod)==TRUE){

        index <- 1

        result.mat.rows <- which(result.mat$modelvec==tempmod & result.mat$hetvec==FALSE)

        for(i in 1:length(d_groups)){
          for(j in 1:ncol(d_groups[[i]])){
            {
              group <- as.numeric(as.numeric(lgf)%in%d_groups[[i]][,j])
              tempscheme1 <- paste(unique(lgf[group==1]),collapse=",")
              tempscheme0 <- paste(unique(lgf[group==0]),collapse=",")
              tempscheme <- paste0("{",paste(c(tempscheme1,tempscheme0),collapse="}{"),"}")

              tempfit.scheme <- tempfit[[index]]
              tempSSResid <- sum(tempfit.scheme$residuals^2)
              tempP <- length(tempfit.scheme$coefficients)
              templogPY <- (-N*(1-b)/2)*log(pi)+((N*b-1)/2)*log(b)+
                (lgamma((N-tempP)/2))-(lgamma((N*b-tempP)/2))+
                ((-N*(1-b)/2)*log(tempSSResid))

              result.mat$marginalvec[result.mat.rows[index]] <- templogPY

              index <- index+1

              all.coefs[[which(result.mat$modelvec==tempmod&
                                 result.mat$schemevec==tempscheme&
                                 result.mat$hetvec==FALSE)]] <- tempfit.scheme$coefficients
              all.vars[[which(result.mat$modelvec==tempmod&
                                result.mat$schemevec==tempscheme&
                                result.mat$hetvec==FALSE)]] <- summary(tempfit.scheme)$sigma^2
            }
          }
        }

      }

      #Heteroscedastic models, group-based fixed effects
      if(grepl("group",tempmod)==TRUE & het[m]==TRUE){

        index <- 1

        result.mat.rows <- which(result.mat$modelvec==tempmod & result.mat$hetvec==TRUE)

        for(i in 1:length(d_groups)){
          for(j in 1:ncol(d_groups[[i]])){
            {
              #Create the distinct group effect
              MM <- model.matrix(tempfit[[index]])
              tempdf <- group.dfs[[index]]

              group <- as.numeric(as.numeric(lgf)%in%d_groups[[i]][,j])
              tempscheme1 <- paste(unique(lgf[group==1]),collapse=",")
              tempscheme0 <- paste(unique(lgf[group==0]),collapse=",")
              tempscheme <- paste0("{",paste(c(tempscheme1,tempscheme0),collapse="}{"),"}")

              level1 <- tempscheme1
              y1 <- tempdf[tempdf$group==level1,which(names(tempdf)==response)]
              n1 <- sum(tempdf$group==level1)
              level2 <- tempscheme0
              y2 <- tempdf[tempdf$group==level2,which(names(tempdf)==response)]
              n2 <- sum(tempdf$group==level2)

              tempP <- ncol(MM)

              f3 =function(arg){
                why=y
                bee=1
                gam1=arg[1]
                gam2=arg[2]
                Sigmad=rep(gam1, N)
                Sigmad[group==0]=gam2
                Sigma=diag(Sigmad)
                d=det(bee*t(MM)%*%Sigma%*%MM)^-.5
                s=(gam1^((bee*n1/2)-1))*(gam2^((bee*n2/2)-1))
                e1=bee*t(why)%*%Sigma%*%why
                e2=bee*t(why)%*%Sigma%*%MM%*%solve(t(MM)%*%Sigma%*%MM)%*%t(MM)%*%Sigma%*%why
                e=exp(-.5*as.numeric(e1-e2))
                return(d*s*e)
              }
              f3b=function(arg){
                why=y
                bee=b
                gam1=arg[1]
                gam2=arg[2]
                Sigmad=rep(gam1, N)
                Sigmad[group==0]=gam2
                Sigma=diag(Sigmad)
                d=det(bee*t(MM)%*%Sigma%*%MM)^-.5
                s=(gam1^((bee*n1/2)-1))*(gam2^((bee*n2/2)-1))
                e1=bee*t(why)%*%Sigma%*%why
                e2=bee*t(why)%*%Sigma%*%MM%*%solve(t(MM)%*%Sigma%*%MM)%*%t(MM)%*%Sigma%*%why
                e=exp(-.5*as.numeric(e1-e2))
                return(d*s*e)
              }

              lvf3 =function(arg){
                phi1=arg[1]
                phi2=arg[2]
                gam1=exp(-phi1)
                gam2=exp(-phi2)
                J=exp(-(phi1+phi2))
                why=y
                bee=1
                Sigmad=rep(gam1, N)
                Sigmad[group==0]=gam2
                Sigma=diag(Sigmad)
                d=det(bee*t(MM)%*%Sigma%*%MM)^-.5
                s=(gam1^((bee*n1/2)-1))*(gam2^((bee*n2/2)-1))
                e1=bee*t(why)%*%Sigma%*%why
                e2=bee*t(why)%*%Sigma%*%MM%*%solve(t(MM)%*%Sigma%*%MM)%*%t(MM)%*%Sigma%*%why
                e=exp(-.5*as.numeric(e1-e2))
                return(d*s*e*J)
              }
              lvf3b=function(arg){
                phi1=arg[1]
                phi2=arg[2]
                gam1=exp(-phi1)
                gam2=exp(-phi2)
                J=exp(-(phi1+phi2))
                why=y
                bee=b
                Sigmad=rep(gam1, N)
                Sigmad[group==0]=gam2
                Sigma=diag(Sigmad)
                d=det(bee*t(MM)%*%Sigma%*%MM)^-.5
                s=(gam1^((bee*n1/2)-1))*(gam2^((bee*n2/2)-1))
                e1=bee*t(why)%*%Sigma%*%why
                e2=bee*t(why)%*%Sigma%*%MM%*%solve(t(MM)%*%Sigma%*%MM)%*%t(MM)%*%Sigma%*%why
                e=exp(-.5*as.numeric(e1-e2))
                return(d*s*e*J)
              }

              lvlogf3 =function(arg){
                phi1=arg[1]
                phi2=arg[2]
                gam1=exp(-phi1)
                gam2=exp(-phi2)
                # J=exp(-(phi1+phi2))
                logJ=-(phi1+phi2)
                why=y
                bee=1
                Sigmad=rep(gam1, N)
                Sigmad[group==0]=gam2
                Sigma=diag(Sigmad)
                # d=det(bee*t(MM)%*%Sigma%*%MM)^-.5
                logd=-.5*log(det(bee*t(MM)%*%Sigma%*%MM))
                # s=(gam1^((bee*n1/2)-1))*(gam2^((bee*n2/2)-1))
                logs=((bee*n1/2)-1)*log(gam1) + ((bee*n2/2)-1)*log(gam2)
                e1=bee*t(why)%*%Sigma%*%why
                e2=bee*t(why)%*%Sigma%*%MM%*%solve(t(MM)%*%Sigma%*%MM)%*%t(MM)%*%Sigma%*%why
                e=exp(-.5*as.numeric(e1-e2))
                return(logd+logs+(-.5*as.numeric(e1-e2))+logJ)
              }
              lvlogf3b=function(arg){
                phi1=arg[1]
                phi2=arg[2]
                gam1=exp(-phi1)
                gam2=exp(-phi2)
                J=exp(-(phi1+phi2))
                logJ=-(phi1+phi2)
                why=y
                bee=b
                Sigmad=rep(gam1, N)
                Sigmad[group==0]=gam2
                Sigma=diag(Sigmad)
                # d=det(bee*t(MM)%*%Sigma%*%MM)^-.5
                logd=-.5*log(det(bee*t(MM)%*%Sigma%*%MM))
                # s=(gam1^((bee*n1/2)-1))*(gam2^((bee*n2/2)-1))
                logs=((bee*n1/2)-1)*log(gam1) + ((bee*n2/2)-1)*log(gam2)
                e1=bee*t(why)%*%Sigma%*%why
                e2=bee*t(why)%*%Sigma%*%MM%*%solve(t(MM)%*%Sigma%*%MM)%*%t(MM)%*%Sigma%*%why
                e=exp(-.5*as.numeric(e1-e2))
                return(logd+logs+(-.5*as.numeric(e1-e2))+logJ)
              }

              lv1hat <- log(sum(tempfit[[index]]$residuals[tempdf$group==level1]^2)/(N-tempP))
              lv2hat <- log(sum(tempfit[[index]]$residuals[tempdf$group==level2]^2)/(N-tempP))

              mode3=optim(c(lv1hat,lv2hat),lvlogf3,
                          control=list(warn.1d.NelderMead=FALSE, fnscale=-1),
                          method="Nelder-Mead")$par

              mode3b=optim(c(lv1hat,lv2hat),lvlogf3b,
                           control=list(warn.1d.NelderMead=FALSE, fnscale=-1),
                           method="Nelder-Mead")$par

              H3=det(-1*optim(mode3,lvlogf3,
                              control=list(warn.1d.NelderMead=FALSE, fnscale=-1),hessian=TRUE,
                              method="Nelder-Mead")$hessian)^-.5
              H3b=det(-1*optim(mode3b,lvlogf3b,
                               control=list(warn.1d.NelderMead=FALSE, fnscale=-1),hessian=TRUE,
                               method="Nelder-Mead")$hessian)^-.5

              templogPY <- log(((2*pi)^(-(N-tempP)/2)))-
                log(((2*pi)^(-(N*b-tempP)/2)))+
                # log(H3*(lvf3(mode3))/(H3b*lvf3b(mode3b)))
                log(H3) + lvlogf3(mode3) - log(H3b) - lvlogf3b(mode3b)

              result.mat$marginalvec[which(result.mat$modelvec==tempmod&
                                             result.mat$hetvec==TRUE)[index]] <- templogPY


              tempfit.scheme <- tempfit[[index]]

              all.coefs[[which(result.mat$modelvec==tempmod&
                                 result.mat$schemevec==tempscheme&
                                 result.mat$hetvec==TRUE)]] <- tempfit.scheme$coefficients
              all.vars[[which(result.mat$modelvec==tempmod&
                                result.mat$schemevec==tempscheme&
                                result.mat$hetvec==TRUE)]] <- exp(c(mode3[2], mode3[1]))
              index=index+1
            }
          }
        }
      }

    }
  }

  #total sum of squares
  SST <- sum((y-mean(y))^2)
  Z <- rep(1,N)

  #Compute marginal model probabilities - Zellner-Siow mixture g-prior case

  if(prior=="zs"){
    for(tempmod in unlist(usermodels)){

      m <- which(tempmod==names(fitted.models))
      tempfit <- fitted.models[[m]]

      # y <- y-mean(y) # added 7/24/20, removed 8/12

      #done 3/23
      #Homoscedsatic models, group-based effects
      if(grepl("group",tempmod)==FALSE){

        # tempSSResid <- sum(tempfit$residuals^2)
        tempP <- sum(!is.na(tempfit$coefficients))-1
        R2A <- summary(tempfit)$r.squared

        #Mode of additive marginal distribution
        faddmode <- function(gee){
          Q <- 1-R2A
          value <- -Q*(tempP+3)*(gee^3)+
            (N-tempP-4-2*(1-R2A))*(gee^2)+
            ((N*(2-R2A)-3)*gee)+N
          return(value)
        }
        #Mode of fractional additive marginal distribution
        faddbmode=function(gee){
          Q <- 1-R2A
          value=-Q*(b^2)*(tempP+3)*(gee^3)+
            (b*(N*b-tempP-4)-2*Q)*(gee^2)+
            (N*b*(2-R2A)-3)*gee+N
          return(value)
        }

        addmode <- uniroot(faddmode,c(0,1e9),check.conv=TRUE,tol=1e-10)$root
        addbmode <- uniroot(faddbmode,c(1e-9,1e9),check.conv=TRUE,tol=1e-10)$root

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

        gamhatadd <- 1/(summary(tempfit)$sigma^2)

        templogPY <- qaddLA - qaddbLA

        result.mat$marginalvec[which(result.mat$modelvec==tempmod&
                                       result.mat$hetvec==FALSE)] <- templogPY
        all.coefs[[which(result.mat$modelvec==tempmod&result.mat$schemevec=="None")]] <- tempfit$coefficients
        all.vars[[which(result.mat$modelvec==tempmod&result.mat$schemevec=="None")]] <- summary(tempfit)$sigma^2
        all.gs[[which(result.mat$modelvec==tempmod&result.mat$schemevec=="None")]] <- addmode
      }

      #Heteroscedastic models with global fixed effects
      if(grepl("group",tempmod)==FALSE & het[m]==TRUE){

        index <- 1

        #done 3/24
        for(i in 1:length(d_groups)){
          for(j in 1:ncol(d_groups[[i]])){
            {
              #Create the distinct group effect
              MM <- column_centerer(model.matrix(tempfit))[,-1] #column centerer added 7/20
              # MM <- model.matrix(tempfit)[,-1] # changed back 8/12
              tempdf <- group.dfs[[index]]
              tempP <- ncol(MM) #-1 already done
              R2A <- summary(tempfit)$r.squared
              # gamhat <- 1/(summary(tempfit)$sigma^2)
              gamhat <- gamhatadd
              lvhat <- log(1/gamhat)
              level1 <- levels(tempdf$group)[1]
              y1 <- tempdf[tempdf$group==level1,which(names(tempdf)==response)]
              n1 <- sum(tempdf$group==level1)
              level2 <- levels(tempdf$group)[2]
              y2 <- tempdf[tempdf$group==level2,which(names(tempdf)==response)]
              n2 <- sum(tempdf$group==level2)
              scheme <- as.numeric(tempdf$group==level1)
              tempscheme1 <- levels(tempdf$group)[1]
              tempscheme0 <- levels(tempdf$group)[2]

              lvlogf3 =function(arg){
                gee=arg[1]
                u=arg[2]
                v=arg[3]
                gam1=exp(-u)
                gam2=exp(-v)
                J=exp(-(u+v))
                logJ <- -(u+v)
                gamvec=rep(NA,N)
                gamvec[scheme==1]=gam1
                gamvec[scheme==0]=gam2
                S=diag(gamvec)
                Z_S=S%*%Z%*%solve(t(Z)%*%S%*%Z)%*%t(Z)%*%S

                value=(n1/2-1)*log(gam1)+(n2/2-1)*log(gam2)-((tempP/2)*log(gee))+
                  .5*log(det(t(MM)%*%S%*%MM))-
                  .5*log(det(t(Z)%*%S%*%Z))-
                  .5*log(det(((gee+1)/gee)*t(MM)%*%S%*%MM-t(MM)%*%Z_S%*%MM))-
                  .5*(t(y)%*%S%*%y-t(y)%*%Z_S%*%y-
                        t(y)%*%(S-Z_S)%*%MM%*%solve(
                          ((gee+1)/gee)*t(MM)%*%S%*%MM-
                            t(MM)%*%Z_S%*%MM)%*%t(MM)%*%(S-Z_S)%*%y)-
                  1.5*log(gee)-(N/(2*gee))+logJ
                return(value)
              }
              lvlogf3b=function(arg){
                gee=arg[1]
                u=arg[2]
                v=arg[3]
                gam1=exp(-u)
                gam2=exp(-v)
                J=exp(-(u+v))
                logJ <- -(u+v)
                gamvec=rep(NA,N)
                gamvec[scheme==1]=gam1
                gamvec[scheme==0]=gam2
                S=diag(gamvec)
                Z_S=S%*%Z%*%solve(t(Z)%*%S%*%Z)%*%t(Z)%*%S
                bg=b*gee

                value=(n1*b/2-1)*log(gam1)+(n2*b/2-1)*log(gam2)-((tempP/2)*log(gee))+
                  .5*log(det(t(MM)%*%S%*%MM))-
                  .5*log(det(t(Z)%*%S%*%Z))-
                  .5*log(det(((bg+1)/bg)*t(MM)%*%S%*%MM-t(MM)%*%Z_S%*%MM))-
                  .5*b*(t(y)%*%S%*%y-t(y)%*%Z_S%*%y-
                          t(y)%*%(S-Z_S)%*%MM%*%solve(
                            ((bg+1)/bg)*t(MM)%*%S%*%MM-
                              t(MM)%*%Z_S%*%MM)%*%t(MM)%*%(S-Z_S)%*%y)-
                  1.5*log(gee)-(N/(2*gee))+logJ
                return(value)
              }

              lvmode3 <- optim(c(N,lvhat,lvhat),
                               lvlogf3,
                               control=list(warn.1d.NelderMead=FALSE, fnscale=-1),
                               method="Nelder-Mead")$par
              lvmode3b <- optim(c(N,lvhat,lvhat),
                                lvlogf3b,
                                control=list(warn.1d.NelderMead=FALSE, fnscale=-1),
                                method="Nelder-Mead")$par

              H3 <- (-1*det(numDeriv::hessian(lvlogf3,lvmode3)))^-.5
              logH3 <- -.5*determinant(-1*numDeriv::hessian(lvlogf3,lvmode3), logarithm=TRUE)$modulus
              H3b <- (-1*det(numDeriv::hessian(lvlogf3b,lvmode3b)))^-.5
              logH3b <- -.5*determinant(-1*numDeriv::hessian(lvlogf3b,lvmode3b), logarithm=TRUE)$modulus

              q3LA <- (-(N-1)/2)*log(2*pi) + .5*log(N/2) - lgamma(.5) +
                1.5*log(2*pi) + logH3 + lvlogf3(lvmode3)
              q3bLA <- (-(N*b-1)/2)*log(2*pi) + (-(tempP+1)/2)*log(b) +
                .5*log(N/2) - lgamma(.5) + 1.5*log(2*pi) +
                logH3b + lvlogf3b(lvmode3b)

              templogPY <- q3LA - q3bLA

              result.mat$marginalvec[which(result.mat$modelvec==tempmod&
                                             result.mat$hetvec==TRUE)[index]] <- templogPY

              tempscheme <- names(group.dfs[index])

              all.coefs[[which(result.mat$modelvec==tempmod&
                                 result.mat$schemevec==tempscheme&
                                 result.mat$hetvec==TRUE)]] <- tempfit$coefficients
              all.vars[[which(result.mat$modelvec==tempmod&
                                result.mat$schemevec==tempscheme&
                                result.mat$hetvec==TRUE)]] <- exp(c(lvmode3[2], lvmode3[1]))
              all.gs[[which(result.mat$modelvec==tempmod&
                              result.mat$schemevec==tempscheme&
                              result.mat$hetvec==TRUE)]] <- lvmode3[1]
              index=index+1
            }
          }
        }
      }

      #Group-based fixed effects
      if(grepl("group",tempmod)==TRUE){

        index <- 1

        result.mat.rows <- which(result.mat$modelvec==tempmod & result.mat$hetvec==FALSE)

        #done 3/23
        #homoscedastic group-based effects
        for(i in 1:length(d_groups)){
          for(j in 1:ncol(d_groups[[i]])){
            {
              {
                tempscheme <- gsub("scheme=","",names(tempfit)[index])

                tempfit.scheme <- tempfit[[index]]
                tempSSResid <- sum(tempfit.scheme$residuals^2)
                tempP <- sum(!is.na(tempfit.scheme$coefficients))-1 #-1 added July 17

                R2A=1-(tempSSResid/SST) #correct
                Q <- 1-R2A

                #Mode of additive marginal distribution
                faddmode <- function(gee){
                  Q <- 1-R2A
                  value <- -Q*(tempP+3)*(gee^3)+
                    (N-tempP-4-2*(1-R2A))*(gee^2)+
                    ((N*(2-R2A)-3)*gee)+N
                  return(value)
                }
                #Mode of fractional additive marginal distribution
                faddbmode=function(gee){
                  Q <- 1-R2A
                  value=-Q*(b^2)*(tempP+3)*(gee^3)+
                    (b*(N*b-tempP-4)-2*Q)*(gee^2)+
                    (N*b*(2-R2A)-3)*gee+N
                  return(value)
                }

                addmode <- uniroot(faddmode,c(0,1e9),check.conv=TRUE,tol=1e-10)$root
                addbmode <- uniroot(faddbmode,c(1e-9,1e9),check.conv=TRUE,tol=1e-10)$root

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

                result.mat$marginalvec[result.mat.rows[index]] <- templogPY

                all.coefs[[which(result.mat$modelvec==tempmod&
                                   result.mat$schemevec==tempscheme&
                                   result.mat$hetvec==FALSE)]] <- tempfit.scheme$coefficients
                all.vars[[which(result.mat$modelvec==tempmod&
                                  result.mat$schemevec==tempscheme&
                                  result.mat$hetvec==FALSE)]] <- summary(tempfit.scheme)$sigma^2
                all.gs[[which(result.mat$modelvec==tempmod&
                                result.mat$schemevec==tempscheme&
                                result.mat$hetvec==FALSE)]] <- addmode

                index <- index+1

              }
            }
          }

        }

        #Heteroscedastic models, group-based fixed effects
        if(grepl("group",tempmod)==TRUE & het[m]==TRUE){

          index <- 1

          result.mat.rows <- which(result.mat$modelvec==tempmod & result.mat$hetvec==TRUE)

          for(i in 1:length(d_groups)){
            for(j in 1:ncol(d_groups[[i]])){
              {
                tempscheme <- gsub("scheme=","",names(tempfit)[index])

                tempfit.scheme <- tempfit[[index]]
                tempSSResid <- sum(tempfit.scheme$residuals^2)
                # tempP <- sum(!is.na(tempfit.scheme$coefficients))-1

                MM <- model.matrix(tempfit.scheme)[,!is.na(tempfit.scheme$coefficients)]
                MM <- column_centerer(MM)[,-1] # removed 8/12
                tempP <- ncol(MM) # added 8/12
                tempdf <- group.dfs[[index]]
                R2A <- summary(tempfit.scheme)$r.squared
                # gamhat <- 1/(summary(tempfit.scheme)$sigma^2)
                gamhat <- gamhatadd
                lvhat <- log(1/gamhat)
                level1 <- levels(tempdf$group)[1]
                y1 <- tempdf[tempdf$group==level1,which(names(tempdf)==response)]#-ybar #added 8/12
                n1 <- sum(tempdf$group==level1)
                level2 <- levels(tempdf$group)[2]
                y2 <- tempdf[tempdf$group==level2,which(names(tempdf)==response)]#-ybar
                n2 <- sum(tempdf$group==level2)
                scheme <- as.numeric(tempdf$group==level1)

                lvlogf3  <- function(arg){
                  gee <- arg[1]
                  u <- arg[2]
                  v <- arg[3]
                  gam1 <- exp(-u)
                  gam2 <- exp(-v)
                  J <- exp(-(u+v))
                  logJ <- -(u+v)
                  gamvec <- rep(NA,N)
                  gamvec[scheme==1] <- gam1
                  gamvec[scheme==0] <- gam2
                  S <- diag(gamvec)
                  Z_S <- S%*%Z%*%solve(t(Z)%*%S%*%Z)%*%t(Z)%*%S

                  value <- (n1/2-1)*(-u)+(n2/2-1)*(-v)-((tempP/2)*log(gee))+
                    .5*determinant(t(MM)%*%S%*%MM, logarithm=TRUE)$modulus -
                    .5*determinant(t(Z)%*%S%*%Z, logarithm=TRUE)$modulus -
                    .5*determinant(((gee+1)/gee)*t(MM)%*%S%*%MM-t(MM)%*%Z_S%*%MM,
                                   logarithm=TRUE)$modulus -
                    .5*(t(y)%*%S%*%y-t(y)%*%Z_S%*%y-
                          t(y)%*%(S-Z_S)%*%MM%*%solve(
                            ((gee+1)/gee)*t(MM)%*%S%*%MM-
                              t(MM)%*%Z_S%*%MM)%*%t(MM)%*%(S-Z_S)%*%y)-
                    1.5*log(gee)-(N/(2*gee))+logJ
                  return(value)
                }
                lvlogf3b <- function(arg){
                  gee <- arg[1]
                  u <- arg[2]
                  v <- arg[3]
                  gam1 <- exp(-u)
                  gam2 <- exp(-v)
                  J <- exp(-(u+v))
                  logJ <- -(u+v)
                  gamvec <- rep(NA,N)
                  gamvec[scheme==1] <- gam1
                  gamvec[scheme==0] <- gam2
                  S <- diag(gamvec)
                  Z_S <- S%*%Z%*%solve(t(Z)%*%S%*%Z)%*%t(Z)%*%S
                  bg <- b*gee

                  value <- (n1*b/2-1)*(-u)+(n2*b/2-1)*(-v)-((tempP/2)*log(gee))+
                    .5*determinant(t(MM)%*%S%*%MM, logarithm=TRUE)$modulus -
                    .5*determinant(t(Z)%*%S%*%Z, logarithm=TRUE)$modulus -
                    .5*determinant(((bg+1)/bg)*t(MM)%*%S%*%MM-t(MM)%*%Z_S%*%MM,
                                   logarithm=TRUE)$modulus -
                    .5*b*(t(y)%*%S%*%y-t(y)%*%Z_S%*%y-
                            t(y)%*%(S-Z_S)%*%MM%*%solve(
                              ((bg+1)/bg)*t(MM)%*%S%*%MM-
                                t(MM)%*%Z_S%*%MM)%*%t(MM)%*%(S-Z_S)%*%y)-
                    1.5*log(gee)-(N/(2*gee))+logJ
                  return(value)
                }

                lvmode3 <- optim(c(N,lvhat,lvhat),
                                 lvlogf3,
                                 control=list(warn.1d.NelderMead=FALSE, fnscale=-1),
                                 method="Nelder-Mead")$par
                lvmode3b <- optim(c(N,lvhat,lvhat),
                                  lvlogf3b,
                                  control=list(warn.1d.NelderMead=FALSE, fnscale=-1),
                                  method="Nelder-Mead")$par

                H3 <- (-1*det(numDeriv::hessian(lvlogf3,lvmode3)))^-.5
                logH3 <- -.5*determinant(-1*numDeriv::hessian(lvlogf3,lvmode3), logarithm=TRUE)$modulus
                H3b <- (-1*det(numDeriv::hessian(lvlogf3b,lvmode3b)))^-.5
                logH3b <- -.5*determinant(-1*numDeriv::hessian(lvlogf3b,lvmode3b), logarithm=TRUE)$modulus

                q3LA <- (-(N-1)/2)*log(2*pi) + .5*log(N/2) - lgamma(.5) +
                  1.5*log(2*pi) + logH3 + lvlogf3(lvmode3)
                q3bLA <- (-(N*b-1)/2)*log(2*pi) + (-(tempP+1)/2)*log(b) +
                  .5*log(N/2) - lgamma(.5) + 1.5*log(2*pi) +
                  logH3b + lvlogf3b(lvmode3b)

                templogPY <- q3LA - q3bLA

                result.mat$marginalvec[which(result.mat$modelvec==tempmod&
                                               result.mat$hetvec==TRUE)[index]] <- templogPY
                tempscheme <- names(group.dfs[index])

                all.coefs[[which(result.mat$modelvec==tempmod&
                                   result.mat$schemevec==tempscheme&
                                   result.mat$hetvec==TRUE)]] <- tempfit.scheme$coefficients
                all.vars[[which(result.mat$modelvec==tempmod&
                                  result.mat$schemevec==tempscheme&
                                  result.mat$hetvec==TRUE)]] <- exp(c(lvmode3[2], lvmode3[1]))
                all.gs[[which(result.mat$modelvec==tempmod&
                                result.mat$schemevec==tempscheme&
                                result.mat$hetvec==TRUE)]] <- lvmode3[1]

                index=index+1
              }
            }
          }
        }
      }

    }
  }


  #remove fix model if necessary
  if(fix==1){
    delete <- which(as.character(result.mat$modelvec)==paste0(response,"~group"))
    result.mat <- result.mat[-delete,]
    row.names(result.mat) <- 1:nrow(result.mat)
    result.mat$classvec <- droplevels(result.mat$classvec)
    result.mat$classvec <- as.factor(rep(1:length(levels(result.mat$classvec)), times=table(result.mat$classvec)))
  }

  modevs <- exp(result.mat$marginalvec+abs(max(result.mat$marginalvec)))
  # modevs <- exp(result.mat$marginalvec-(max(result.mat$marginalvec)))

  #Create uniform prior by class
  modpriors <- rep(NA,length(modevs))
  modprobs <- rep(NA,length(modevs))

  modpriors[result.mat$scheme=="None"] <- 1/length(levels(as.factor(unique(result.mat$classvec))))
  modpriors[result.mat$scheme!="None"] <- 1/(ngroups*length(levels(as.factor(unique(result.mat$classvec)))))

  for(i in 1:length(modevs)){
    modprobs[i]=(modevs[i]*modpriors[i])/(modevs%*%modpriors)
  }

  result.mat$modpriors <- modpriors
  result.mat$modprobs <- round(modprobs,8)

  #Aggregate probabilities by grouping
  #scheme to obtain class probabilities
  classprobs <- rep(NA,length(unique(result.mat$classvec)))
  names(classprobs) <- levels(result.mat$classvec)
  for(c in as.numeric(as.character(unique(result.mat$classvec)))){
    classprobs[c] <- sum(modprobs[as.numeric(as.character(result.mat$classvec))==c])
  }

  result.mat$mle.index <- seq(1:nrow(result.mat))
  sorted <- result.mat[order(-modprobs),]
  mle.index <- sorted$mle.index
  sorted$cumuprob <- cumsum(sorted$modprobs)

  sorted <- sorted[,-c(3,8)]

  variance <- rep(NA, length(result.mat$hetvec))
  variance[sorted$hetvec==TRUE] <- c("Heterosk")
  variance[sorted$hetvec==FALSE] <- c("Homosk")

  sorted$hetvec <- variance

  sorted$df.Index <- rep(NA, nrow(result.mat)) #nrow(result.mat)
  for(i in 1:nrow(result.mat)){
    if(sorted$schemevec[i]!="None"){
      sorted$df.Index[i] <- which(names(group.dfs)==sorted$schemevec[i])
    }
    if(sorted$schemevec[i]=="None"){
      sorted$df.Index[i] <- which(names(group.dfs)==sorted$modelvec[i])
    }
  }

  sorted$mle.index <- mle.index
  sorted$model.index <- 1:nrow(sorted)
  sorted$class <- paste0(sorted$modelvec, ", ", sorted$hetvec)
  names(sorted) <- c("Model", "Scheme", "Variance", "logFlik",
                     "Mod.Prior", "Fmodprob", "Cumulative", "dataf.Index", "mle.index", "Model.Index", "Class")
  row.names(sorted) <- NULL

  classprobs <- rep(NA, length(levels(as.factor(sorted$Class))))
  names(classprobs) <- levels(as.factor(sorted$Class))
  for(class in levels(as.factor(sorted$Class))){
    tempprobs <- sum(as.numeric(as.character(sorted$Fmodprob[sorted$Class==class])))
    classprobs[which(names(classprobs)==class)] <- tempprobs
  }

  schemeprobs <- rep(NA, length(levels(as.factor(sorted$Scheme))))
  names(schemeprobs) <- levels(as.factor(sorted$Scheme))
  for(scheme in levels(as.factor(sorted$Scheme))){
    tempprobs <- sum(as.numeric(as.character(sorted$Fmodprob[sorted$Scheme==scheme])))
    schemeprobs[which(names(schemeprobs)==scheme)] <- tempprobs
  }

  schemeprobs <- data.frame("Scheme.Prob"=sort(schemeprobs, decreasing=TRUE))
  classprobs <- data.frame("Class.Prob"=sort(classprobs, decreasing=TRUE))
  # print(c(lb, ub))
  if(prior=="flat"){
    out <- list("results"=sorted, "group.datafs"=group.dfs, "class.Probs"=round(classprobs, 8),
                "scheme.Probs"=schemeprobs,
                "coefficients"=all.coefs, "variances"=all.vars, "b"=b)
  }
  if(prior=="zs"){
    out <- list("results"=sorted, "group.datafs"=group.dfs, "class.Probs"=round(classprobs, 8),
                "scheme.Probs"=schemeprobs,
                "coefficients"=all.coefs, "variances"=all.vars, "gs"=all.gs, "b"=b)
  }
  class(out) <- "slgf"
  return(out)
}


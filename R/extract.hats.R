#' Obtain the concentrated maximum likelihood estimators from a specified model.
#'
#' @importFrom stats var
#'
#' @description \code{extract.hats} Returns the concentrated maximum likelihood estimators from a specified model.
#'
#' @param slgf.obj output from \code{ms.slgf}
#'
#' @param model.index the model index of the model for which estimates are desired.
#'
#' @return \code{extract.hats} returns a list with the following elements:\cr
#' 1) \code{model}, the model desired\cr
#' 2) \code{scheme}, the scheme associated with the model desired\cr
#' 3) \code{coef}, the regression coefficients associated with the model desired\cr
#' 4) \code{sigsq}, the error variance(s) assicoated with the model desired\cr
#' 5) \code{g}, the g estimate, if \code{prior="zs"}\cr
#'
#' @export
#' @examples
#' # Obtain the concentrated maximum likelihood estimates
#' # for the second-most probable model.
#'
#' library(numDeriv)
#'
#' set.seed(314159)
#' test.data <- data.frame("y"=c(rnorm(10,0,1), rnorm(10,3,1), rnorm(10,5,3)),
#'                         "x"=c(rep("A",10), rep("B",10), rep("C",10)))
#' test.models <- list("y~1", "y~x", "y~group")
#' test.models
#' test.out <- ms.slgf(dataf=test.data, response="y", lgf="x",
#'                     usermodels=test.models,
#'                     prior="flat", het=c(1,1,1), min.levels=1)
#' extract.hats(test.out, 2)
extract.hats <- function(slgf.obj, model.index=NULL){
  result <- slgf.obj$results
  if(is.null(model.index)){
    model.index <- which(result$Model==model & result$Scheme==scheme & result$Variance==var)
  }
  if(result$Variance[model.index]=="Homosk"){
    model <- as.character(result$Model[model.index])
    scheme <- as.character(result$Scheme[model.index])
    mle.index <- result$mle.index[model.index]
    mle.coef <- unlist(unname(slgf.obj$coefficients[mle.index]))
    mle.var <- unlist(unname(slgf.obj$variances[mle.index]))
    return.list <- list(model=model, scheme=scheme, coef=mle.coef, sigsq=mle.var)
    names(return.list)[1] <- "model"
    return(return.list)
    # break
  }
  if(result$Variance[model.index]=="Heterosk"){
    model <- as.character(result$Model[model.index])
    scheme <- as.character(result$Scheme[model.index])
    scheme1 <- paste0(strsplit(scheme, "}")[[1]][1],"}")
    scheme2 <- paste0(strsplit(scheme, "}")[[1]][2],"}")
    mle.index <- result$mle.index[model.index]
    mle.coef <- unlist(unname(slgf.obj$coefficients[mle.index]))
    mle.var1 <- unlist(unname(slgf.obj$variances[mle.index][[1]][1]))
    mle.var2 <- unlist(unname(slgf.obj$variances[mle.index][[1]][2]))

    return.list <- list(model=model, scheme=scheme, coef=mle.coef,
                        sigsq.1=mle.var1, sigsq.2=mle.var2)
    names(return.list)[4] <- paste0("sigsq",".",scheme1)
    names(return.list)[5] <- paste0("sigsq",".",scheme2)

    return(return.list)
    # break
  }

}

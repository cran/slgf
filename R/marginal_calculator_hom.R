.marginal_calculator_hom <- function(tempfit, N, b){
  tempSSResid <- sum(tempfit$residuals^2)
  tempP <- length(tempfit$coefficients)
  templogPY <- (-N*(1-b)/2)*log(pi) +
               ((N*b-1)/2)*log(b) +
               (lgamma((N-tempP)/2))-(lgamma((N*b-tempP)/2))+
               ((-N*(1-b)/2)*log(tempSSResid))
  return(templogPY)
}


## ---- include = FALSE---------------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)

## ----setup, include=FALSE, warning=FALSE--------------------------------------
knitr::opts_chunk$set(echo = TRUE)
#knitr::opts_chunk$set(fig.width=8, fig.height=6)
library(slgf)
library(rcrossref)
library(captioner)
library(formatR)
library(numDeriv)
# setwd("..")

## -----------------------------------------------------------------------------
data(smell)

## ---- fig.cap="O'Brien and Heft (1995) studied olfactory function by age ($y$-axis), where age was divided into five categories ($x$-axis). The data suggests potential heteroscedasticity, with latent grouping structure {1,2,3}{4,5}.", warning=FALSE, message=F, tidy=TRUE, cache=FALSE, tidy.opts=list(width.cutoff=60), fig.align="center", fig.pos="h", fig.width=5, fig.height=3----
boxplot(smell$olf ~ smell$agecat, pch=16, 
        xlab="Age Category", ylab="Olfactory Score", xaxt="n", 
        col=c("gray30","gray30","gray30","white","white"),
        border=c("black","black","black","red","red"),
        main="O'Brien and Heft (1995) Smell Data")
axis(1, labels=c("1\n(n=38)","2\n(n=35)","3\n(n=21)","4\n(n=43)","5\n(n=42)"),at=1:5, lwd.tick=0)

## ---- cache=FALSE-------------------------------------------------------------
smell.models <- list("olf~1", "olf~agecat", "olf~group")
smell.het <- c(0, 0, 1)

## ---- cache=FALSE-------------------------------------------------------------
smell.out <- ms.slgf(dataf=smell, response="olf", lgf="agecat", usermodels=smell.models, 
                     prior="flat", het=smell.het, min.levels=1)

## -----------------------------------------------------------------------------
smell.out$results[1:5, c(1:3,6,11)]

## ---- cache=FALSE-------------------------------------------------------------
smell.hats <- extract.hats(smell.out, model.index=1)
print(smell.hats$`sigsq.{4,5}`)
print(smell.hats$`sigsq.{1,2,3}`)

## ---- cache=FALSE-------------------------------------------------------------
smell.out <- ms.slgf(dataf=smell, response="olf", lgf="agecat", 
                     usermodels=smell.models, prior="flat", 
                     het=c(1,1,1), min.levels=1)

smell.out$results[1:5, c(1:3,6,11)]

## -----------------------------------------------------------------------------
smell.out$scheme.Probs
smell.out$class.Probs

## ---- echo=FALSE--------------------------------------------------------------
data(lymphoma)
plot(c(0.6,2.5),c(min(lymphoma)-.05,max(lymphoma))+.05,col="white", ylab="Gene Expression", xlab="Tissue Source", main="Franck et. al. (2013)", xaxt="n", cex.lab=1, cex.axis=1, cex.main=1)
axis(1,labels=c("Tumor","Normal"), at=c(1,2))
lines(1:2,lymphoma[1,],lwd=2,col="gray70")
lines(1:2,lymphoma[2,],lwd=2,col="gray70")
lines(1:2,lymphoma[3,],lwd=2,col="gray70",lty=1)
lines(1:2,lymphoma[4,],lwd=2,col="gray70",lty=1)
lines(1:2,lymphoma[5,],lwd=2,col="gray70")
lines(1:2,lymphoma[6,],lwd=2,col="gray70",lty=1)
text(rep(.85,6), c(9.3278,9.5158,8.7835,8.6372,9.4581,8.7222), c("Dog 1","Dog 2","Dog 3","Dog 4","Dog 5","Dog 6"),cex=1)


## ---- cache=FALSE, echo=TRUE, message=FALSE, results='hide', warning=FALSE----
lymphoma.tall <- maketall(lymphoma)
lymphoma.tall <- data.frame("gene"=lymphoma.tall[,1], "dog"=lymphoma.tall[,2], 
                            "tissue"=lymphoma.tall[,3])

lymphoma.models <- list("gene~dog+tissue", "gene~dog+group:tissue")
lymphoma.out <- ms.slgf(dataf=lymphoma.tall, response="gene", lgf="dog", 
                        usermodels=lymphoma.models, prior="zs", 
                        het=c(1,1), min.levels=2)

## -----------------------------------------------------------------------------
lymphoma.out$results[1:5, c(1:3,6)]

## -----------------------------------------------------------------------------
lymphoma.out$class.Probs
head(lymphoma.out$scheme.Probs, 5)

## -----------------------------------------------------------------------------
lymphoma.out$results[1,c(1:3,6,8)]

## -----------------------------------------------------------------------------
lymphoma.out$group.datafs[[18]]


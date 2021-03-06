% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/column_centerer.R
\name{column_centerer}
\alias{column_centerer}
\title{Column centerer for a design matrix.}
\usage{
column_centerer(mm)
}
\arguments{
\item{mm}{a model matrix.}
}
\value{
\code{column_centerer} centers the columns of a design matrix by each column's mean.
}
\description{
\code{column_centerer} Centers the columns of a matrix by the column mean.
}
\examples{
set.seed(314159)
test.data <- data.frame("y"=c(rnorm(10,0,1), rnorm(10,3,1), rnorm(10,5,3)),
                        "x1"=c(rep("A",10), rep("B",10), rep("C",10)),
                        "x2"=rnorm(30,0,1))
m <- lm(y~x1+x2, data=test.data)
mm <- model.matrix(m)
column_centerer(mm)

}

% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/Errors.R
\name{ClusteredErrors}
\alias{ClusteredErrors}
\title{This gives the covariance matrix when clustered errors are assumed.}
\usage{
ClusteredErrors(Model, cluster)
}
\arguments{
\item{Model}{A linear model from the \code{\link{lm}} command.}

\item{cluster}{A vector specifying the clustering variable.}
}
\value{
The covariance matrix for the input linear model with clustered standard errors as specified by the \code{cluster} parameter.
}
\description{
This gives the covariance matrix when clustered errors are assumed.
}
\examples{
x = sin(1:100)
y = 1 + x + c(rnorm(50), rnorm(50, sd = 3))
z = as.factor(c(rep("A", 50), rep("B", 50) ))
dd = data.frame(x = sin(1:100), y = y, z = z)
form = "y ~ x"
rm(x,y,z)
Model  = lm(form = form, data = dd)
CovarianceMatrix = ClusteredErrors(Model, dd$z)
}
\seealso{
\code{\link{lm}}; \code{\link{vcovHC}}
}


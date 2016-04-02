#' Find the covariance matrix for an input linear regression model
#' @export CovMatrix
#' @param Model A linear model from the \code{\link{lm}} command.
#' @param errors The type of standard errors. The options are c("homoskedastic", "robust", "white", "HC0" , "HC2" , "HC3" , "HC4", "HC4m", "HC5") where the last ones are as described by \code{\link{vcovHC}}
#' @param cluster The clustering variable. If input errors (above) are ignored.
#'
#' @return The covariance matrix for the input linear model with the standard errors as specified by the errors parameter.
#' @seealso \code{\link{lm}}; \code{\link{vcovHC}}
#' @examples
#' x = sin(1:100)
#' y = 1 + x + c(rnorm(50), rnorm(50, sd = 3))
#' z = as.factor(c(rep("A", 50), rep("B", 50) ))
#' dd = data.frame(x = sin(1:100), y = y, z = z)
#' form = "y ~ x"
#' rm(x,y,z)

#' M1  = lm(form = form, data = dd)
#' Cov = CovMatrix(M1, cluster = dd$z)

CovMatrix = function(Model, errors = c("homoskedastic", "robust", "white", "HC0" , "HC2" , "HC3" , "HC4", "HC4m", "HC5"), cluster = FALSE){
  if (cluster[[1]] == FALSE){
    SandwichErrorCodes        = c("const"        , "HC1"   ,  "HC"  , "HC0" , "HC2" , "HC3" , "HC4", "HC4m", "HC5")
    names(SandwichErrorCodes) = c("homoskedastic", "robust", "white", "HC0" , "HC2" , "HC3" , "HC4", "HC4m", "HC5")
    CovarianceMatrix = sandwich::vcovHC(Model, type = SandwichErrorCodes[errors][[1]] )
  } else {
    CovarianceMatrix = ClusteredErrors(Model, cluster)
  }
}

#' This gives the covariance matrix when clustered errors are assumed.
#' @export ClusteredErrors
#' @param Model A linear model from the \code{\link{lm}} command.
#' @param cluster A vector specifying the clustering variable.
#'
#' @return The covariance matrix for the input linear model with clustered standard errors as specified by the \code{cluster} parameter.
#' @seealso \code{\link{lm}}; \code{\link{vcovHC}}
#' @examples
#' x = sin(1:100)
#' y = 1 + x + c(rnorm(50), rnorm(50, sd = 3))
#' z = as.factor(c(rep("A", 50), rep("B", 50) ))
#' dd = data.frame(x = sin(1:100), y = y, z = z)
#' form = "y ~ x"
#' rm(x,y,z)

#' Model  = lm(form = form, data = dd)
#' CovarianceMatrix = ClusteredErrors(Model, dd$z)

ClusteredErrors = function(Model, cluster) {
  # http://stats.stackexchange.com/questions/10017/standard-error-clustering-in-r-either-manually-or-in-plm
  M <- length(unique(cluster))
  N <- length(cluster)
  dfc <- (M/(M-1))*((N-1)/(N-Model$rank))
  u <- apply(sandwich::estfun(Model),2, function(x) tapply(x, cluster, sum))
  vcovCL <-dfc * sandwich::sandwich(Model, meat = crossprod(u)/N) * 1
  return(vcovCL)
}

#' This gives regression output for the econreg model. Anova is output along with a coefficient test (with the specified errors), $R^2$, adjusted $R^2$ and residual standard error. Note that the residual standard error is done assuming homoskedastic standard errors. This is called by the \code{\link{reg}} function if \code{quietly = FALSE} or can be called manually.
#' @export summary.econreg
#' @param econreg An econreg object as created by the \code{\link{reg}} function.
#' @param ... Any other options for the stargazer package. See  \code{\link{stargazer}} for descriptions of these.
#'
#' @return This does not return anything. A summary is printed to the console.
#' @seealso \code{\link{anova}} \code{\link{summary.lm}} \code{\link{coeftest}}
#' @examples
#' x = sin(1:100)
#' y = 1 + x + c(rnorm(50), rnorm(50, sd = 3))
#' z = as.factor(c(rep("A", 50), rep("B", 50) ))
#' dd = data.frame(x = sin(1:100), y = y, z = z)
#' form = "y ~ x"
#' rm(x,y,z)

#' M1 = reg(form = form, data = dd, cluster = "z")
#' summary(M1)

#' Probit binary choice model
#' @export probit
#' @param form The regression formula. See \code{\link{lm}} for details on writing this formula.
#' @param data The data.frame for the regression. As subsetting is not currently supported you should drop observations you do not want to include from the data.frame.
#' @param errors The type of standard errors. The options are c("homoskedastic", "robust", "white", "HC0" , "HC2" , "HC3" , "HC4", "HC4m", "HC5") where the last ones are as described by \code{\link{sandwich::vcovHC}}
#' @param cluster The clustering variable. If input errors (above) are ignored.
#' @param quietly Output a summary? FALSE by default.
#' @param warnings Warn about biased and inconsistent parameters in presence of heteroskedasticity?
#'
#' @description Note that robust errors are allowed in this case as is allowed in Stata.
#' It should be noted however that in the presence of heteroskedasticity parameters of a probit model are biased and inconsistent and no standard errors corrections will fix this problem. See \code{\link{http://davegiles.blogspot.com.es/2013/05/robust-standard-errors-for-nonlinear.html}} for more discussion. If there is heteroskedasticity you can consider using the \code{gmlx} package which allows it to be modelled such that it won't bias estimates.
#' The above warning will be displayed by default when robust errors are used without heteroskedasticity being modelled. To disable warnings set warnings = FALSE.
#' The results here will never match Stata's output to all decimal places. The reason is because Stata uses a different numerical optimisation algorithm.
#'
#' @return An econreg object which contains a list:
#'  \enumerate{
#'  \item a linear model
#'  \item the covariance matrix
#'  \item the type of errors
#'  \item Class - The class of the model, probit
#'  }
#' @seealso \code{\link{lm}}; \code{\link{sandwich::vcovHC}}
#' @examples
#' x = sin(1:100)
#' y = round((x - 50)/50 + c(rnorm(50), rnorm(50, sd = 3)))
#' y[y>0] = 1 ; y[y<0] = 0
#' z = as.factor(c(rep("A", 50), rep("B", 50) ))
#' dd = data.frame(x = x, y = y, z = z)
#' form = "y ~ x"
#' rm(x,y,z)

#' M4 = probit(form = form, data = dd)
#' M5 = probit(form = form, data = dd, errors = "robust")

probit = function(form, data, errors = c("homoskedastic", "robust", "white", "HC0" , "HC2" , "HC3" , "HC4", "HC4m", "HC5"), cluster = FALSE, quietly = FALSE, warnings = TRUE){
  errors = errors[1]
  if (errors[1] != "homoskedastic" & warnings){warning("In the presence of heteroskedasticity the coefficients from a probit model are biased and inconsistent. Whilst the input errors have been applied you should check errors are homoskedastic or model the form of heteroskedasticity. See helpfile ?probit for more discussion")}
  Modl = do.call("glm", list(formula = form, family = binomial (link = "probit"), data = substitute(data)))

  if((class(cluster[[1]]) == "character") & length(cluster) == 1){
    cluster = data[,cluster]
    errors = "Clustered"
  }
  CovarianceMatrix = CovMatrix(Modl, errors, cluster)
  if (errors[1] == "homoskedastic") {CovarianceMatrix = summary(Modl)$cov.scaled}

  ModelWithCov = list(Model = Modl, CovarianceMatrix = CovarianceMatrix, errors = errors, Class = "logit" )
  class(ModelWithCov) = "econprobit"

  if (quietly == FALSE){
    summary(ModelWithCov)
  }

  return(ModelWithCov)
}

#' Logit binary choice model
#' @export logit
#' @param form The regression formula. See \code{\link{lm}} for details on writing this formula.
#' @param data The data.frame for the regression. As subsetting is not currently supported you should drop observations you do not want to include from the data.frame.
#' @param errors The type of standard errors. The options are c("homoskedastic", "robust", "white", "HC0" , "HC2" , "HC3" , "HC4", "HC4m", "HC5") where the last ones are as described by \code{\link{sandwich::vcovHC}}
#' @param cluster The clustering variable. If input errors (above) are ignored.
#' @param quietly Output a summary? FALSE by default.
#' @param warnings Warn about biased and inconsistent parameters in presence of heteroskedasticity?
#'
#' @description Note that robust errors are allowed in this case as is allowed in Stata.
#' It should be noted however that in the presence of heteroskedasticity parameters of a logit model are biased and inconsistent and no standard errors corrections will fix this problem. See \code{\link{http://davegiles.blogspot.com.es/2013/05/robust-standard-errors-for-nonlinear.html}} for more discussion. If there is heteroskedasticity you can consider using the \code{gmlx} package which allows it to be modelled such that it won't bias estimates.
#' The above warning will be displayed by default when robust errors are used without heteroskedasticity being modelled. To disable warnings set warnings = FALSE.
#' The results here will never match Stata's output to all decimal places. The reason is because Stata uses a different numerical optimisation algorithm.
#'
#' @return An econreg object which contains a list:
#'  \enumerate{
#'  \item a linear model
#'  \item the covariance matrix
#'  \item the type of errors
#'  \item Class - The class of the model, logit
#'  }
#' @seealso \code{\link{lm}}; \code{\link{sandwich::vcovHC}}
#' @examples
#' x = sin(1:100)
#' y = round((x - 50)/50 + c(rnorm(50), rnorm(50, sd = 3)))
#' y[y>0] = 1 ; y[y<0] = 0
#' z = as.factor(c(rep("A", 50), rep("B", 50) ))
#' dd = data.frame(x = x, y = y, z = z)
#' form = "y ~ x"
#' rm(x,y,z)

#' M4 = logit(form = form, data = dd)
#' M5 = logit(form = form, data = dd, errors = "robust")

logit = function(form, data, errors = c("homoskedastic", "robust", "white", "HC0" , "HC2" , "HC3" , "HC4", "HC4m", "HC5"), cluster = FALSE, quietly = FALSE, warnings = TRUE){
  errors = errors[1]
  if (errors[1] != "homoskedastic" & warnings){warning("In the presence of heteroskedasticity the coefficients from a probit model are biased and inconsistent. Whilst the input errors have been applied you should check errors are homoskedastic or model the form of heteroskedasticity. See helpfile ?probit for more discussion")}
  Modl = do.call("glm", list(formula = form, family = binomial (link = "logit"), data = substitute(data)))

  if((class(cluster[[1]]) == "character") & length(cluster) == 1){
    cluster = data[,cluster]
    errors = "Clustered"
  }
  CovarianceMatrix = CovMatrix(Modl, errors, cluster)
  if (errors[1] == "homoskedastic") {CovarianceMatrix = summary(Modl)$cov.scaled}

  ModelWithCov = list(Model = Modl, CovarianceMatrix = CovarianceMatrix, errors = errors, Class = "logit" )
  class(ModelWithCov) = "econlogit"

  if (quietly == FALSE){
    summary(ModelWithCov)
  }

  return(ModelWithCov)
}



#' This gives regression output for the econprobit model. Anova is output along with a coefficient test (with the specified errors), $R^2$, adjusted $R^2$ and residual standard error. Note that the residual standard error is done assuming homoskedastic standard errors. This is called by the \code{\link{reg}} function if \code{quietly = FALSE} or can be called manually.
#' @export summary.econprobit
#' @param econprobit An econprobit object as created by the \code{\link{probit}} function.
#' @return This does not return anything. A summary is printed to the console.
#' @seealso \code{\link{anova}} \code{\link{summary.glm}} \code{\link{coeftest}}
#' @examples
#'

summary.econprobit = function(econreg){
  cat(paste0("Summary of probit model with ", econreg$errors, " standard errors", "\n"))
  #cat(paste0(econreg$Model["call"], "\n"))
  #print(lmtest::waldtest(econreg$Model, vcov = econreg$CovarianceMatrix))
  #print(anova(econreg$Model))
  #Sm = summary.lm(econreg$Model)
  #cat(paste0("Residual standard error (calculated with homoskedastic errors): ",  format(Sm$sigma, digits = 3) , " on ", format(Sm$df[2], digits = 3) ," degrees of freedom \nMultiple R-squared: ", format(Sm$r.squared, digits = 3), ",	Adjusted R-squared: ", format(Sm$adj.r.squared, digits = 3)))
  print(lmtest::coeftest(x = econreg$Model, vcov. = econreg$CovarianceMatrix))
}

#' This gives regression output for the econlogit model. Anova is output along with a coefficient test (with the specified errors), $R^2$, adjusted $R^2$ and residual standard error. Note that the residual standard error is done assuming homoskedastic standard errors. This is called by the \code{\link{reg}} function if \code{quietly = FALSE} or can be called manually.
#' @export summary.econlogit
#' @param econlogit An econlogit object as created by the \code{\link{logit}} function.
#' @return This does not return anything. A summary is printed to the console.
#' @seealso \code{\link{anova}} \code{\link{summary.glm}} \code{\link{coeftest}}
#' @examples
#'

summary.econlogit = function(econlogit){
  cat(paste0("Summary of logit model with ", econlogit$errors, " standard errors", "\n"))
  #cat(paste0(econreg$Model["call"], "\n"))
  #print(lmtest::waldtest(econreg$Model, vcov = econreg$CovarianceMatrix))
  #print(anova(econreg$Model))
  #Sm = summary.lm(econreg$Model)
  #cat(paste0("Residual standard error (calculated with homoskedastic errors): ",  format(Sm$sigma, digits = 3) , " on ", format(Sm$df[2], digits = 3) ," degrees of freedom \nMultiple R-squared: ", format(Sm$r.squared, digits = 3), ",	Adjusted R-squared: ", format(Sm$adj.r.squared, digits = 3)))
  print(lmtest::coeftest(x = econlogit$Model, vcov. = econlogit$CovarianceMatrix))
}

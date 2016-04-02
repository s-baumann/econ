#' Linear Regression
#' @export reg
#' @param form The regression formula. See \code{\link{lm}} for details on writing this formula.
#' @param data The data.frame for the regression. As subsetting is not currently supported you should drop observations you do not want to include from the data.frame.
#' @param errors The type of standard errors. The options are c("homoskedastic", "robust", "white", "HC0" , "HC2" , "HC3" , "HC4", "HC4m", "HC5") where the last ones are as described by \code{\link{sandwich::vcovHC}}
#' @param cluster The clustering variable. If input errors (above) are ignored.
#' @param quietly Output a summary? FALSE by default.
#'
#' @return An econreg object which contains a list:
#'  \enumerate{
#'  \item{"a linear model"}
#'  \item{"the covariance matrix"}
#'  \item{"the type of errors"}
#'  \item{"Class - The class of the model, lm"}
#'  }
#' @seealso \code{\link{lm}}; \code{\link{sandwich::vcovHC}}
#' @examples
#' x = sin(1:100)
#' y = 1 + x + c(rnorm(50), rnorm(50, sd = 3))
#' z = as.factor(c(rep("A", 50), rep("B", 50) ))
#' dd = data.frame(x = x, y = y, z = z)
#' form = "y ~ x"
#' rm(x,y,z)

#' M1 = reg(form = form, data = dd, cluster = "z")
#' M2 = reg(form = form, data = dd, errors = "homoskedastic")
#' M3 = reg(form = form, data = dd, errors = "robust")

reg = function(form, data, errors = c("homoskedastic", "robust", "white", "HC0" , "HC2" , "HC3" , "HC4", "HC4m", "HC5"), cluster = FALSE, quietly = FALSE){
  errors = errors[1]
  Modl = do.call("lm", list(formula = form, data = substitute(data)))

  if((class(cluster[[1]]) == "character") & length(cluster) == 1){
    cluster = data[,cluster]
    errors = "Clustered"
    }
  CovarianceMatrix = CovMatrix(Modl, errors, cluster)

    ModelWithCov = list(Model = Modl, CovarianceMatrix = CovarianceMatrix, errors = errors, Class = class(Modl))
  class(ModelWithCov) = "econreg"

  if (quietly == FALSE){
    summary(ModelWithCov)
  }

  return(ModelWithCov)
}

summary.econreg = function(econreg){
  cat(paste0("Summary of linear model with ", econreg$errors, " standard errors", "\n"))
  cat(paste0(econreg$Model["call"], "\n"))
  print(lmtest::waldtest(econreg$Model, vcov = econreg$CovarianceMatrix))
  print(anova(econreg$Model))
  Sm = summary.lm(econreg$Model)
  cat(paste0("Residual standard error (calculated with homoskedastic errors): ",  format(Sm$sigma, digits = 3) , " on ", format(Sm$df[2], digits = 3) ," degrees of freedom \nMultiple R-squared: ", format(Sm$r.squared, digits = 3), ",	Adjusted R-squared: ", format(Sm$adj.r.squared, digits = 3)))
  print(lmtest::coeftest(x = econreg$Model, vcov. = econreg$CovarianceMatrix))
}




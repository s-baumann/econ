#' Panel DataLinear Regression
#' @export xtreg
#' @param form The regression formula. See \code{\link{xlm}} in the plm package for details on writing this formula.
#' @param data The data.frame for the regression. As subsetting is not currently supported you should drop observations you do not want to include from the data.frame.
#' @param errors The type of standard errors. The options are c("homoskedastic", "robust", "white", "HC0" , "HC2" , "HC3" , "HC4", "HC4m", "HC5") where the last ones are as described by \code{\link{sandwich::vcovHC}}
#' @param cluster The clustering variable. If input errors (above) are ignored.
#' @param idtime The group id variable name and time variable entered as a vector. ie. c("individual", "time"). If omitted first two variables in data.frame are used for this.
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
#' @description If you get a names attribute error this can be caused by trying to do a within transformation on data where there is no withingroup variation in y.
#' @examples
#' x = sin(1:100)
#' y = 1 + x + c(rnorm(50), rnorm(50, sd = 3))
#' z = as.factor(c(rep("A", 50), rep("B", 50) ))
#' t = c(1:50, 1:50)
#' dd = data.frame(z = z, t = t, x = x, y = y)
#' form = "y ~ x"
#' rm(x,y,z)

#' M1 = reg(form = form, data = dd, cluster = "z")
#' M2 = reg(form = form, data = dd, errors = "homoskedastic")
#' M3 = reg(form = form, data = dd, errors = "robust")

xtreg = function(form, data, errors = c("homoskedastic", "robust", "white", "HC0" , "HC2" , "HC3" , "HC4", "HC4m", "HC5"), cluster = FALSE, quietly = FALSE, idtime = NULL){
  errors = errors[1]
  if (is.null(idtime)){data = plm::pdata.frame(data)} else {data = plm::pdata.frame(data, index = idtime)}
  Modl = do.call(plm::plm, list(formula = form, data = substitute(data)))

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

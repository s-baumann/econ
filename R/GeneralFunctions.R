#' This is a function to estimate the tests of model significance with nonstandard standard errors.
#' @export SignificanceTests
#' @param Model The model.
#' @param CovarianceMatrix The Covariance Matrix
#' @param Class The class of model. Could be "lm", "probit".
#' @param digits How many digits to output the f statistic.
#'
#' @return Rreturns an F Statistic with stars for significance.
#' @seealso \code{\link{stargazer}}
#' @examples
#' x = sin(1:100)
#' y = 1 + x + c(rnorm(50), rnorm(50, sd = 3))
#' z = as.factor(c(rep("A", 50), rep("B", 50) ))
#' dd = data.frame(x = sin(1:100), y = y, z = z)
#' form = "y ~ x"
#' rm(x,y,z)

#' M1 = reg(form = form, data = dd, cluster = "z")
#' SignificanceTests(M1$Model, M1$CovarianceMatrix, M1$Class, digits = 4)

SignificanceTests = function(Model, CovarianceMatrix, Class, digits = 4){
  if (Class == "lm"){
    FFF = lmtest::waldtest(Model, vcov = CovarianceMatrix)
    FF  = format(FFF$F[2], digits = digits)
    ProbF = FFF$`Pr(>F)`[2]
    FF[0.05 < ProbF & ProbF < 0.1]   = paste0(FF[0.05 <ProbF & ProbF < 0.1], "*")
    FF[0.01 < ProbF & ProbF <= 0.05] = paste0(FF[0.01 <ProbF & ProbF <= 0.05], "**")
    FF[ProbF <= 0.01] = paste0(FF[ProbF <= 0.01], "***")
    return(FF)
  }
  if (Class == "probit" | Class == "logit"){
    FFF =  aod::wald.test(b=coef(Model), Sigma=CovarianceMatrix, Terms = which(names(coef(Model)) != "(Intercept)"))
    FFF = FFF$result$chi2
    FF  = format(FFF["chi2"][[1]], digits = digits)
    ProbF = FFF["P"][[1]]
    FF[0.05 < ProbF & ProbF < 0.1]   = paste0(FF[0.05 <ProbF & ProbF < 0.1], "*")
    FF[0.01 < ProbF & ProbF <= 0.05] = paste0(FF[0.01 <ProbF & ProbF <= 0.05], "**")
    FF[ProbF <= 0.01] = paste0(FF[ProbF <= 0.01], "***")
    FF = NA
  }
}

#' This is a wrapper function for the stargazer package so that econreg objects can easily be passed to it.
#' @export gazer
#' @param models A econreg object or a list of econreg objects.
#' @param ... Any other options for the stargazer package. See  \code{\link{stargazer}} for descriptions of these.
#'
#' @return This does not return anything. The table is either printed to the console or saved on the filing system depending on what options are input.
#' @seealso \code{\link{stargazer}}
#' @examples
#' x = sin(1:100)
#' y = 1 + x + c(rnorm(50), rnorm(50, sd = 3))
#' z = as.factor(c(rep("A", 50), rep("B", 50) ))
#' dd = data.frame(x = sin(1:100), y = y, z = z)
#' form = "y ~ x"
#' rm(x,y,z)

#' M1 = reg(form = form, data = dd, cluster = "z")
#' gazer(M1)
#' M2 = reg(form = form, data = dd, errors = "homoskedastic")

#' models = list(M1, M2)
#' gazer(models)

gazer = function(models,  type = "text", digits = 7,...){
  if (is.null(models$Model)){
    # A list of econregs is input
    mods = sapply(models, `[`, "Model")
    CovarianceMatrices = sapply(models, `[`, "CovarianceMatrix")
    Classes = sapply(models, `[`, "Class")
    errors = lapply(CovarianceMatrices , function(x) sqrt(diag(x)))
    FF = unlist(sapply(1:length(models), function(x) SignificanceTests(mods[[x]], CovarianceMatrices[[x]], Classes[[x]], digits = digits)))
  } else {
    # A single econreg is input
    mods = models$Model
    CovarianceMatrices = models$CovarianceMatrix
    Classes = models$Class
    errors = sqrt(diag(models$CovarianceMatrix))
    FF = SignificanceTests(mods, CovarianceMatrices, Classes , digits = digits)
  }
  stargazer::stargazer(mods, se = errors, type = type, omit.stat  = c("f", "aic","ll") , add.lines = list(c("F Statistic", FF)) , digits = digits, ... )
}

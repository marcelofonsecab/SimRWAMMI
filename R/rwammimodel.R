#' @export rwammi.model

rwammi.model <- function(dataframe, weight = NULL, Ncomp = 2){
  data = transform_usable_data(dataframe, median, "dataframe")
  Ngen = nlevels(data$gen)
  Nenv = nlevels(data$env)
  weight = c(weight)

  model = rlm(yield ~ gen + env, data = data, w = weight)
  residual.matrix = matrix(model$residuals, ncol = Nenv, nrow = Ngen)
  fitted.values = matrix(model$fitted.values, ncol = Nenv, nrow = Ngen)

  SVD.aux = WeightedSvdRes(residual.matrix, weight = weight, Ncomp = Ncomp)
  colnames(SVD.aux) = colnames(residual.matrix)
  rownames(SVD.aux) = rownames(residual.matrix)
  SVD.red = SVD.aux

  return(list(residual = residual.matrix, matrix.fitted = fitted.values,
              reduced.SVD = SVD.red))
}


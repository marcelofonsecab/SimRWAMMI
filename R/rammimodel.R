#' @export rammi.model

rammi.model <- function(dataframe, weight = NULL, Ncomp = 2){
  data = transform_usable_data(dataframe, median, "dataframe")
  Ngen = nlevels(data$gen)
  Nenv = nlevels(data$env)
  weight = c(weight)

  if(is.null(weight)){
    model = rlm(yield ~ gen + env, data = data)
  } else {
    model = rlm(yield ~ gen + env, data = data, w = weight)
  }
  residual.matrix = matrix(model$residuals, ncol = Nenv, nrow = Ngen)
  fitted.values = matrix(model$fitted.values, ncol = Nenv, nrow = Ngen)

  SVD = rSVDdpd(residual.matrix, alpha = 0.3, maxiter = 100,
                initu= svd(residual.matrix)$u, initv = svd(residual.matrix)$v,
                eps = 1e-4)
  SVD.red = SVD$u[,1:Ncomp]%*%diag(SVD$d[1:Ncomp])%*%t(SVD$v[,1:Ncomp])

  return(list(residual = residual.matrix, matrix.fitted = fitted.values,
              reduced.SVD = SVD.red, SVD = SVD))
}


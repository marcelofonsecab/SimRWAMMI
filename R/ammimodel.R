#' @export ammi.model

ammi.model <- function(dataframe, Ncomp = 2){
  df.mean = transform_usable_data(dataframe, mean, type = "dataframe")
  Ngen = nlevels(df.mean$gen)
  Nenv = nlevels(df.mean$env)
  env.names = levels(df.mean$env)
  gen.names = levels(df.mean$gen)

  #Modelo simples AMMI
  model.GEI = lm(yield ~ gen + env, data = df.mean)
  matrix.GEI = matrix(residuals(model.GEI), ncol = Nenv, nrow = Ngen)
  svd.GEI = svd(matrix.GEI)
  svd.u = svd.GEI$u
  svd.d = diag(svd.GEI$d)
  svd.v = svd.GEI$v

  #Colocando a parte multiplicativa na matriz dependendo de quantos NComp foi colocado na fun??o
  #Caso coloque Ncomp = 0, a matriz ser? apenas a soma dos efeitos aditivos
  aux = matrix(model.GEI$fitted.values, ncol = Nenv)
  if(Ncomp >= 1){
    for(i in 1:Ncomp){
      aux = aux + (svd.u[,i]*svd.d[i,i])%*%t(svd.v[,i])
    }
  }
  colnames(aux) = env.names
  rownames(aux) = gen.names
  matrix.fitted = matrix(model.GEI$fitted.values, ncol = Nenv)

  SVD.red = 0
  if(Ncomp == 1){
    SVD.red = 0
  } else {
    SVD.red = svd.u[, 1:Ncomp]%*%(svd.d[1:Ncomp, 1:Ncomp])%*%t(svd.v[, 1:Ncomp])
  }
  colnames(SVD.red) = env.names
  rownames(SVD.red) = gen.names

  return(list(estimated.data = aux, anovatable = anova(model.GEI), SVD = svd.GEI,
              residual = matrix.GEI, matrix.fitted = matrix.fitted, reduced.SVD = SVD.red))
}


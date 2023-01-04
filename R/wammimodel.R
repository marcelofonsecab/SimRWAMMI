#' @export wammi.model

wammi.model <- function(data, weight = NULL, Ncomp = 2){
  dataframe = transform_usable_data(data, mean, type = "dataframe")
  Nenv = nlevels(dataframe$env)
  if(is.null(weight)){
    #Se n?o colocar o vetor de peso a fun??o ir? parar
    stop("Essa fun??o precisa de um vetor de pesos")
  }

  weight = c(weight) #Vetor de pesos
  weighted.lm = lm(yield ~ gen + env, weights = weight, data = dataframe)

  #valores preditos na regress?o ponderada
  fitted.values = matrix(weighted.lm$fitted.values, ncol = Nenv)

  #Modelo linear com pesos

  residual.matrix = matrix(weighted.lm$residuals, ncol = Nenv)
  #Matriz de residuos do modelo ponderados
  colnames(residual.matrix) = levels(dataframe$env)
  rownames(residual.matrix) = levels(dataframe$gen)

  #W-SVD
  SVD.aux = WeightedSvdRes(residual.matrix, weight = weight, Ncomp = Ncomp)
  colnames(SVD.aux) = colnames(residual.matrix)
  rownames(SVD.aux) = rownames(residual.matrix)
  SVD.red = SVD.aux

  dataframe = cbind(dataframe, "W.residuals" = melt(residual.matrix)$value)
  return(list(dataframe = dataframe,
              residual = residual.matrix, matrix.fitted = fitted.values, reduced.SVD = SVD.red))
}


#' @export Data_Contamination

Data_Contamination = function(data, porcentagem = 5, seed = 1,
                              tipo = "shift", k = 5, c = 10){
  set.seed(seed)
  tipos_aux = c("shift", "pointmass", "highvariance")
  match.arg(tipo, tipos_aux)
  percentage_aux = ifelse(porcentagem > 1, porcentagem/100, porcentagem)
  df_aux = data[data$rep==1, ]
  matrix_aux = tapply(data$yield,
                      data[,c("gen","env")],
                      mean)
  Medias_Desvios = cbind(media = apply(matrix_aux, 2, mean),
                         desvio = apply(matrix_aux, 2, sd))

  Nenv = levels(df_aux$env) %>% length()
  Ngen = levels(df_aux$gen) %>% length()
  Nrow = nrow(df_aux)

  sample_aux = sample(1:Nrow, Nrow * percentage_aux, replace = F)
  df_tmp = df_aux[sample_aux,]

  for(i in 1:Nenv){
    sample_tmp = which(df_tmp$env == rownames(Medias_Desvios)[i])
    df_tmp[sample_tmp, "yield"] = rnorm(length(sample_tmp),
                                        mean = Medias_Desvios[i,1] + ifelse(tipo == "shift",
                                                                            k,
                                                                            0) * Medias_Desvios[i,2],
                                        sd = Medias_Desvios[i,2] / ifelse(tipo == "pointmass" |
                                                                            tipo == "highvariance",
                                                                          sqrt(c),
                                                                          1))
  }

  df_aux[sample_aux,] = df_tmp
  df_aux = rbind(df_aux, data[data$rep != 1,])
  data.full = list(data.sem.cont = data, data.cont = df_aux)

  return(data.full)
}

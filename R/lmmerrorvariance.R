#' @export LMM.Error_variance

LMM.Error_variance = function(data){
  quant.gen = length(unique(data$gen))
  quant.env = length(unique(data$env))
  if(quant.env == 1){
    aux = lmer(yield ~ rep + (1|gen), data = data)
    tmp = attr(VarCorr(aux), "sc")^2
  } else if(quant.gen == 1){
    aux = lmer(yield ~ rep + (1|env), data = data)
    tmp = attr(VarCorr(aux), "sc")^2
  }
  return(tmp)
}

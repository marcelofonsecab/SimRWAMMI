#' @export Create.Replications

Create.Replications = function(data, reps = 2, sig = 1){

  Nenv = length(levels(data$env))
  Ngen = length(levels(data$gen))
  N = Nenv * Ngen
  data.tmp = data.aux = NULL

  for(i in 1:reps){
    resids = rnorm(N, mean = 0, sd = sqrt(sig))
    data.aux = data
    data.aux$yield = data.aux$yield + resids
    data.aux$rep = i
    data.tmp = rbind(data.tmp, data.aux)
  }


  return(data.tmp)
}

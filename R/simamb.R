#' @export sim.amb

sim.amb <- function(seed = NULL, Ngen = 100, Nenv = 8, Ncomp = 2,
                    effectGlobal = c(mean = 15, sd = sqrt(3)),
                    effectGen = c(mean = 5, sd = 1),
                    effectEnv = c(mean = 8, sd = sqrt(2)),
                    k = rep(1, Nenv)){
  #seed = vetor com valor para a semente para a simula??o
  #Ngen = n?mero de gen?tipos
  #Nenv = n?mero de ambientes
  #Nrep = n?mero de replica??es #Fixado igual a 1
  #Ncomp = n?mero de componentes principais utilizadas
  #effectGlobal = m?dia global da simula??o (vetor contendo a m?dia e o desvio padr?o)
  #effectGen = efeito do gen?tipo (vetor contendo a m?dia e o desvio padr?o)
  #effectEnv = efeito do ambiente (vetor contendo a m?dia e o desvio padr?o)
  seed.aux = ifelse(is.null(seed), sample(1:2^16, 1), seed) #Caso n?o tenha colocado uma semente, gerar? uma semente aleat?ria
  set.seed(seed.aux)

  globalMean = rnorm(1, effectGlobal[1], effectGlobal[2]) #M?dia global
  alpha = rnorm(Ngen, mean = effectGen[1], effectGen[2]) #Efeito do gen?tipo
  beta = rnorm(Nenv, mean = effectEnv[1], effectEnv[2]) #Efeito do ambiente

  rand.mat = matrix(runif(Ngen*Nenv, min = -0.5, max =  0.5), ncol = Nenv)
  rand.svd = svd(rand.mat) #SVD da matriz de valores aleatorios
  rand.u = rand.svd$u #Matriz U
  rand.v = rand.svd$v #Matriz V
  rand.d = diag(rand.svd$d) #Vetor diagonal da matriz D

  simulated.amb = matrix(rep(1, Ngen))%*%t(matrix(rep(1, Nenv)))*globalMean +
    alpha%*%t(matrix(rep(1, Nenv))) +
    matrix(rep(1, Ngen))%*%t(beta) #Matriz apenas com os efeitos
  #Simulando o ambiente com N componentes principais a partir da matriz al?atoria via U(-0.5, 0.5)
  for(j in 1:Ncomp){
    simulated.amb = simulated.amb +
      ((rand.u[,j]*rand.d[j,j])%*%t(rand.v[,j]))*(k[j]) #Adi??o dos efeitos aleatorios
  }

  #Organizando a matriz e transformando em dataframe
  aux = matrix(rep(0, Ngen*Nenv*4), ncol = 4)
  aux = as.data.frame(aux)
  colnames(aux) = c("gen", "env", "rep", "yield")
  aux$gen = rep(paste0("G", sprintf('%0.3d', 1:Ngen)), Nenv)
  aux[order(aux$gen),]$env = rep(paste("E", 1:Nenv, sep=""), Ngen)
  aux$yield = c(simulated.amb)

  #Separando em uma dataframe com todas as replica??es e um dataframe das medias com base nas rep
  finaldf = aux
  finaldf$gen = as.factor(finaldf$gen)
  finaldf$env = as.factor(finaldf$env)
  finaldf$rep = 1
  finaldf$rep = as.factor(finaldf$rep)
  data = as.data.frame(finaldf)

  return(data)
}

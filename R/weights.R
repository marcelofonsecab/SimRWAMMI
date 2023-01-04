Weights = function(data, errorvariances){

  Ngen = nlevels(data$gen)
  Nenv = nlevels(data$env)
  Nrep = max(data$rep)
  rep.weights = tapply(data[,"rep"],
                       data[,c("gen", "env")],
                       max) / Nrep
  lmm.env.weights.tmp = 1/errorvariances$EnvironmentLMM$Error_Variance
  lmm.Weights.env.tmp = matrix(rep(lmm.env.weights.tmp/max(lmm.env.weights.tmp), Ngen),
                               ncol = Nenv, byrow = T)
  lmm.Weights.env = lmm.Weights.env.tmp * rep.weights

  lmm.gen.weight.tmp = 1/errorvariances$GenotypeLMM$Error_Variance
  lmm.Weights.gen.tmp = matrix(rep(lmm.gen.weight.tmp/max(lmm.gen.weight.tmp), Nenv),
                               nrow = Ngen, byrow = F)
  lmm.Weights.gen = lmm.Weights.gen.tmp * rep.weights

  data.rlm = transform_usable_data(data, median, "dataframe")
  model.rlm.tmp = rlm(yield ~ gen + env, data = data.rlm)

  rlm.weights.tmp = model.rlm.tmp$w
  Weights.rlm = matrix(rlm.weights.tmp, ncol = Nenv, byrow = F)

  # Combinations LMM

  lmm.Weights.gendotenv = lmm.Weights.gen.tmp * lmm.Weights.env.tmp * rep.weights
  lmm.Weights.genplusenv = ( (lmm.Weights.gen.tmp + lmm.Weights.env.tmp) / 2 ) * rep.weights
  lmm.Weights.rlmdotenv = Weights.rlm * lmm.Weights.env.tmp * rep.weights
  lmm.Weights.rlmdotgen = Weights.rlm * lmm.Weights.gen.tmp * rep.weights
  lmm.Weights.rlmdotgendotenv = Weights.rlm * lmm.Weights.gen.tmp * lmm.Weights.env.tmp * rep.weights
  lmm.Weights.rlmdotgenplusenv = Weights.rlm * ( (lmm.Weights.gen.tmp + lmm.Weights.env.tmp) / 2 ) * rep.weights

  # Combinations RLMM

  rlmm.env.weights.tmp = 1/errorvariances$EnvironmentRLMM$Error_Variance
  rlmm.Weights.env.tmp = matrix(rep(rlmm.env.weights.tmp/max(rlmm.env.weights.tmp), Ngen),
                                ncol = Nenv, byrow = T)
  rlmm.Weights.env = rlmm.Weights.env.tmp * rep.weights

  rlmm.gen.weight.tmp = 1/errorvariances$GenotypeRLMM$Error_Variance
  rlmm.Weights.gen.tmp = matrix(rep(rlmm.gen.weight.tmp/max(rlmm.gen.weight.tmp), Nenv),
                                nrow = Ngen, byrow = F)
  rlmm.Weights.gen = rlmm.Weights.gen.tmp * rep.weights

  rlmm.Weights.gendotenv = rlmm.Weights.gen.tmp * rlmm.Weights.env.tmp * rep.weights
  rlmm.Weights.genplusenv = ( (rlmm.Weights.gen.tmp + rlmm.Weights.env.tmp) / 2 ) * rep.weights
  rlmm.Weights.rlmdotenv = Weights.rlm * rlmm.Weights.env.tmp * rep.weights
  rlmm.Weights.rlmdotgen = Weights.rlm * rlmm.Weights.gen.tmp * rep.weights
  rlmm.Weights.rlmdotgendotenv = Weights.rlm * rlmm.Weights.gen.tmp * rlmm.Weights.env.tmp * rep.weights
  rlmm.Weights.rlmdotgenplusenv = Weights.rlm * ( (rlmm.Weights.gen.tmp + rlmm.Weights.env.tmp) / 2 ) * rep.weights


  return(list(LMM = list(lmm.Weights.Env = lmm.Weights.env,
                         lmm.Weights.Gen = lmm.Weights.gen,
                         lmm.Weights.GendotEnv = lmm.Weights.gendotenv,
                         lmm.Weights.GenplusEnv = lmm.Weights.genplusenv,
                         Weights.Rlm = Weights.rlm,
                         lmm.Weights.RlmdotEnv = lmm.Weights.rlmdotenv,
                         lmm.Weights.RlmdotGen = lmm.Weights.rlmdotgen,
                         lmm.Weights.RlmdotGendotEnv = lmm.Weights.rlmdotgendotenv,
                         lmm.Weights.RlmdotGenplusEnv = lmm.Weights.rlmdotgenplusenv),
              RLMM = list(rlmm.Weights.Env = rlmm.Weights.env,
                          rlmm.Weights.Gen = rlmm.Weights.gen,
                          rlmm.Weights.GendotEnv = rlmm.Weights.gendotenv,
                          rlmm.Weights.GenplusEnv = rlmm.Weights.genplusenv,
                          Weights.Rlm = Weights.rlm,
                          rlmm.Weights.RlmdotEnv = rlmm.Weights.rlmdotenv,
                          rlmm.Weights.RlmdotGen = rlmm.Weights.rlmdotgen,
                          rlmm.Weights.RlmdotGendotEnv = rlmm.Weights.rlmdotgendotenv,
                          rlmm.Weights.RlmdotGenplusEnv = rlmm.Weights.rlmdotgenplusenv)))
}


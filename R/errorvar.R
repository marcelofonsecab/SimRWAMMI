Error_Var = function(data, cluster = NULL){

  if(!is.null(cluster)){
    print("Usando loop em paralelo")
  }

  tmp1 = data %>%
    group_by(gen) %>%
    group_split() %>%
    pblapply(LMM.Error_variance, cl = cluster) %>%
    unlist() %>%
    tibble(Group = unique(data$gen), Error_Variance = .)

  tmp2 = data %>%
    group_by(env) %>%
    group_split() %>%
    pblapply(LMM.Error_variance, cl = cluster) %>%
    unlist() %>%
    tibble(Group = unique(data$env), Error_Variance = .)

  tmp3 = data %>%
    group_by(gen) %>%
    group_split() %>%
    pblapply(RLMM.Error_variance, cl = cluster) %>%
    unlist() %>%
    tibble(Group = unique(data$gen), Error_Variance = .)

  tmp4 = data %>%
    group_by(env) %>%
    group_split() %>%
    pblapply(RLMM.Error_variance, cl = cluster) %>%
    unlist() %>%
    tibble(Group = unique(data$env), Error_Variance = .)

  tmp = list(GenotypeLMM = tmp1, EnvironmentLMM = tmp2,
             GenotypeRLMM = tmp3, EnvironmentRLMM = tmp4)

  return(tmp)
}

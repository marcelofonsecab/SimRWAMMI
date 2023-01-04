#' @export transform_usable_data

transform_usable_data = function(dataframe, func, type = c("dataframe", "matrix")){
  type_aux = c("dataframe", "matrix")
  match.arg(type, type_aux)

  if(type == "dataframe"){
    df = dataframe %>%
      group_by(gen, env) %>%
      summarise(yield = func(yield)) %>%
      arrange(env)
  } else if(type == "matrix"){
    df = tapply(dataframe[, "yield"], dataframe[, c("gen", "env")], func)
  }
  return(df)
}

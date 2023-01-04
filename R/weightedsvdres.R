#' @export WeightedSvdRes

WeightedSvdRes <- function(data, weight = NULL, Ncomp = 2){
  Ngen = nrow(data)
  Nenv = ncol(data)

  W = c(weight)
  D = data
  X = matrix(0, Ngen, Nenv)

  aux = matrix(1, Ngen, Nenv)
  Xold = Inf*aux
  Err = Inf
  eps = 1e-10
  Xold = X
  wsvd = svd(W*D + (1-W)*X)
  U = wsvd$u
  S = diag(wsvd$d)
  V = wsvd$v
  if(Ncomp + 1 < length(wsvd$d)){
    S[(Ncomp + 1):length(wsvd$d), (Ncomp + 1):length(wsvd$d)] = 0
  }
  X = U%*%S%*%t(V)
  while(Err > eps){
    Xold = X
    wsvd = svd(W*D + (1-W)*X)
    U = wsvd$u
    S = diag(wsvd$d)
    V = wsvd$v
    if(Ncomp + 1 < length(wsvd$d)){
      S[(Ncomp + 1):length(wsvd$d), (Ncomp + 1):length(wsvd$d)] = 0
    }
    X = U%*%S%*%t(V)
    Err = sum(sum((X-Xold)^2))
  }
  return(X)
}

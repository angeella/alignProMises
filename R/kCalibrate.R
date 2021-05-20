#' @title von Mises Fisher Procrustes model
#' @description perform the functional alignment using the von Mises Fisher Procrustes model
#' @usage kCalibrate(data)
#' @param D squared euclidean distance matrix
#' @param p distances to compute bandwith value
#' @author Angela Andreella and Daniela Corbetta
#' @return Returns list of matrices
#' @export
#' @importFrom diffusionMap epsilonCompute 
#' 
#' 

kCalibrate <- function(D, p){
  
  eps <- epsilonCompute(D, p = p)
  #kQ <- exp(D*(1-1/eps)) %*% solve(exp(-D))
  kQ <- exp(-D/eps) 
  return(kQ)  
}
#' @title von Mises Fisher Procrustes model
#' @description perform the functional alignment using the von Mises Fisher Procrustes model
#' @usage kCalibrate(data)
#' @param D squared euclidean distance matrix
#' @param p distances to compute bandwith value
#' @author Angela Andreella
#' @return Returns list of matrices
#' @export
#' @import diffusionMap epsilonCompute 
#' 
#' 

kCalibrate <- function(D, p){
  
  eps <- epsilonCompute(D, p = p)
  kQ <- exp(D*(1-1/eps))
    
  return(kQ)  
}
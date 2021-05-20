#' @title Calibrate parameters 
#' @description Calibrates the parameters of the von Mises-Fisher distribution to be used as prior in the von Mises-Fisher-Procrustes model.
#' @usage kCalibrate(D, p = 0.01)
#' @param D squared euclidean distance matrix between the coordinates of the voxels
#' @param p distances to compute bandwith value. Default value is 0.01
#' @author Angela Andreella and Daniela Corbetta
#' @return Returns directly \code{k*Q}, i.e. the product between the concentration parameter and the location matrix parameter of the prior distribution
#' @export
#' @importFrom diffusionMap epsilonCompute 
#' 
#' 

kCalibrate <- function(D, p = 0.01){
  
  eps <- epsilonCompute(D, p = p)
  #kQ <- exp(D*(1-1/eps)) %*% solve(exp(-D))
  kQ <- exp(-D/eps) 
  return(kQ)  
}
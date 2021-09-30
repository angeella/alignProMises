#' @title ProMises model
#' @description Perform functional alignment of a matrix by von Mises Fisher Procrustes model with known reference matrix
#' @usage GPASub(X, Q = NULL, k, kQ = NULL, ref_ds, 
#' scaling = TRUE, reflection = TRUE, centered = TRUE)
#' @param X data, i.e., matrix with dimension time points - voxels 
#' @param Q value of the location parameter of the prior distribution. It has dimension voxels x voxels, it could be not symmetric.
#' @param k value of the concentration parameter of the prior distribution
#' @param kQ product between the location parameter and the concentration parameter of the prior distribution. Could be supplied instead of \code{k} and \code{Q}
#' @param ref_ds reference matrix to align
#' @param scaling Flag to apply scaling transformation
#' @param reflection Flag to apply reflection transformation
#' @param centered centered data?
#' @author Angela Andreella and Daniela Corbetta
#' @return \code{vMFP} returns a list with two components:
#' \item{\code{Xest}}{the aligned matrix}
#' \item{\code{R}}{the rotation matrix}
#' @export 
#' @importFrom rARPACK svds
#' @useDynLib ProMises
#' @importFrom Rcpp sourceCpp




GPASub <- function(X, Q = NULL, k, kQ = NULL, ref_ds, scaling = TRUE, reflection = TRUE, centered = TRUE){
  
  nc <- min(dim(X)[1:2])
  #nc <- dim(X)[1]
  
  if(is.null(Q)){ Q <- matrix(0, nrow = nc, ncol = nc)
  }
  if(is.null(kQ)){
    out <- svdC(t(t(ref_ds) %*% X + k * t(Q))) 
  }
  else{
    #Put transposes to save memory.
    out <- svdC(t(t(ref_ds) %*% X + t(kQ)))
  }
  s <- out$S
  U <- out$U
  Vt <- t(out$V)
  if(!reflection){
    s_new <- diag(length(s))
    # s_new <- Diagonal(n = length(s), x = 1)
    s_new[nc,nc] <- sign(det(U %*% Vt))
    Tr <- (U %*% s_new) %*% Vt
    scale <- sum(s_new %*% s)
  }else{
    Tr <-  U %*% Vt
    scale <-  sum(s)
  }
  R = Tr
  if(!centered){
    scale <- scale / norm(X %*% R, type = "F")
  }
  if(!scaling){
    Xest <- X %*% R
  }else{
    Xest <- X %*% R * scale
  }
  
  return(list(Xest = Xest, R = R))
}


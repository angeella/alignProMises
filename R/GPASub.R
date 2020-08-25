#' @title Generalized Procrustes Analysis with prior
#' @description perform Generalized Procrustes Analysis with prior
#' @usage GPASub(X, Q = NULL, k, ref_ds, scaling = TRUE, reflection = TRUE)
#' @param X data, i.e., list of matrices with dimension time points - voxels 
#' @param Q maximum number of iteration
#' @param k value of the concentration parameter of the prior distribution
#' @param ref_ds starting matrix to align
#' @param scaling Flag to apply scaling transformation
#' @param reflection Flag to apply reflection transformation
#' @author Angela Andreella
#' @return Returns list of matrices
#' @export 
#' @import svds rARPACK
#' @import Diagonal Matrix
#' 
GPASub <- function(X, Q = NULL, k, ref_ds, scaling = TRUE, reflection = TRUE){
  
  nc <- dim(X)[1]
  
  if(is.null(Q)){ Q <- matrix(0, nrow = nc, ncol = nc)
  }
  #Put transposes to save memory.
  out <- svdC(t(t(ref_ds) %*% X + k * t(Q))) 
  s <- out$d
  U <- out$u
  Vt <- t(out$v)
  if(!reflection){
    s_new <- Diagonal(n=length(s), x=1)
    s_new[nc,nc] <- sign(det(U %*% Vt))
    Tr <- (U %*% s_new) %*% Vt
    scale <- sum(s_new * s)
  }else{
    Tr <-  U %*% Vt
    scale <-  sum(s)
  }
  R = Tr
  if(!scaling){
      Xest <- X %*% t(R)
  }else{
      Xest <- X %*% t(R) * scale
  }
  
  return(list(Xest = Xest, R = R))
}
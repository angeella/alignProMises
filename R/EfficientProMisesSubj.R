#' @title Efficient ProMises model
#' @description Performs the functional alignment using the Efficient ProMises model allowing for different number of columns (voxel or pixel) between matrices.
#' @usage EfficientProMisesSubj(data, maxIt=10, t =.001, k = 0, Q = NULL, 
#' ref_ds = NULL, scaling = T, reflection= T, subj= T, coord = NULL, 
#' singleQ = F, l = NULL)
#' @param data data, i.e., list of matrices with dimension time points - voxels. Matrices can have different number of columns but must have the same numbers of rows.
#' @param maxIt maximum number of iteration
#' @param t the threshold value to be reached as the minimum relative reduction between the mean matrices
#' @param k value of the concentration parameter of the prior distribution
#' @param Q value of the location parameter of the prior distribution. It has dimension time-points x time-points, it could be not symmetric.
#' @param ref_ds starting matrix to align. If \code{NULL}, it is set to the element-wise arithmetic mean of the reduced matrices
#' @param scaling Flag to apply scaling transformation
#' @param reflection Flag to apply reflection transformation
#' @param subj Flag if each subject has his/her own set of voxel after voxel selection step
#' @param coord list with 3-dim or 2-dim coordinates of the variables. Matrices can have different coordinates. If the location parameter \code{Q = NULL}, then \code{coord} is used to compute it
#' @param l number of singular vectors required to obtained the reduced transformation of the matrices. If \code{NULL}, it is set equal to the number of rows. 
#' @param singleQ if \code{TRUE}, the same location parameter is assumed for the prior distribution of all the rotation parameters.
#' @author Angela Andreella and Daniela Corbetta
#' @return \code{EfficientProMisesSubj} returns a list with five components:
#' \item{\code{Xest}}{a list with the aligned matrices}
#' \item{\code{R}}{an array with the rotation matrices}
#' \item{\code{dist}}{a vector with length equal to the number of iterations that contains the distances between a reference matrix and the previous one}
#' \item{\code{count}}{the number of iterations done by the algorithm}
#' \item{\code{M}}{the element-wise mean matrix of the aligned matrices.}
#' @references For the theory on the Efficient ProMises model see: A. Andreella and L. Finos
#' (2022), Procrustes analysis for high-dimensional data, 
#' Psychometrika 87, 1422-1438
#' @export
#' @importFrom plyr aaply 
#' @importFrom foreach foreach
#' @importFrom foreach %dopar%
#' @importFrom rARPACK svds

EfficientProMisesSubj <- function(data, maxIt=10, t =.001, k = 0, Q = NULL, ref_ds = NULL, scaling = T, reflection= T, subj= T, coord = NULL, l = NULL, singleQ = F){
  
  if(!is.list(data)){warnings("Please insert a list of matrices with dimension time points - voxels")}
  
  row <- sapply(data, nrow)[1]
  col <- sapply(data, ncol)
  nsubj <- length(data)
  
  count = 1
  dist = vector()
  dist[1] <- Inf
  
  j <- row
  
  if(is.null(l))
    l <- row
  
  V <- foreach(i = c(1:nsubj)) %dopar% {
    out <- svds(data[[i]], k = l)
    out$v 
    
  }
 
  Xstar <- array(NA, dim=c(j,l,nsubj))

  for (i in 1:nsubj){
    Xstar[,,i] <- data[[i]] %*% V[[i]] 
  }
  
  ref_ds <- colMeans(aperm(Xstar, c(3, 1, 2)))
  
  Xest <-  array(NA, dim(Xstar))
  R <-  array(NA, c(l, l, nsubj))
  
  # Compute Q*=exp(-D)
    if(is.null(Q)){
      if(subj){
        Qstar <- foreach(i = c(1:nsubj)) %dopar% {
          if(!is.null(coord)){
            
            coord_star <- t(V[[i]]) %*% coord[[i]]
            D <- dist(coord_star, method = "euclidean", diag = T, upper = T)
            exp(-as.matrix(D))
          }
          else 
            matrix(0, nrow = l, ncol = l)
        }
        
      }
      
      else{
        Qstar <- foreach(i = c(1:nsubj)) %dopar% {
          if(!is.null(coord)){
            
            coord_star <- t(V[[i]]) %*% coord
            D <- dist(coord_star, method = "euclidean", diag = T, upper = T)
            exp(-as.matrix(D))
          }
          else 
            matrix(0, nrow = l, ncol = l)
        }
      }
    }
    else Qstar <- Q
  
  
  while(dist[count] > t & count < maxIt){
    # Different location parameter for each image
    if(!singleQ){
      out <- foreach(i = c(1:nsubj)) %dopar% {
        
        GPASub(Xstar[,,i], Qstar[[i]], k, kQ = NULL, ref_ds, scaling, reflection, centered=F)
        
      }
    }
    # else single parameter for every image
    else{
      out <- foreach(i = c(1:nsubj)) %dopar% {
        
        GPASub(Xstar[,,i], Qstar, k, kQ = NULL, ref_ds, scaling, reflection, centered=F)
        
      }
    }
    count <- count + 1 
    Xest = array(unlist(sapply(c(1:nsubj), function(x) out[[x]]$Xest,simplify = F)), dim = dim(Xstar))
    R = array(unlist(sapply(c(1:nsubj), function(x) out[[x]]$R,simplify = F)), dim = c(l,l,nsubj))
    ref_ds_old = ref_ds
   
    ref_ds <- colMeans(aperm(Xest, c(3, 1, 2)))
    dist[count] <- norm(ref_ds-ref_ds_old, type="F")
  }
  
  Xest1 <- vector(mode = "list", length = nsubj)
  for (i in 1:nsubj)
    Xest1[[i]] <- Xest[,,i] %*% t(V[[i]])
  
  return(list(Xest = Xest1, R = R, dist = dist,  count = count, M = ref_ds))
}

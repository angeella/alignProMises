#' @title Efficient ProMises model
#' @description Performs the functional alignment using the light version of the von Mises Fisher Procrustes model allowing for different number of voxels between matrices
#' @usage EfficientProMises(data, maxIt=10, t =.001, k = 0, Q = NULL, 
#' ref_ds = NULL, scaling = T, reflection= T, subj= T, coord = NULL)
#' @param data data, i.e., list of matrices with dimension time points - voxels. Matrices can have different number of columns but must have same numbers of rows.
#' @param maxIt maximum number of iteration
#' @param t the threshold value to be reached as the minimum relative reduction between the matrices
#' @param k value of the concentration parameter of the prior distribution
#' @param Q value of the location parameter of the prior distribution. It has dimension time-points x time-points, it could be not symmetric.
#' @param ref_ds starting matrix to align
#' @param scaling Flag to apply scaling transformation
#' @param reflection Flag to apply reflection transformation
#' @param subj Flag if each subject has his/her own set of voxel after voxel selection step
#' @param coord list with 3-dim coordinates of the variables. Matrices can have different coordinates. If If the location parameter \code{Q = NULL}, then \code{coord} is used to compute it
#' @author Angela Andreella and Daniela Corbetta
#' @return \code{EfficientProMises} returns a list with four components:
#' \item{\code{Xest}}{a list with the aligned matrices}
#' \item{\code{R}}{an array with the rotation matrices}
#' \item{\code{dist}}{a vector with length equal to the number of iterations that contains the distances between a reference matrix and the previous one}
#' \item{\code{count}}{the number of iterations done by the algorithm}
#' @export
#' @importFrom plyr aaply 
#' @importFrom foreach foreach
#' @importFrom foreach %dopar%
#' @importFrom rARPACK svds

EfficientProMisesSubj <- function(data, maxIt=10, t =.001, k = 0, Q = NULL, ref_ds = NULL, scaling = T, reflection= T, subj= T, coord = NULL){
  
  if(!is.list(data)){warnings("Please insert a list of matrices with dimension time points - voxels")}
  
  row <- sapply(data, nrow)[1]
  col <- sapply(data, ncol)
  nsubj <- length(data)
  
  count = 1
  dist = vector()
  dist[1] <- Inf
  
  k <- row
  
  V <- foreach(i = c(1:nsubj)) %dopar% {
    out <- svds(matrix(unlist(data[[i]]), nrow = k), k = k)
    out$v 
    
  }
  # preparo l'array che contiene le X*
  Xstar <- array(NA, dim=c(k,k,nsubj))
  dim(Xstar)
  # da cambiare
  for (i in 1:nsubj){
    Xstar[,,i] <- matrix(unlist(data[[i]]), nrow = row) %*% matrix(unlist(V[[i]]), nrow = col[i])
  }
  
  ref_ds <- colMeans(aperm(Xstar, c(3, 1, 2)))
  
  Xest <-  array(NA, dim(Xstar))
  R <-  array(NA, c(row, row, nsubj))
  
  # Calcolo Q*=exp(-D)
    if(is.null(Q)){
      if(subj){
        Qstar <- foreach(i = c(1:nsubj)) %dopar% {
          if(!is.null(coord)){
            
            coord_star <- t(V[[i]]) %*% coord[[i]]
            D <- dist(coord_star, method = "euclidean", diag = T, upper = T)
            as.matrix(exp(-D))
          }
          else 
            matrix(0, nrow = row, ncol = row)
        }
        
      }
      
      else{
        Qstar <- foreach(i = c(1:nsubj)) %dopar% {
          if(!is.null(coord)){
            
            coord_star <- t(V[[i]]) %*% coord
            D <- dist(coord_star, method = "euclidean", diag = T, upper = T)
            as.matrix(exp(-D))
          }
          else 
            matrix(0, nrow = row, ncol = row)
        }
      }
    }
    else Qstar <- Q
  
  
  while(dist[count] > t & count < maxIt){
    
    out <- foreach(i = c(1:nsubj)) %dopar% {

        GPASub(Xstar[,,i], Qstar[[i]], k, kQ = NULL, ref_ds, scaling, reflection, centered=F)
      
    }
    count <- count + 1 
    Xest = array(unlist(sapply(c(1:nsubj), function(x) out[[x]]$Xest,simplify = F)), dim = dim(Xstar))
    R = array(unlist(sapply(c(1:nsubj), function(x) out[[x]]$R,simplify = F)), dim = c(row,row,nsubj))
    ref_ds_old = ref_ds
    # ref_ds = aaply(Xest, c(1,2), mean)
    ref_ds <- colMeans(aperm(Xest, c(3, 1, 2)))
    dist[count] <- norm(ref_ds-ref_ds_old, type="F")
  }
  
  Xest1 <- vector(mode = "list", length = nsubj)
  for (i in 1:nsubj)
    Xest1[[i]] <- Xest[,,i] %*% t(matrix(unlist(V[[i]]), nrow = col[i]))
  # for (i in 1:nsubj) Xest1[,,i] <- Xest[,,i] %*% t(V)
  
  return(list(Xest = Xest1, R = R, dist = dist,  count = count))
}




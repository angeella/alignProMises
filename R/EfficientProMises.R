#' @title Efficient ProMises model
#' @description Performs the functional alignment using the Efficient ProMises model. Matrices must have the same number of rows.
#' @usage EfficientProMises(data, maxIt=10, t =.001, k = 0, Q = NULL, 
#' ref_ds = NULL, scaling = T, reflection= T, subj= F, centered = T, coord = NULL)
#' @param data data, i.e., array of matrices with dimension time points - voxels or list of matrices with dimension time points - voxels 
#' @param maxIt maximum number of iteration
#' @param t the threshold value to be reached as the minimum relative reduction between the mean matrices
#' @param k value of the concentration parameter of the prior distribution
#' @param Q value of the location parameter of the prior distribution. It has dimension time-points x time-points, it could be not symmetric.
#' @param ref_ds starting matrix to align
#' @param scaling Flag to apply scaling transformation
#' @param reflection Flag to apply reflection transformation
#' @param subj Flag if each subject has his/her own set of voxel after voxel selection step
#' @param centered center data?
#' @param coord 3-dim or 2-dim coordinates of the variables. If \code{subj = F}  then coord is a matrix with dimensions voxels x 2/3, if \code{subj = T} then coord is a list of matrices 
#' with dimensions voxels x 2/3. If the location parameter \code{Q = NULL}, then \code{coord} is used to compute it 
#' @author Angela Andreella and Daniela Corbetta
#' @return \code{EfficientProMises} returns a list with four components:
#' \item{\code{Xest}}{an array with the aligned matrices}
#' \item{\code{R}}{an array with the rotation matrices}
#' \item{\code{dist}}{a vector with length equal to the number of iterations that contains the distances between a reference matrix and the previous one}
#' \item{\code{count}}{the number of iterations done by the algorithm}
#' @references For the theory on the ProMises model see: A. Andreella and L. Finos
#' (2022), Procrustes analysis for high-dimensional data, 
#' Psychometrika 87, 1422-1438
#' @export
#' @importFrom plyr aaply 
#' @importFrom foreach foreach
#' @importFrom foreach %dopar%
#' @importFrom rARPACK svds

EfficientProMises <- function(data, maxIt=10, t =.001, k = 0, Q = NULL, ref_ds = NULL, scaling = T, reflection= T, subj= F, centered = T, coord = NULL){
  
  if(!is.array(data) & !is.list(data)){warnings("Please insert an array or a list of matrices with dimension time points - voxels")}
  if(is.list(data)){
    data <- array(as.numeric(unlist(data)), dim=c(nrow(data[[1]]), ncol(data[[1]]), length(data)))
  }
  
  
  row <- dim(data)[1] 
  col <- dim(data)[2] 
  nsubj <- dim(data)[3]
  
  count = 1
  dist = vector()
  dist[1] <- Inf
  
  # M <- aaply(data, c(1,2), mean)
  M <- colMeans(aperm(data, c(3, 1, 2)))
  if(centered){
    datas_centered <- aaply(data, 3, function(x) x - M)
    X <- aaply(datas_centered, 1, function(x) x/norm(x,type="F"))
    X <- aperm(X,c(2,3,1))
  }else{
    X <- data
  }
  
  if(is.null(ref_ds)){
    ref_ds <- M 
  }
  
  out <- svds(ref_ds, k = nrow(ref_ds))
  V <- out$v 
  Xstar <- array(NA, dim=c(nrow(X), ncol(V), nsubj))
  
  Xstar[] <- apply(X, 3, function(x) x%*%V)
  ref_ds <- ref_ds %*% V
  
  Xest <-  array(NA, dim(Xstar))
  R <-  array(NA, c(row, row, nsubj))
  
  # Calcolo Qstar
  if(is.null(Q)){
    if(subj){
      Qstar <- foreach(i = c(1:nsubj)) %dopar% {
          if(!is.null(coord)){
            
            coord_star <- t(V) %*% coord[[i]]
            D <- dist(coord_star, method = "euclidean", diag = T, upper = T)
            exp(-as.matrix(D))
          }
          else 
            matrix(0, nrow = row, ncol = row)
        }
        
    }
    
    else{
        if(!is.null(coord)){
          
          coord_star <- t(V) %*% coord
          D = dist(coord_star, method = "euclidean", diag = T, upper = T)
          Qstar <- exp(-as.matrix(D))
        }
        else 
          Qstar <- matrix(0, nrow = row, ncol = row)
      }
  }
  else Qstar <- Q
  
  
  
  while(dist[count] > t & count < maxIt){
    
    out <- foreach(i = c(1:nsubj)) %dopar% {
      
      if(subj){
        
        GPASub(Xstar[,,i], Qstar[[i]], k, kQ = NULL, ref_ds, scaling, reflection, centered)
        #vMFP(X[,,i], k, Q[,,i], ref_ds, scaling, reflection)
      }
      
      else{
        
        GPASub(Xstar[,,i], Qstar, k, kQ = NULL, ref_ds, scaling, reflection, centered)
        #vMFP(X[,,i], k, Q, ref_ds, scaling, reflection) 
      }
      
    }
    count <- count + 1 
    Xest = array(unlist(sapply(c(1:nsubj), function(x) out[[x]]$Xest,simplify = F)), dim = dim(Xstar))
    R = array(unlist(sapply(c(1:nsubj), function(x) out[[x]]$R,simplify = F)), dim = c(row,row,nsubj))
    ref_ds_old = ref_ds
    # ref_ds = aaply(Xest, c(1,2), mean)
    ref_ds <- colMeans(aperm(Xest, c(3, 1, 2)))
    dist[count] <- norm(ref_ds-ref_ds_old, type="F")
  }
  
  Xest1 <- array(NA, c(row, col, nsubj))
  for (i in 1:nsubj) Xest1[,,i] <- Xest[,,i] %*% t(V)
  
  return(list(Xest = Xest1, R = R, dist = dist, count = count))
}

#' @title Efficient ProMises model
#' @description perform the functional alignment using the von Mises Fisher Procrustes model
#' @usage vMFPmodelLight(data, maxIt=10, t=.001, k = 0, Q = NULL, ref_ds = NULL, scaling= T, reflection= T, subj= F, coord = NULL)
#' @param data data, i.e., array of matrices with dimension time points - voxels 
#' @param maxIt maximum number of iteration
#' @param t the threshold value to be reached as the minimum relative reduction between the matrices
#' @param k value of the concentration parameter of the prior distribution
#' @param Q value of the location parameter of the prior distribution. It has dimension voxels x voxels, it could be not symmetric.
#' @param ref_ds starting matrix to align
#' @param scaling Flag to apply scaling transformation
#' @param reflection Flag to apply reflection transformation
#' @param subj Flag if each subject has his/her own set of voxel after voxel selection step
#' @param centered centered data?
#' @param coord 3-dim coordinates of the variabiles
#' @author Angela Andreella
#' @return Returns list of matrices
#' @export
#' @importFrom plyr aaply 
#' @importFrom foreach foreach
#' @importFrom foreach %dopar%
#' @importFrom rARPACK svds

EfficientProMises <- function(data, maxIt=10, t =.001, k = 0, Q = NULL, ref_ds = NULL, scaling = T, reflection= T, subj= F, centered = T, coord = NULL){
  
  if(!is.array(data)){warnings("Please insert an array of matrices with dimension time points - voxels")}
  
  row <- dim(data)[1] 
  col <- dim(data)[2] 
  nsubj <- dim(data)[3]
  
  count = 1
  dist = vector()
  dist[1] <- Inf
  
  M <- aaply(data, c(1,2), mean)
  if(centered){
    datas_centered <- aaply(data, 3, function(x) x - M)
    X <- aaply(datas_centered, 1, function(x) x/norm(x,type="F"))
    X<-aperm(X,c(2,3,1))
  }else{
    X<- data
  }
  if(is.null(Q)){ Q <- matrix(0, nrow = col, ncol = col)
  }
  
  if(is.null(ref_ds)){
    ref_ds <- M
  }
  while(dist[count] > t | count < maxIt){
    Xest <-  array(NA, dim(X))
    R <-  array(NA, c(col,col, nsubj))
    
    
    out <-foreach(i = c(1:nsubj)) %dopar% {
      out <- svds(ref_ds, k = nrow(ref_ds)) 
      V <- out$v
      Xest <- X[,,i] %*% V
      if(subj){

        if(!is.null(coord)){
          coord_star <- t(V) %*% coord
          D = dist(coord_star, method = "euclidean")
          Q[,,i] <- exp(-D)
        }else{
          Q[,,i] <- matrix(0, nrow = row, ncol = row)
        }

        vMFP(X[,,i], k, Q[,,i], ref_ds, scaling, reflection)
      }else{
        if(!is.null(coord)){
          coord_star <- t(V) %*% coord
          D = dist(coord_star, method = "euclidean")
          Q <- exp(-D)
        }else{
          Q <- matrix(0, nrow = row, ncol = row)
        }
        vMFP(X[,,i], k, Q, ref_ds, scaling, reflection) 
      }
      
    }
    count <- count + 1 
    Xest = array(unlist(sapply(c(1:nsubj), function(x) out[[x]]$Xest,simplify = F)), dim = dim(X))
    R = array(unlist(sapply(c(1:nsubj), function(x) out[[x]]$R,simplify = F)), dim = c(col,col,nsubj))
    ref_ds_old = ref_ds
    ref_ds = aaply(Xest, c(1,2), mean)
    dist[count] <- norm(ref_ds-ref_ds_old,type="F")
  }
  
  Xest = sapply(c(1:nsubj), function(x) Xest[x] %*% V)
  
  return(list(Xest = Xest, R = R, dist = dist, count = count))
}
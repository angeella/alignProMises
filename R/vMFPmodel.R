#' @title von Mises Fisher Procrustes model
#' @description perform the functional alignment using the von Mises Fisher Procrustes model
#' @usage vMFPmodel(data, maxIt=10, t=.001, k = 0, Q = NULL, ref_ds = NULL, scaling= T, reflection= T, subj= F)
#' @param data data, i.e., array of matrices with dimension time points - voxels or list of matrices with dimension time points - voxels
#' @param maxIt maximum number of iteration
#' @param t the threshold value to be reached as the minimum relative reduction between the matrices
#' @param k value of the concentration parameter of the prior distribution
#' @param Q value of the location parameter of the prior distribution. It has dimension voxels x voxels, it could be not symmetric.
#' @param ref_ds starting matrix to align
#' @param scaling Flag to apply scaling transformation
#' @param reflection Flag to apply reflection transformation
#' @param subj Flag if each subject has his/her own set of voxel after voxel selection step
#' @param centered centered data?
#' @author Angela Andreella and Daniela Corbetta
#' @return Returns list of matrices
#' @export
#' @importFrom plyr aaply 
#' @importFrom foreach foreach
#' @importFrom foreach %dopar%

vMFPmodel <- function(data, maxIt=10, t =.001, k = 0, Q = NULL, ref_ds = NULL, scaling = T, reflection= T, subj= F, centered = T){
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
    X<-aperm(X,c(2,3,1))
  }else{
    X<- data
  }
  if(is.null(Q)){ Q <- matrix(0, nrow = col, ncol = col)
  }
  if(is.null(ref_ds)){
    ref_ds <- M
  }
  Xest <-  array(NA, dim(X))
  R <-  array(NA, c(col,col, nsubj))

  while(dist[count] > t | count < maxIt){
    
    out <-foreach(i = c(1:nsubj)) %dopar% {
      if(subj){
       # GPASub(X[,,i], Q[,,i], k, ref_ds, scaling, reflection)
        vMFP(X[,,i], Q[,,i], k, ref_ds, scaling, reflection)
      }else{
      #  GPASub(X[,,i], Q, k, ref_ds, scaling, reflection) 
        vMFP(X[,,i], Q, k, ref_ds, scaling, reflection) 
      }
     
    }
    count <- count + 1 
    Xest = array(unlist(sapply(c(1:nsubj), function(x) out[[x]]$Xest,simplify = F)), dim = dim(X))
    R = array(unlist(sapply(c(1:nsubj), function(x) out[[x]]$R,simplify = F)), dim = c(col,col,nsubj))
    ref_ds_old = ref_ds
    #ref_ds = aaply(Xest, c(1,2), mean)
    ref_ds <- colMeans(aperm(Xest, c(3, 1, 2)))
    dist[count] <- norm(ref_ds-ref_ds_old,type="F")
  }

  return(list(Xest = Xest, R = R, dist = dist, count = count))
}

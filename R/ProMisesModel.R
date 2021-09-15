#' @title ProMises model
#' @description Performs the functional alignment using the von Mises Fisher Procrustes model with unknown reference matrix
#' @usage ProMisesModel(data, maxIt=10, t=.001, k = 0, Q = NULL, 
#'        ref_ds = NULL, scaling= T, reflection= T, 
#'        subj= F, centered=T, kCalibrate = F, D = NULL, p = 0.01
#'        ind = 2)
#' @param data data, i.e., array of matrices with dimension time points - voxels or list of matrices with dimension time points - voxels
#' @param maxIt maximum number of iterations
#' @param t the threshold value to be reached as the minimum relative reduction between the matrices
#' @param k value of the concentration parameter of the prior distribution
#' @param Q value of the location parameter of the prior distribution. It has dimension voxels x voxels, it could be not symmetric.
#' @param ref_ds starting reference matrix to align. If \code{NULL}, at the first iteration it is set equal to the element-wise mean of the data matrices
#' @param scaling Flag to apply scaling transformation
#' @param reflection Flag to apply reflection transformation
#' @param subj Flag if each subject has his/her own set of voxel after voxel selection step
#' @param centered if \code{TRUE}, data matrices are centered before the algorithm starts
#' @param kCalibrate if \code{TRUE}, the parameters of the prior distribution are computed with the function \code{kCalibrate}
#' @param D squared euclidean distance matrix between the coordinates of the voxels. Necessary only when \code{kCalibrate = TRUE}
#' @param p distances to compute bandwith value for the \code{kCalibrate} function. Default value is 0.01
#' @param ind (only required in the case with two matrices) the index of the matrix to be used as reference, must be 1 or 2
#' @author Angela Andreella and Daniela Corbetta
#' @return In the case with more than two matrices, \code{ProMisesModel} returns a list with four components:
#' \item{\code{Xest}}{an array with the aligned matrices}
#' \item{\code{R}}{an array with the rotation matrices}
#' \item{\code{dist}}{a vector with length equal to the number of iterations that contains the distances between a reference matrix and the previous one}
#' \item{\code{count}}{the number of iterations done by the algorithm}
#' In the case with only two matrices, \code{ProMisesModel} returns a list with two components:
#' \item{\code{Xest}}{the aligned matrix}
#' \item{\code{R}}{the rotation matrix}
#' @references For the theory on the von Mises-Fisher-Procrustes model see: A. Andreella and L. Finos
#' (2020), The von Mises-Fisher Procrustes model in functional Magnetic Resonance Imaging data, 
#' arXiv: 2008.04631
#' @export
#' @importFrom plyr aaply 
#' @importFrom foreach foreach
#' @importFrom foreach %dopar%

ProMisesModel <- function(data, maxIt=10, t =.001, k = 0, Q = NULL, ref_ds = NULL, scaling = T, reflection= T, subj= F, centered = T, kCalibrate = FALSE, D = NULL, p = 0.01, ind = 2){
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
    X<- data
  }
  if(is.null(Q)){ Q <- matrix(0, nrow = col, ncol = col)
  }
  if(is.null(ref_ds) & nsubj != 2){
    ref_ds <- M
  }
  if(kCalibrate){
    kQ <- kCalibrate(D = D, p = p)
  }
  else{
    kQ <- NULL
  }
  
  # 2 images: explicit solution
  if (nsubj == 2 & is.null(ref_ds)){
    if(ind != 1 & ind != 2){warnings("ind must be 1 or 2")}
    ref_ds <- X[,,ind]
    out <- ProMises::GPASub(X[,,-ind], Q, k, kQ, ref_ds, scaling, reflection, centered)
    Xest <- array(NA, dim = dim(X))
    Xest[,,ind] <- data[,,ind]
    Xest[,,-ind] <- out$Xest
    R <- array(NA, dim = c(col, col, 2))
    R[,,ind] <- diag(col)
    R[,,-ind] <- out$R
    return(list(Xest = Xest, R = R))
  }
  
  else{
  
  Xest <-  array(NA, dim(X))
  R <-  array(NA, c(col,col, nsubj))
  
  while(dist[count] > t & count < maxIt){
    
    out <-foreach(i = c(1:nsubj)) %dopar% {
      if(subj){
        # GPASub(X[,,i], Q[,,i], k, ref_ds, scaling, reflection)
        ProMises::GPASub(X[,,i], Q[,,i], k, kQ, ref_ds, scaling, reflection, centered)
      }else{
        #  GPASub(X[,,i], Q, k, ref_ds, scaling, reflection) 
        ProMises::GPASub(X[,,i], Q, k, kQ, ref_ds, scaling, reflection, centered) 
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
}

#' @title ProMises model
#' @description Performs the functional alignment of input data matrices using the ProMises model with unknown reference matrix.
#' @param data data, i.e., array of matrices with dimension time points - voxels or list of matrices with dimension time points - voxels.
#' @param maxIt maximum number of iterations.
#' @param t the threshold value to be reached as the minimum relative reduction between the mean matrices.
#' @param k value of the concentration parameter of the prior distribution of the rotation parameter.
#' @param Q value of the location parameter of the prior distribution of the rotation parameter. It has dimension voxels x voxels, it could be not symmetric. If \code{Q=NULL}, it is set equal to a matrix of 0 (i.e. no regularization, classic Generalized Procrustes Analysis).
#' @param ref_ds starting reference matrix to align. If \code{NULL}, at the first iteration it is set equal to the element-wise mean of the data matrices.
#' @param scaling Flag to apply scaling transformation.
#' @param reflection Flag to apply reflection transformation.
#' @param subj Flag if each subject has his/her own set of voxel after voxel selection step. If \code{subj=TRUE}, then the model considers a different location parameter for the prior distribution of each rotation parameter and \code{Q} has to be an array of dimensions voxels x voxels x number of matrices.
#' @param center Flag to apply centering transformation to data prior to alignment step.
#' @param kCalibrate if \code{TRUE}, the parameters of the prior distribution are computed with the function \code{kCalibrate}.
#' @param D squared euclidean distance matrix between the coordinates of the voxels. Necessary only when \code{kCalibrate = TRUE}.
#' @param p distances to compute bandwith value for the \code{kCalibrate} function. Default value is 0.01.
#' @param ind (only required in the case with two matrices) the index of the matrix to be used as reference, must be 1 or 2.
#' @author Angela Andreella and Daniela Corbetta
#' @return In the case with more than two matrices, \code{ProMisesModel} returns a list with four components:
#' \item{\code{Xest}}{an array with the aligned matrices}
#' \item{\code{R}}{an array with the rotation matrices}
#' \item{\code{dist}}{a vector with length equal to the number of iterations that contains the distances between a reference matrix and the previous one}
#' \item{\code{count}}{the number of iterations done by the algorithm}
#' In the case with only two matrices, \code{ProMisesModel} returns a list with two components:
#' \item{\code{Xest}}{an array with the aligned matrices}
#' \item{\code{R}}{an array with the rotation matrices (one is the identity)}
#' @details 
#' Let \eqn{X_1,\dots, X_N \in \mathbb{R}^{n\times m}} be the \eqn{N} matrices to be aligned. 
#' The ProMises model assumes  that each matrix can be represented as a random perturbation of a common
#' reference matrix plus an error term:
#' \deqn{X_i=(M+E_i)R_i^\top,} where
#' \itemize{
#' \item \eqn{M \in \mathbb{R}^{n\times m}} is the common reference matrix,
#' \item \eqn{E_i \in \mathbb{R}^{n\times m}} is a matrix error term following a matrix normal distribution,
#' \item \eqn{R_i \in \mathbb{R}^{m\times m}} is the rotation parameter and follows a von Mises-Fisher
#' distribution with location parameter \eqn{Q} and concentration parameter \eqn{k}.
#' }
#' 
#' The rotation parameters \eqn{R_i} are estimated using the maximum a posteriori estimate,
#' given by \eqn{\hat{R}_i = U_iV_i^\top}, where \eqn{U_i} and \eqn{V_i} are obtained
#' from the singular value decomposition of \eqn{Q^* = X_i^\top M + kQ}.
#' 
#' The common reference matrix \eqn{M} represents the shared structure among all matrices,
#' serving as a baseline for the alignment process. Since it is unknown, it is iteratively updated during the
#' algorithm execution. In particular, the algorithm starts with an initial guess for \eqn{M} 
#' as the element-wise mean matrix of the input matrices \eqn{X_i}. Then, in each iteration, \eqn{M} is updated based on the mean of the
#' aligned matrices until convergence (Frobenius distance between two subsequent 
#' estimates of \eqn{M} less than \code{t}) or until the maximum number of iterations is reached.
#'
#' 
#' If the number of columns \eqn{m} is large or if each matrix has a different number of column,
#' then use the \link[alignProMises]{EfficientProMises} function.
#' 
#' If there are only two matrices, then one matrix is projected onto the other. The solution 
#' is explicit (see Green, 1952) and is equal to \eqn{\hat{R}=UV^\top},
#' where \eqn{U} and \eqn{V} come from the SVD of \eqn{X_1^\top X_2}.
#' @references For the theory on the ProMises model and more details on the algorithm see: 
#' A. Andreella and L. Finos (2022), Procrustes analysis for high-dimensional data, 
#' Psychometrika 87, 1422-1438
#' 
#' For the two matrix case see: B. Green (1952), The orthogonal approximation of an oblique structure in
#' factor analysis, Psychometrika 17, 429-440.
#' @examples{
#' ## create a random array of 3 matrices with 24576 time-points and 60 voxels
#' data<- array(rnorm(24576*60*3), dim = c(24576,60,3))
#' ## Align the three matrices setting the location parameter equal to the 
#' ## identity and the concentration parameter equal to 1
#' out <- ProMisesModel(data, maxIt = 20, t = 1/100, k = 1, Q = diag(1,60), 
#' scaling = FALSE, reflection = FALSE, subj = FALSE, center = FALSE)
#' }
#' @export
#' @importFrom plyr aaply 
#' @importFrom foreach foreach
#' @importFrom foreach %dopar%

ProMisesModel <- function(data, maxIt=10, t =.001, k = 0, Q = NULL, ref_ds = NULL, scaling = TRUE, reflection= TRUE, subj = FALSE, center = TRUE, kCalibrate = FALSE, D = NULL, p = 0.01, ind = 2){
  if(!is.array(data) & !is.list(data)){
    stop("Please insert an array or a list of matrices with dimension time points - voxels")
  }
  if(is.list(data)){
    data <- array(as.numeric(unlist(data)), dim=c(nrow(data[[1]]), ncol(data[[1]]), length(data)))
  }
  
  row <- dim(data)[1] 
  col <- dim(data)[2] 
  nsubj <- dim(data)[3]
  
  if(subj & (is.na(dim(Q)[3]) | dim(Q)[3]!=nsubj)){
    stop("If subj=T then please provide location parameter for each matrix (i.e. Q must be an array)")
  }
  
  if(col>=4000){
    warning("High dimensions, consider switching to EfficientProMises")
  }
  
  count = 1
  dist = vector()
  dist[1] <- Inf
  
  M <- colMeans(aperm(data, c(3, 1, 2)))
  if(center){
    datas_centered <- aaply(data, 3, function(x) x - M)
    X <- aaply(datas_centered, 1, function(x) x/norm(x,type="F"))
    X <- aperm(X,c(2,3,1))
  } 
  else{
    X<- data
  }
  if(is.null(Q)){ 
    Q <- matrix(0, nrow = col, ncol = col)
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
    out <- GPASub(X[,,-ind], Q, k, kQ, ref_ds, scaling, reflection, center)
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
  i <- NULL
  
  while(dist[count] > t & count < maxIt){
    
    out <- foreach(i = c(1:nsubj)) %dopar% {
      if(subj){
        GPASub(X[,,i], Q[,,i], k, kQ, ref_ds, scaling, reflection, center)
      }else{
        GPASub(X[,,i], Q, k, kQ, ref_ds, scaling, reflection, center) 
      }
      
    }
    count <- count + 1 
    Xest = array(unlist(sapply(c(1:nsubj), function(x) out[[x]]$Xest,simplify = F)), dim = dim(X))
    R = array(unlist(sapply(c(1:nsubj), function(x) out[[x]]$R,simplify = F)), dim = c(col,col,nsubj))
    ref_ds_old = ref_ds
    ref_ds <- colMeans(aperm(Xest, c(3, 1, 2)))
    dist[count] <- norm(ref_ds-ref_ds_old,type="F")
  }
  
  return(list(Xest = Xest, R = R, dist = dist, count = count))
  }
}

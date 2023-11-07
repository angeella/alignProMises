#' @title Efficient ProMises model
#' @description Performs the functional alignment using the Efficient ProMises model allowing for different number of columns (voxel or pixel) between matrices.
#' @param data data, i.e., list of matrices with dimension time points - voxels. Matrices can have different number of columns but must have the same numbers of rows.
#' @param maxIt maximum number of iteration
#' @param t the threshold value to be reached as the minimum relative reduction between the mean matrices
#' @param k value of the concentration parameter of the prior distribution
#' @param Q value of the location parameter of the prior distribution. It has dimension time-points x time-points, it could be not symmetric.
#' @param ref_ds starting matrix to align. If \code{NULL}, it is set to the element-wise arithmetic mean of the reduced matrices
#' @param scaling Flag to apply scaling transformation
#' @param reflection Flag to apply reflection transformation
#' @param subj Flag if each subject has his/her own set of voxel after voxel selection step. If \code{TRUE}, then \code{coord} has to be a list with each subjects' coordinates and the model will assume a different location parameter for the prior distribution of the rotation parameters. 
#' @param coord list with 3-dim or 2-dim coordinates of the voxels of dimensions voxels x 2/3. Matrices can have different coordinates. If the location parameter \code{Q = NULL}, then \code{coord} is used to compute it (see Details)
#' @param l number of singular vectors required to obtained the reduced transformation of the matrices. If \code{NULL}, it is set equal to the number of rows. 
#' @author Angela Andreella and Daniela Corbetta
#' @details
#' \bold{Location parameter calculation:} 
#' If \code{Q=NULL} and \code{coord=NULL}, then \code{Q} is set to a matrix of zeros (i.e. no regularization, standard Generalized Procrustes Analysis).
#' If \code{Q=NULL} but \code{coord} is not null, then \code{Q} is calculated as follows:
#' \itemize{
#' \item For each input matrix \eqn{X_i \in \mathbb{R^{n\times m}}}, compute the thin-SVD \eqn{X_i=U_iD_iV_i^\top}, with \eqn{V_i\in \mathbb{R}^{m\times n}}
#' \item pre-multiply each coordinates matrix \eqn{C_i} (same for all the matrices if \code{subj=F}) by \eqn{V_i^\top}: \eqn{C_i^*=V_i^\top C_i}
#' \item compute the euclidean distance matrix \eqn{D} among points using the coordinates in \eqn{C_i^*}
#' \item set \eqn{Q=\exp\{-D\}}}
#' @return \code{EfficientProMisesSubj} returns a list with five components:
#' \item{\code{Xest}}{a list with the aligned matrices}
#' \item{\code{R}}{an array with the rotation matrices}
#' \item{\code{dist}}{a vector with length equal to the number of iterations that contains the distances between a reference matrix and the previous one}
#' \item{\code{count}}{the number of iterations done by the algorithm}
#' \item{\code{M}}{the element-wise mean matrix of the aligned matrices.}
#' @examples{
#' ## Create random list of matrices with different number of columns
#' X <- list(matrix(rnorm(100*4000), nrow=100),
#'           matrix(rnorm(100*3669), nrow=100),
#'           matrix(rnorm(100*3500), nrow=100))
#' ## Align the matrices with the Efficient ProMises model. 
#' ## Since subj = F and Q is a single matrix, all the rotation parameters 
#' ## will have the same location parameter
#' out <- EfficientProMisesSubj(data = X, maxIt = 10, t = 1, k = 1,
#'                              Q = diag(1,100), subj = FALSE, scaling = FALSE,
#'                              reflection = FALSE, ref_ds = NULL, coord = NULL, 
#'                              l = NULL)
#'                              
#' ## create random coordinates
#' C <- list(cbind(sample(1:4000), sample(1:4000)),
#'           cbind(sample(1:3669), sample(1:3669)),
#'           cbind(sample(1:3500), sample(1:3500)))
#'
#' ## Align matrices considering a different location parameter for every 
#' ## matrix, calculated from the coordinates as specified in details
#' out1 <- EfficientProMisesSubj(data = X, maxIt = 10, t = 1, k = 1,
#'                               Q = NULL, subj = TRUE, scaling = FALSE,
#'                               reflection = FALSE, ref_ds = NULL, coord = C, 
#'                               l = NULL)   
#'                               
#' ## Extract only the first 15 singular vectors to compute the low-dimensional
#' ## representation of the matrices
#' out2 <- EfficientProMisesSubj(data = X, maxIt = 10, t = 1, k = 1,
#'                               Q = NULL, subj = TRUE, scaling = FALSE,
#'                               reflection = FALSE, ref_ds = NULL, coord = C, 
#'                               l = 15)                         
#' }
#' @references For the theory on the Efficient ProMises model see: A. Andreella 
#' and L. Finos (2022), Procrustes analysis for high-dimensional data, 
#' Psychometrika 87, 1422-1438
#' @export
#' @importFrom plyr aaply 
#' @importFrom foreach foreach
#' @importFrom foreach %dopar%
#' @importFrom rARPACK svds

EfficientProMisesSubj <- function(data, maxIt=10, t =.001, k = 0, Q = NULL, ref_ds = NULL, scaling = TRUE, reflection= TRUE, subj= TRUE, coord = NULL, l = NULL){
  
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
      if(subj){ #different Q for every rotation parameter
        Qstar <- foreach(i = c(1:nsubj)) %dopar% {
          if(!is.null(coord)){
            # if Q is not pre-specified and we have coordinates, set Q*=exp(-dist(Q^t C))
            coord_star <- t(V[[i]]) %*% coord[[i]]
            D <- dist(coord_star, method = "euclidean", diag = T, upper = T)
            exp(-as.matrix(D))
          }
          else 
            # else set it to 0, no regularization
            matrix(0, nrow = l, ncol = l)
        }
        
      }
      
      else{ #same coordinates for all the matrices, only one location parameter
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
    if(subj){
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

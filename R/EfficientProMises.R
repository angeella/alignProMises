#' @title Efficient ProMises model
#' @description Performs the functional alignment of input data matrices using the Efficient ProMises model. This model has a reduced computational load compared to the full ProMises model and it allows input matrices to have different number of columns.
#' @param data data, i.e., list or array of matrices with dimension time points - voxels. Matrices can have different number of columns but must have the same numbers of rows.
#' @param maxIt maximum number of iterations.
#' @param t the threshold value to be reached as the minimum relative reduction between the mean matrices.
#' @param k value of the concentration parameter of the prior distribution of the rotation parameter.
#' @param Q value of the location parameter of the prior distribution of the rotation parameter. It has dimension time-points x time-points, it could be not symmetric.
#' @param ref_ds starting reference matrix to align. If \code{NULL}, it is set to the element-wise arithmetic mean of the reduced matrices.
#' @param scaling Flag to apply scaling transformation.
#' @param reflection Flag to apply reflection transformation.
#' @param center Flag to apply centering transformation to reduced data prior to alignment step.
#' @param subj Flag if each subject has his/her own set of voxel after voxel selection step. If \code{TRUE}, then \code{coord} has to be a list with each subjects' coordinates and the model will assume a different location parameter for the prior distribution of the rotation parameters. 
#' @param coord list with 3-dim or 2-dim coordinates of the voxels of dimensions voxels x 2/3. Matrices can have different coordinates. If the location parameter \code{Q = NULL}, then \code{coord} is used to compute it (see Details).
#' @param version which version of the Efficient ProMises model to use. If \code{version = "M"},
#' then it decomposes the element-wise mean matrix (default), if \code{version = "X"}, it decomposes the
#' single input matrices. See Details for the description of the two versions.
#' @param l number of singular vectors required to obtained the reduced transformation of the matrices. If \code{NULL}, it is set equal to the number of rows. 
#' @author Angela Andreella and Daniela Corbetta
#' @details
#' Let \eqn{X_1,\dots, X_N \in \mathbb{R}^{n\times m}} be the \eqn{N} matrices to be aligned.
#' The Efficient ProMises model computes a low dimensional representation of the \eqn{X_i} prior
#' to the alignment step, achieving computational efficacy without loss of information
#' (see Andreella and Finos, 2022 for details).
#' 
#' The function allows for two versions of dimensionality reduction:
#' \itemize{
#' \item \code{version = "M"}: If all the matrices have the same dimensions, compute \eqn{\hat{M}=\sum_{i=1}^N X_i / N},
#' derive the thin-SVD of \eqn{\hat{M}=UDV^\top}, with \eqn{V \in \mathbb{R}^{m\times n}},
#' and multiply each matrix by \eqn{V}. Then apply the ProMises model to \eqn{X_i^*=X_iV}.
#' This version is the default since it saves time.
#' \item \code{version = "X"}: If the matrices have a different number of columns, compute the thin-SVD of each \eqn{X_i=U_iD_iV_i^\top}
#' and multiply each \eqn{X_i} by \eqn{V_i}. Then apply the ProMises model to \eqn{X_i^*=X_iV_i}.
#' }
#' In both cases, at the end of the alignment step, the aligned matrices are projected
#' back to the original space with the inverse transformation \eqn{V^\top} (or \eqn{V_i^\top}).
#' 
#' \bold{Prior location parameter calculation:} 
#' 
#' If \code{Q=NULL} and \code{coord=NULL}, then \code{Q} is set to a matrix of zeros (i.e. no regularization, standard Generalized Procrustes Analysis).
#' If \code{Q=NULL} but \code{coord} is not null, then \code{Q} is calculated as follows
#' when \code{version="X"}:
#' \itemize{
#' \item For each input matrix \eqn{X_i \in \mathbb{R^{n\times m}}}, compute the thin-SVD \eqn{X_i=U_iD_iV_i^\top}, with \eqn{V_i\in \mathbb{R}^{m\times n}};
#' \item pre-multiply each coordinates matrix \eqn{C_i} (same coordinates matrix for all the matrices if \code{subj=F}) by \eqn{V_i^\top}: \eqn{C_i^*=V_i^\top C_i};
#' \item compute the euclidean distance matrix \eqn{D_i} among points using the coordinates in \eqn{C_i^*};
#' \item set \eqn{Q_i=\exp\{-D_i\}}.}
#' If \code{version="M"}, the coordinates matrix are multiplied by \eqn{V} that comes from the thin-SVD of \eqn{\hat{M}}.
#' @return \code{EfficientProMisesSubj} returns a list with five components:
#' \item{\code{Xest}}{a list with the aligned matrices}
#' \item{\code{R}}{an array with the rotation matrices}
#' \item{\code{dist}}{a vector with length equal to the number of iterations that contains the distances between a reference matrix and the previous one}
#' \item{\code{count}}{the number of iterations done by the algorithm}
#' \item{\code{M}}{the element-wise mean matrix of the reduced aligned matrices}
#' @examples{
#' ## Create random list of matrices with different number of columns
#' X <- list(matrix(rnorm(100*4000), nrow=100),
#'           matrix(rnorm(100*3669), nrow=100),
#'           matrix(rnorm(100*3500), nrow=100))
#' ## Align the matrices with the Efficient ProMises model. 
#' ## Since subj = F and Q is a single matrix, all the rotation parameters 
#' ## will have the same location parameter
#' out <- EfficientProMises(data = X, maxIt = 10, t = 1, k = 1,
#'                          Q = diag(1,100), subj = FALSE, scaling = FALSE,
#'                          center = FALSE, reflection = FALSE, 
#'                          ref_ds = NULL, coord = NULL, version="X", 
#'                          l = NULL)
#'                              
#' ## create random coordinates
#' C <- list(cbind(sample(1:4000), sample(1:4000)),
#'           cbind(sample(1:3669), sample(1:3669)),
#'           cbind(sample(1:3500), sample(1:3500)))
#'
#' ## Align matrices considering a different location parameter for every 
#' ## matrix, calculated from the coordinates as specified in details
#' out1 <- EfficientProMises(data = X, maxIt = 10, t = 1, k = 1,
#'                           Q = NULL, subj = TRUE, scaling = FALSE, 
#'                           center = FALSE, reflection = FALSE, 
#'                           ref_ds = NULL, coord = C, 
#'                           version="X", l = NULL)   
#'                               
#' ## Extract only the first 15 singular vectors to compute the low-dimensional
#' ## representation of the matrices
#' out2 <- EfficientProMises(data = X, maxIt = 10, t = 1, k = 1,
#'                           Q = NULL, subj = TRUE, scaling = FALSE, 
#'                           center = FALSE, reflection = FALSE, 
#'                           ref_ds = NULL, coord = C, 
#'                           version="X", l = 15)                         
#' }
#' @references For the theory on the Efficient ProMises model see: A. Andreella 
#' and L. Finos (2022), Procrustes analysis for high-dimensional data, 
#' Psychometrika 87, 1422-1438
#' @export
#' @importFrom plyr aaply 
#' @importFrom foreach foreach
#' @importFrom foreach %dopar%
#' @importFrom rARPACK svds

EfficientProMises <- function(data, maxIt = 10, t = 0.001, k = 0, Q = NULL, ref_ds = NULL, scaling = TRUE, reflection= TRUE, center = TRUE, subj = FALSE, coord = NULL, version = "M", l = NULL){
  
  if(!is.list(data) & !is.array(data)){
    stop("Please insert an array or a list of matrices with dimension time points - voxels")
  }
  if(subj & is.null(Q) & !is.list(coord)){
    stop("If subj=T please provide coordinates for each matrix")
  }
  if(!subj & is.null(Q) & is.list(coord)){
    stop("If subj=F please provide a unique coordinate matrix")
  }
  if(is.list(data)){
    if(version=="M" & !all(sapply(data, function(mat) ncol(mat) == ncol(data[[1]])))){
      stop("Input matrices do not have the same number of columns. Use version = 'X' or subset input matrices")
    }
  }
  if(!is.null(Q) & !is.null(coord)){
    warning("Q is provided, coord will be ignored")
  }
  
  # If version=M work with arrays
  if (is.list(data) & version=="M"){
    data <- array(as.numeric(unlist(data)), dim=c(nrow(data[[1]]), ncol(data[[1]]), length(data)))
  } 
  
  # If version=X work with lists
  if (is.array(data) & version=="X"){
    data <- lapply(1:dim(data)[3], function(i) data[, , i])
  }
  
  if (is.list(data)){
    row <- sapply(data, nrow)[1]
    col <- sapply(data, ncol)
    nsubj <- length(data)
  }
  
  if (is.array(data)){
    row <- dim(data)[1] 
    col <- dim(data)[2] 
    nsubj <- dim(data)[3]
  }
  
  j <- row
  
  # If not otherwise specified, consider full Efficient model 
  if(is.null(l))
    l <- row
  
  Xstar <- array(NA, dim=c(j,l,nsubj))
  
  # if version=M, compute mean matrix, decompose it and obtain reduced Xi*
  if(version=="M"){
    if(is.null(ref_ds)){
      M <- colMeans(aperm(data, c(3, 1, 2)))
    }
    else M <- ref_ds
    # suppress Warning messages:
    #  In fun(A, k, nu, nv, opts, mattype = "matrix") :
    #  all singular values are requested, svd() is used instead
    out <- suppressWarnings(svds(M, k = l))
    V <- out$v
    Xstar[] <- apply(data, 3, function(x) x%*%V)
    # ref_ds <- ref_ds %*% V
    
    # Compute Q*=exp(-D)
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
  }
  
  # else decompose single matrices
  if(version=="X"){
    V <- foreach(i = c(1:nsubj)) %dopar% {
      # suppress Warning messages:
      #  In fun(A, k, nu, nv, opts, mattype = "matrix") :
      #  all singular values are requested, svd() is used instead
      out <- suppressWarnings(svds(data[[i]], k = l))
      out$v 
    }
    for (i in 1:nsubj){
      Xstar[,,i] <- data[[i]] %*% V[[i]] 
    }
    # ref_ds <- colMeans(aperm(Xstar, c(3, 1, 2)))
    
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
      
      else{ #same coordinates for all the matrices
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
  }
  
  ref_ds <- colMeans(aperm(Xstar, c(3, 1, 2))) # Mean of reduced matrices as starting M matrix
  
  # Centering transformation
  if(center){
    datas_centered <- aaply(Xstar, 3, function(x) x - ref_ds)
    X <- aaply(datas_centered, 1, function(x) x/norm(x, type="F"))
    Xstar <- aperm(X, c(2,3,1))
  }
  
  Xest <-  array(NA, dim(Xstar))
  R <-  array(NA, c(l, l, nsubj))
  count = 1
  dist = vector()
  dist[1] <- Inf
  
  singleQ <- ifelse(!is.list(Qstar), TRUE, FALSE)

  while(dist[count] > t & count < maxIt){
    # Different location parameter for each image
    if(!singleQ){
      out <- foreach(i = c(1:nsubj)) %dopar% {
        ProMises(Xstar[,,i], k=k, Q=Q[[i]], ref_ds=ref_ds, scaling=scaling, 
                 reflection=reflection, centered=center)
      }
    }
    # else single parameter for every image
    else{
      out <- foreach(i = c(1:nsubj)) %dopar% {
        ProMises(Xstar[,,i], k=k, Q=Q, ref_ds=ref_ds, scaling=scaling, 
                 reflection=reflection, centered=center) 
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
  {
    if(version=="M") Xest1[[i]] <- Xest[,,i] %*% t(V)
    else Xest1[[i]] <- Xest[,,i] %*% t(V[[i]])
  }
  
  return(list(Xest = Xest1, R = R, dist = dist,  count = count, M = ref_ds))
}

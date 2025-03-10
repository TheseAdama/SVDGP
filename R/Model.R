#' SVDGP: Singular Value Decomposition and Gaussian Process Modeling of Spatio-Temporal Functions
#'
#' This function performs Singular Value Decomposition (SVD) on a centered spatio-temporal observation matrix and independently models the spatial components using Gaussian Processes.
#'
#' @param D A matrix representing spatial data (design matrix).
#' @param FD A matrix representing spatio-temporal observations. FD must have the same number of rows as D.
#' @param formula An optional object of class \code{"formula"} specifying the linear trend for the kriging model (see \code{lm}).
#' The formula should include only the input variables. The default is \code{~1}, which defines a constant trend.
#' @param K The number of components to retain. If not specified, the default is automatically determined to preserve 99% of the total variance.
#' @param spcovtype A vector specifying the covariance function types used for Gaussian Process modeling of the \code{K} spatial components (default: \code{rep("matern5_2", K)}).
#' @param typekrig The type of kriging applied (e.g., "SK" for Simple Kriging; default: "SK").
#' @param ... Additional parameters to be passed to the \code{\link[DiceKriging]{km}} function from the \code{DiceKriging} package.
#' @importFrom DiceKriging km
#' @seealso \code{\link[DiceKriging]{km}}
#' @return A list with the following components:
#' \item{GP}{A list of Gaussian Process models, one for each spatial component, fitted using \code{DiceKriging}.}
#' \item{U}{The matrix containing the left singular vectors (spatial modes) from the SVD. Each column represents a spatial component.}
#' \item{V}{The matrix containing the right singular vectors (temporal modes) from the SVD. Each column represents a temporal component.}
#' \item{Sigma}{The diagonal matrix of singular values from the SVD, indicating the contribution of each component to the decomposition.}
#' \item{Err}{The relative reconstruction error between the original spatio-temporal observations and the approximated matrix obtained through the decomposition.}
#' \item{K}{The number of components retained, either provided by the user or determined based on the 99% variance threshold.}
#' \item{type}{The type of kriging used for Gaussian Process modeling (e.g., "SK" for Simple Kriging).}
#' \item{FDmean}{The column-wise average of the spatio-temporal observation matrix.}
#'
#' @import DiceKriging
#' @export
#' @examples
#' # Test function
#' fxt <- function(x, t) { sin(2 * pi * t * (x - t)) }
#'
#' # Design matrix
#' D <- matrix(seq(0, 1, length = 5), ncol = 1)
#'
#' # Discretized time points
#' tt <- seq(0, 1, length = 15)
#'
#' # Spatio-temporal observations matrix
#' FD <- apply(array(tt), 1, function(t) {
#'   apply(D, 1, function(x) { fxt(x, t) })
#' })
#'
#' # SVD and Gaussian Process modeling
#' model <- svdgppmodel(D, FD, formula = ~1, K = 4, spcovtype = "matern5_2",  typekrig = "SK")
#'
#' # Extract Gaussian Process models
#' model$GP
#'
#' # Plot spatial and temporal modes
#' plot(D[, 1], model$U[, 1], type = 'l', col = "blue", lwd = 2,
#'      xlab = 'x', ylab = "U1")
#' plot(tt, model$V[, 1], type = 'l', col = "red", lwd = 2,
#'      xlab = 't', ylab = "V1")
#'
#' # Print reconstruction error
#' print(model$Err)
svdgppmodel <- function(D, FD,  formula = ~1, K=NULL,
                        spcovtype=NULL,
                        typekrig = "SK", ...){

  if(!is.matrix(D)) stop("D must be a matrix!")
  if(!is.matrix(FD)) stop("FD must be a matrix!")
  if(nrow(FD)!=nrow(D)) stop("FD must have the same number of row as D!")
  kmcovtype <- c("gauss", "matern5_2", "matern3_2", "exp", "powexp")
  if (!inherits(formula, "formula")) stop("`formula` must be an object of class 'formula'!")

  FDmu = 0#colMeans(FD)
  FDcentred = FD - FDmu
  R = svd(FDcentred)
  Sigma = diag(R$d)
  Lambda = diag(Sigma)^2
  Lnorm = (Lambda - min(Lambda)) / (max(Lambda) - min(Lambda))

  if(is.null(K)){
    percentage=0.99
    K = min(which(cumsum(Lnorm) > percentage))
  }

  if(is.null(spcovtype)){
    spcovtype = rep("matern5_2", K)
  }else{
    if(length(spcovtype)==1){
      spcovtype = rep(spcovtype, K)
    }
  }

  for(k in 1:K){
    if(!(spcovtype[k] %in% kmcovtype)){
      stop('All spcovtype must be in c("gauss", "matern5_2", "matern3_2", "exp", "powexp")')
      }
  }

  Sigma = diag(R$d)[1:K, 1:K]
  U = matrix(R$u[,1:K], ncol=K, byrow=FALSE)
  V = matrix(R$v[,1:K], ncol=K, byrow=FALSE)
  FDapp = U %*% Sigma %*% t(V) + FDmu
  Err = norm(FD - FDapp, type = "F") / norm(FD, type = "F")

  Yk = U %*% Sigma
  GP = list()
  for (k in 1:K) {
    YD = matrix(Yk[,k], ncol=1)
    GP[[k]] <- DiceKriging::km(
      formula = formula,
      design = data.frame(D=D),
      response = data.frame(YD=YD),
      covtype = spcovtype[k],
      ...)
  }

  R = list(FDmean=FDmu, GP = GP, U = U, V = V, Sigma = Sigma,
           Err = Err, K = K, type = typekrig)
  return(R)
}

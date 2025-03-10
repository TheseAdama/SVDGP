#' SVDpredict: Spatio-Temporal Prediction Using SVD-GP Model
#'
#' This function predicts spatio-temporal values using a trained SVD-GP model. It provides predictions for the spatio-temporal field, spatial components, temporal components, and associated covariance matrices.
#'
#' @param model A trained SVD-GP model obtained using the \code{svdgppmodel()} function.
#' @param Dpred A matrix of new spatial locations where predictions are to be made.
#' @param computecov Logical. If \code{TRUE}, computes the posterior covariance matrix for the predictions (default: \code{FALSE}).
#'
#' @return A list containing the following elements:
#' \item{Mn}{A matrix of predicted spatio-temporal values. Each row corresponds to a spatial prediction location, and each column represents a temporal point.}
#' \item{Kn}{The posterior covariance matrix for the predicted spatio-temporal values (only included if \code{computecov = TRUE}).}
#' \item{Coefpred}{Predictions of the coefficients for the \code{K} spatial components modeled by Gaussian Processes. Each column corresponds to one spatial component.}
#' \item{spcov}{A list of covariance matrices for each predicted spatial component. Each entry corresponds to one of the \code{K} components.}
#' \item{tempcov}{The temporal covariance matrix derived from the singular value decomposition. This matrix captures the temporal relationships between the modes.}
#' @import DiceKriging
#' @export
#' @examples
#' # Example usage:
#' # Test function
#' fxt <- function(x, t) { sin(2 * pi * t * (x - t)) }
#'
#' # Design matrix
#' D <- matrix(seq(0, 1, length.out = 5), ncol = 1)
#'
#' # Discretized time points
#' tt <- seq(0, 1, length.out = 15)
#'
#' # Spatio-temporal observations matrix
#' FD <- apply(array(tt), 1, function(t) {
#'   apply(D, 1, function(x) { fxt(x, t) })
#' })
#'
#' # Perform SVD and Gaussian Process modeling
#' model <- svdgppmodel(D, FD, formula = ~1, K = 4, spcovtype = "matern5_2",typekrig = "SK")
#' Dpred <- matrix(runif(20), ncol=1)  # New spatial locations
#' Yp <- svdgppredict(model, Dpred, computecov=TRUE)
#' print(Yp$Mn)
#' print(Yp$Kn)
svdgppredict <- function(model, Dpred, computecov=FALSE){

  GPm = model$GPmu
  GPmodel = model$GP
  Sigma = model$Sigma
  V = model$V
  K = model$K
  type = model$type

  Sigmat = V %*% t(V)
  Ypred = matrix(0, nrow=nrow(Dpred), ncol=K)
  Kx = list()
  for(k in 1:K){
    Yp = DiceKriging::predict(GPmodel[[k]], newdata = data.frame(D=Dpred),
                              type = type, se.compute = TRUE, cov.compute = TRUE, checkNames = FALSE)
    Ypred[,k] <- Yp$mean
    Kx[[k]] = Yp$cov
  }

  Mn = Ypred %*% t(V) + model$FDmean

  if(computecov){
    n = nrow(Dpred)
    nt = ncol(Sigmat)
    Kn = matrix(0, ncol=nt*n, nrow=n*nt)
    for(k in 1:K){
      Kn = Kn + kronecker(Kx[[k]], Sigmat)
    }
  }else{
    Kn = NULL
  }

  R = list(Mn = Mn, Kn=Kn, Coefpred = Ypred, spcov = Kx, tempcov = Sigmat)
  return(R)
}

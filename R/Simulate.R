#' Simulate spatio-temporal observations from a SVDGP model
#'
#' This function simulates spatio-temporal data based on a SVDGP model.
#' It draws multiple realizations from the Gaussian Process for each spatial component and reconstructs
#' the full spatio-temporal data using the decomposition.
#'
#' @param model A SVDGP object, typically the output from the `svdgppmodel()` function. It should contain
#' the Gaussian Process models, singular value decomposition matrices (`U`, `V`, and `Sigma`),
#' and the number of components (`K`).
#' @param D A matrix of new spatial data for which to simulate spatio-temporal observations. Each row
#' represents a point in space.
#' @param L An integer specifying the number of realizations (simulations) to generate. Default is 100.
#' @return An array of simulated spatio-temporal data with the following structure:
#' \item{Ysim}{An array of dimensions \code{(L, nrow(D), nt)}, where:
#' \describe{
#'   \item{\code{L}}{The number of realizations (simulations) generated.}
#'   \item{\code{nrow(D)}}{The number of spatial locations for the new data points.}
#'   \item{\code{nt}}{The number of temporal points.}
#' }}
#' @import DiceKriging mvtnorm
#' @export
#' @examples
#' # Example usage:
#' # Test function
#' fxt <- function(x, t) {sin(2 * pi * t * (x - t))}
#' # Design of experiments
#' D <- matrix(seq(0, 1, length = 5), ncol = 1)
#' # Discretized time points
#' tt <- seq(0, 1, length = 15)
#' # Spatio-temporal observations matrix
#' FD <- apply(array(tt), 1, function(t) {apply(D, 1, function(x) {fxt(x, t)})})
#'
#' # SVD decomposition and GP modelling
#' model <- svdgppmodel(D, FD, formula = ~1, K = 4, spcovtype = "matern5_2",  typekrig = "SK")
#' D <- matrix(runif(20), ncol=1)  # New spatial locations
#' Ysim <- simsvdgp(model, D, L=100)
simsvdgp <- function(model, D, L){
  if(!is.list(model)) stop("model must be a valid KLGP model list.")
  if(!is.matrix(D)) stop("D must be a matrix representing new spatial locations.")
  if(!is.numeric(L) || L <= 0) stop("L must be a positive integer representing the number of simulations.")

  GP <- model$GP
  Sigma <- model$Sigma
  V <- model$V
  K <- model$K
  type <- model$type
  eps=1e-10
  Sigmat <- V %*% (t(Sigma) %*% Sigma) %*% t(V)
  US <- array(data = 0, dim=c(L, nrow(D), K))
  for(k in 1:K){
    Yp <- DiceKriging::predict(GP[[k]], newdata = data.frame(D=D), type = type, se.compute = TRUE, cov.compute = TRUE, checkNames = FALSE)
    US[, , k] <- mvtnorm::rmvnorm(L, mean = Yp$mean,
                                  sigma = Yp$cov+diag(eps, nrow(Yp$cov)),
                                  method = "chol")
  }

  Ysim <- array(data = 0, dim=c(L, nrow(D), nrow(Sigmat)))
  for(l in 1:L){
    Ysim[l,,] <- model$FDmean + US[l, ,] %*% t(V)
  }

  return(Ysim)
}


#' Update an Existing SVDGP Model with New Data
#'
#' This function updates an existing KLGP model with new design points (spatial data) and corresponding spatio-temporal observations.
#'
#' @param model A trained SVDGP model, typically created using the \code{\link{svdgpmodel}} function.
#' The model contains Gaussian Process (GP) components and the Karhunen-Loève decomposition matrices.
#' @param Dnew A matrix representing the new design points (spatial data).
#' @param FDnew A matrix of new spatio-temporal observations.
#' The rows should correspond to the rows in \code{Dnew}, and the columns represent temporal data.
#' @param ... Additional arguments passed to the \code{\link[DiceKriging]{update}} function
#' in the DiceKriging package, such as \code{cov.reestim} or \code{trend.reestim}.
#' @return A list containing the updated SVDGP model, with the Gaussian Process components retrained
#' on the original and new data.
#' \item{model}{The updated SVDGP model, with retrained Gaussian Process components for each spatial mode.}
#' @importFrom DiceKriging update
#' @seealso \code{\link{svdgpmodel}}, \code{\link[DiceKriging]{update}}
#' @examples
#' # Example usage:
#' # Test function
#' fxt <- function(x, t) {sin(2 * pi * t * (x - t))}
#' # Design of experiments
#' D <- matrix(seq(0, 1, length = 50), ncol = 1)
#' # Discretized time points
#' tt <- seq(0, 1, length = 15)
#' # Spatio-temporal observations matrix
#' FD <- apply(array(tt), 1, function(t) {apply(D, 1, function(x) {fxt(x, t)})})
#'
#' # SVD and GP modelling
#' model <- svdgppmodel(D, FD, formula = ~1, K = 4, spcovtype = "matern5_2",  typekrig = "SK")
#' Dnew <- matrix(runif(10), ncol = 1)
#' FDnew <- apply(array(tt), 1, function(t) {apply(Dnew, 1, function(x) {fxt(x, t)})})
#' model <- updatesvdgp(model, Dnew, FDnew)
#' @export
updatesvdgp <- function(model, Dnew, FDnew, ...){

  if(!is.matrix(Dnew)) stop("Dnew must be a matrix")
  if(!is.matrix(FDnew)) stop("FDnew must be a matrix")
  if(nrow(FDnew)!=nrow(Dnew)) stop("FDnew must have the same number of row as Dnew!")

  FDm = mean(FDnew)
  FDc = FDnew - FDm
  K = model$K
  Ynew = FDc%*%model$V
  for (k in 1:K) {
    model$GP[[k]] <- DiceKriging::update(model$GP[[k]], newX=data.frame(D=Dnew), newy=data.frame(Ynew[,k]),
                                         cov.reestim = TRUE, trend.reestim = TRUE,...)
  }

  model$FDmean = (nrow(model$U)*model$FDmean + nrow(FDnew)*FDm)/(nrow(model$U)+nrow(FDnew))

  return(model)
}

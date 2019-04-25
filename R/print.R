
#' Summary method for ICCier objects.
#'
#' Summaries of parameter estimates and MCMC convergence diagnostics.
#'
#' @param object An ICCier object.
#' @param prob Probability interval (default: .95)
#' @param ... Currently not used.
#'
#' @return TBD
#' @export
#'
summary.ICCier <- function(object,prob=.95,...){

}

.get_beta <- function(object,prob=.95,...){

}

.get_gamma <- function(object,prob=.95,...){

}

.get_eta <- function(object,prob=.95,...){

}

.get_icc <- function(object,prob=.95,...){

}

.posterior_mean <- function(object,pars){
  samps <- as.matrix(object$fit, pars)
  colMeans(samps)
}

.posterior_sd <- function(object,pars){
  samps <- as.matrix(object$fit,pars)
  apply(samps,2,sd)
}

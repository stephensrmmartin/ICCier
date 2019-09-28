
#' ICCier Method for posterior intervals.
#'
#' @param object ICCier object.
#' @param prob A number p (0 < p < 1) indicating the desired probability mass to include in the intervals.
#' @param pars Parameter name.
#' @param ... Not used.
#'
#' @return Matrix with two columns (lower and upper) and one parameter per row.
#'
#' @keywords internal
posterior_interval.ICCier <- function(object, prob = .95, pars, ...){
  L <- (1-prob)/2
  U <- 1 - L
  samps <- as.matrix(object$fit,pars=pars)
  out <- t(apply(samps,2,function(x){quantile(x,probs=c(L,U))}))
  colnames(out) <- paste0(c(L,U)*100,'%')
  out
}

#' Extract pointwise log-likelihood
#'
#' @param object ICCier object
#' @param ... Not used.
#'
#' @return S by N matrix of log likelihoods.
#' @export
#' @importFrom rstantools log_lik
#' @export log_lik
#' @method log_lik ICCier
#' @aliases log_lik
#'
#' @keywords internal
log_lik.ICCier <- function(object,...){
  samps <- extract_log_lik(object$fit,...)
  samps
}

#' Extract the number of posterior samples.
#'
#' @param object ICCier object.
#' @param ... Not used.
#'
#' @return Real value. Number of posterior samples stored (post-warmup).
#' @export
#' @importFrom rstantools nsamples
#' @export nsamples
#' @method nsamples ICCier
#' @aliases nsamples
#' @keywords internal
nsamples.ICCier <- function(object, ...){
  n.samps <- (object$fit@sim$iter - object$fit@sim$warmup)*object$fit@sim$chains
  return(n.samps)
}

#' Compute leave-one-out (LOO).
#'
#' @param x ICCier object.
#'
#' @return Loo object. See \code{\link[loo]{loo}}.
#' @export
#' @import loo
#' @export loo
#' @method loo ICCier
#' @aliases loo
loo.ICCier <- function(x,...){
  object <- x
  LL_array <- log_lik(object,merge_chains=FALSE)
  r_eff <- relative_eff(exp(LL_array))
  loo(x=LL_array,r_eff=r_eff,...)
}

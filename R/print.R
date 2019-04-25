
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

print.ICCier <- function(object,...){

}

.get_beta <- function(object,prob=.95,...){
  mu <- .posterior_mean(object,'beta0')
  mu_group <- .posterior_mean(object, 'mu_group')
  mu.ci <- posterior_interval(object,prob=prob,pars='beta0')
  mu_group.ci <- posterior_interval(object,prob=prob,pars='mu_group')

  K <- object$stan_data$K
  mu_group <- matrix(mu_group,nrow=K)
  mu.L <- mu.ci[,1]
  mu.U <- mu.ci[,2]
  mu_group.L <- matrix(mu_group.ci[,1],nrow=K)
  mu_group.U <- matrix(mu_group.ci[,2],nrow=K)

  return(mget(c('mu','mu_group','mu.L','mu.U','mu_group.L','mu_group.U')))
}

.get_gamma <- function(object,prob=.95,...){
  gamma <- .posterior_mean(object, 'gamma')
  gamma_group <- .posterior_mean(object, 'gamma_group')
  gamma.ci <- posterior_interval(object, prob=prob, pars = 'gamma')
  gamma_group.ci <- posterior_interval(object, prob=prob, pars = 'gamma_group')

  # Format as matrix
  K <- object$stan_data$K
  P_l1 <- object$stan_data$P_l1
  P_l2 <- object$stan_data$P_l2
  gamma <- matrix(gamma,nrow=P_l2)
  gamma_group <- matrix(gamma_group,nrow=K)
  gamma.L <- matrix(gamma.ci[,1],nrow=P_l2)
  gamma.U <- matrix(gamma.ci[,2],nrow=P_l2)
  gamma_group.L <- matrix(gamma_group.ci[,1],nrow=K)
  gamma_group.U <- matrix(gamma_group.ci[,2],nrow=K)

  return(mget(c('gamma','gamma_group','gamma.L','gamma.U','gamma_group.L','gamma_group.U')))
}

.get_eta <- function(object,prob=.95,...){
  eta <- .posterior_mean(object, 'eta')
  eta.ci <- posterior_interval(object,prob=prob,pars = 'eta')

  P_l1 <- object$stan_data$P_l1
  P_l2 <- object$stan_data$P_l2
  eta <- matrix(eta, nrow = P_l2)
  eta.L <- matrix(eta.ci[,1],nrow=P_l2)
  eta.U <- matrix(eta.ci[,2],nrow=P_l2)

  return(mget(c('eta','eta.L','eta.U')))
}

.posterior_mean <- function(object,pars){
  samps <- as.matrix(object$fit, pars)
  colMeans(samps)
}

.posterior_sd <- function(object,pars){
  samps <- as.matrix(object$fit,pars)
  apply(samps,2,sd)
}


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
  beta <- .get_beta(object,prob,...)
  gamma <- .get_gamma(object,prob,...)
  eta <- .get_eta(object,prob,...)
  prob <- prob
  formula <- object$formula

  out <- mget(c('formula','prob','beta','gamma','eta'))
  class(out) <- 'summary.ICCier'
  out
}

print.ICCier <- function(object,...){
  cat('Formula:',deparse(object$formula),'\n')
  cat('\n')


  cat('Coefficients:','\n')
  cat('Mean: \t',format(.get_beta(object)$mu,...),'\n\n')
  cat('L1 Scale: \n'); print(t(.get_gamma(object)$gamma),...); cat('\n')
  cat('L2 Scale: \n'); print(t(.get_eta(object)$eta),...); cat('\n')

  invisible(object)
}

print.summary.ICCier <- function(object){
  cat('Formula:',deparse(object$formula))

  invisible(object)
}

.diag_rhat <- function(object){

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

  out <- mget(c('mu','mu.L','mu.U','mu_group','mu_group.L','mu_group.U'))

  return(out)
}

.get_gamma <- function(object,prob=.95,...){
  fnames <- .get_formula_names(object)

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

  rownames(gamma) <- rownames(gamma.L) <- rownames(gamma.U) <- fnames$l2
  colnames(gamma) <- colnames(gamma.L) <- colnames(gamma.U) <- fnames$l1

  colnames(gamma_group) <- colnames(gamma_group.L) <- colnames(gamma_group.U) <- fnames$l1

  out <- mget(c('gamma','gamma.L','gamma.U','gamma_group','gamma_group.L','gamma_group.U'))

  return(out)
}

.get_eta <- function(object,prob=.95,...){
  fnames <- .get_formula_names(object)
  eta <- .posterior_mean(object, 'eta')
  eta.ci <- posterior_interval(object,prob=prob,pars = 'eta')

  P_l1 <- object$stan_data$P_l1
  P_l2 <- object$stan_data$P_l2
  eta <- matrix(eta, nrow = P_l2)
  eta.L <- matrix(eta.ci[,1],nrow=P_l2)
  eta.U <- matrix(eta.ci[,2],nrow=P_l2)

  rownames(eta) <- rownames(eta.L) <- rownames(eta.U) <- fnames$l2
  colnames(eta) <- colnames(eta.L) <- colnames(eta.U) <- c('(Intercept.L)',fnames$l1)

  out <- mget(c('eta','eta.L','eta.U'))

  return(out)
}

.posterior_mean <- function(object,pars){
  samps <- as.matrix(object$fit, pars)
  colMeans(samps)
}

.posterior_sd <- function(object,pars){
  samps <- as.matrix(object$fit,pars)
  apply(samps,2,sd)
}

.get_formula_names <- function(object){
  l1 <- colnames(object$stan_data$x_sca_l1)
  l2 <- colnames(object$stan_data$x_sca_l2)
  outcome <- colnames(model.part(object$formula,object$data,lhs=1))
  grouping <- colnames(model.part(object$formula,object$data,lhs=2))
  return(mget(c('l1','l2','outcome','grouping')))
}

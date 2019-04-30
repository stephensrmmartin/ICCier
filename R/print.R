
#' Summary method for ICCier objects.
#'
#' Summaries of parameter estimates and MCMC convergence diagnostics.
#'
#' @param object An ICCier object.
#' @param prob Probability interval (default: .95)
#' @param ... Currently not used.
#'
#' @return ICCier summary object. List containing:
#' \describe{
#' \item{formula}{Model formula}
#' \item{prob}{The specified prob argument}
#' \item{estimate}{List of matrices of mean posterior point estimates}
#' \item{ci}{List containing the lower (L) and upper (U) interval estimates for each parameter matrix.}
#' \item{meta}{Meta-data including diagnostics, number of iterations, number of chains, total number of observations, and number of groups (persons)}
#' }
#' @export
#'
summary.ICCier <- function(object,prob=.95,...){
  beta <- .get_beta(object,prob,...)
  gamma <- .get_gamma(object,prob,...)
  eta <- .get_eta(object,prob,...)
  omega <- .get_omega(object,prob,...)
  prob <- prob
  formula <- object$formula

  chains <- object$fit@sim$chains
  iter <- list(iter=object$fit@sim$iter)
  iter$post <- iter$iter - object$fit@sim$warmup
  iter$total <- iter$post*chains
  K <- object$stan_data$K
  N <- object$stan_data$N
  meta <- list(diagnostics=object$diagnostics,iter=iter,chains=chains,N=N,K=K)

  estimate <- list(beta = beta$mu,
                gamma = t(gamma$gamma),
                eta = t(eta$eta),
                omega=omega$omega,
                beta_group = beta$beta_group,
                gamma_group = gamma$gamma_group
                )
  ci.L <- list(beta = beta$mu.L,
                gamma = t(gamma$gamma.L),
                eta = t(eta$eta.L),
                omega=omega$omega.L,
                beta_group = beta$beta_group.L,
                gamma_group = gamma$gamma_group.L
                )
  ci.U <- list(beta = beta$mu.U,
                gamma = t(gamma$gamma.U),
                eta = t(eta$eta.U),
                omega=omega$omega.U,
                beta_group = beta$beta_group.U,
                gamma_group = gamma$gamma_group.U
                )
  ci <- list(L=ci.L, U=ci.U)

  out <- list(formula=formula,
              prob=prob,
              estimate=estimate,
              ci=ci,
              meta=meta
              )
  class(out) <- 'summary.ICCier'
  out
}

#' Print method for ICCIer object.
#'
#' @param object ICCier object.
#' @param ... Further arguments to \code{print}.
#'
#' @export
#'
print.ICCier <- function(object,...){
  cat('Formula:',deparse(object$formula),'\n')
  cat('\n')


  cat('Coefficients:','\n')
  cat('Mean: \t',format(.get_beta(object)$mu,...),'\n\n')
  cat('L1 Scale: \n'); print(t(.get_gamma(object)$gamma),...); cat('\n')
  cat('L2 Scale: \n'); print(t(.get_eta(object)$eta),...); cat('\n')
  cat('L1 Cor: \n'); print(.get_omega(object)$omega,...); cat('\n')

  invisible(object)
}

#' Print method for ICCier summaries.
#'
#' @param object Output of \code{summary(ICCierObject)}.
#' @inheritParams print.ICCier
#'
#' @return Invisibly returns summary object.
#' @export
#'
print.summary.ICCier <- function(object,...){

  cat('Formula:',deparse(object$formula),'\n')
  cat('Number of observations:', object$meta$N,'\n')
  cat('Number of groups:',object$meta$K,'\n')
  cat('\n--------------------\n')
  .print_diagnostics(object$meta$diagnostics)
  cat('\n--------------------\n')

  # colnames(object$ci$L$eta) <- rep(paste0((1-object$prob)/2*100,'%'),ncol(object$ci$L$eta))
  # colnames(object$ci$U$eta) <- rep(paste0((1 - (1-object$prob)/2)*100,'%'),ncol(object$ci$U$eta))
  # colnames(object$ci$L$gamma) <- rep(paste0((1-object$prob)/2*100,'%'),ncol(object$ci$L$gamma))
  # colnames(object$ci$U$gamma) <- rep(paste0((1 - (1-object$prob)/2)*100,'%'),ncol(object$ci$U$gamma))
  # colnames(object$ci$L$omega) <- rep(paste0((1-object$prob)/2*100,'%'),ncol(object$ci$L$omega))
  # colnames(object$ci$U$omega) <- rep(paste0((1 - (1-object$prob)/2)*100,'%'),ncol(object$ci$U$omega))

  beta.sum <- unname(c(format(object$estimate$beta,...),paste0('[',format(object$ci$L$beta,...),' ',format(object$ci$U$beta,...),']')))
  eta.sum <- cbind(format(object$estimate$eta,...),matrix(paste0('[',format(object$ci$L$eta,...),' ',format(object$ci$U$eta,...),']'),nrow=nrow(object$estimate$eta)))
  gamma.sum <- cbind(format(object$estimate$gamma,...),matrix(paste0('[',format(object$ci$L$gamma,...),' ',format(object$ci$U$gamma,...),']'),nrow=nrow(object$estimate$gamma)))
  names(beta.sum) <- c('',paste0(object$prob*100,'%'))
  colnames(eta.sum)[colnames(eta.sum) == ''] <- paste0(object$prob*100,'%')
  colnames(gamma.sum)[colnames(gamma.sum) == ''] <- paste0(object$prob*100,'%')

  cat('Coefficients:','\n\n')
  cat('Mean: \n');print(beta.sum,quote=FALSE,...);cat('\n')
  cat('L1 Scale: \n'); print(gamma.sum,quote=FALSE); cat('\n')
  cat('L2 Scale: \n'); print(eta.sum,quote=FALSE); cat('\n')
  cat('L1 Cor: \n'); print(object$estimate$omega,...); cat('\n')

  invisible(object)
}
## TODO: Add digits to summary() as undocumented parameter. Store in meta. Pull into print, so can just run summary(object,digits=2).
## TODO: Add option to extract gamma_group and beta_group. Not in summary, but in coef().

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

.get_omega <- function(object,prob=.95,...){
  fnames <- .get_formula_names(object)
  omega <- .posterior_mean(object,'Omega')
  omega.ci <- posterior_interval(object,prob=prob,pars='Omega')

  D <- object$stan_data$P_l1 + 1
  omega <- matrix(omega,nrow=D)
  omega.L <- matrix(omega.ci[,1],D)
  omega.U <- matrix(omega.ci[,2],D)

  rownames(omega) <- colnames(omega) <- rownames(omega.L) <- colnames(omega.L) <-rownames(omega.U) <- colnames(omega.U) <- c('(Intercept.L)',fnames$l1)

  out <- mget(c('omega','omega.L','omega.U'))

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

.print_diagnostics <- function(diagnostics){
  diags <- diagnostics
  cat('Diagnostics:','\n')
  if(sum(diags$rhats > 1.1)){
    cat('\t Rhats: Failed\n')
    cat('\t Some Rhats > 1.1. Do not interpret results! The largest 10 are:\n')
    print(head(sort(diags$rhats,decreasing=TRUE),10))
  } else{
    cat('\t Rhats: Passed\n')
  }

  if(diags$div > 0){
    cat('\t Divergent Transitions: Failed -',diags$div,'divergent transitions detected. \n')
  } else {
    cat('\t Divergent transitions: Passed \n')
  }

  if(diags$tree.max > 0){
    cat('\t Max treedepth hit:',diags$tree.max,'\n')
  }
  if(any(diags$bfmi < .2)){
    cat('\t Low E-BFMI detected in chains',which(diags$bfmi < .2), '\n')
  }

}

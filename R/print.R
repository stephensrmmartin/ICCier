
#' Summary method for ICCier objects.
#'
#' Summaries of parameter estimates and MCMC convergence diagnostics.
#'
#' @param object An ICCier object.
#' @inheritParams posterior_interval.ICCier
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
  dots <- list(...)
  if(is.null(dots$digits)){
    digits <- 3
  } else {
    digits <- dots$digits
    dots$digits <- NULL
  }
  beta <- .get_beta(object,prob,...)
  gamma <- .get_gamma(object,prob,...)
  eta <- .get_eta(object,prob,...)
  omega <- .get_omega(object,prob,...)
  icc_mean <- .get_icc_mean(object,prob,...)
  icc_sd <- .get_icc_sd(object,prob,...)
  prob <- prob
  formula <- object$formula
  type <- .get_model_type(object)

  chains <- object$fit@sim$chains
  iter <- list(iter=object$fit@sim$iter)
  iter$post <- iter$iter - object$fit@sim$warmup
  iter$total <- iter$post*chains
  K <- object$stan_data$K
  N <- object$stan_data$N
  meta <- list(diagnostics=object$diagnostics,iter=iter,chains=chains,N=N,K=K,digits=digits,type=type)

  estimate <- list(beta = t(beta$mu),
                gamma = t(gamma$gamma),
                eta = t(eta$eta),
                omega=omega$omega,
                beta_group = beta$mu_group,
                gamma_group = gamma$gamma_group,
                icc_mean = icc_mean$icc_mean,
                icc_sd = icc_sd$icc_sd
                )
  ci.L <- list(beta = t(beta$mu.L),
                gamma = t(gamma$gamma.L),
                eta = t(eta$eta.L),
                omega=omega$omega.L,
                beta_group = beta$mu_group.L,
                gamma_group = gamma$gamma_group.L,
                icc_mean = icc_mean$icc_mean.L,
                icc_sd = icc_sd$icc_sd.L
                )
  ci.U <- list(beta = t(beta$mu.U),
                gamma = t(gamma$gamma.U),
                eta = t(eta$eta.U),
                omega=omega$omega.U,
                beta_group = beta$mu_group.U,
                gamma_group = gamma$gamma_group.U,
                icc_mean = icc_mean$icc_mean.U,
                icc_sd = icc_sd$icc_sd.U
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
#' @param x ICCier object.
#' @param ... Further arguments to \code{print}.
#'
#' @export
#' @keywords internal
print.ICCier <- function(x,...){
  object <- x
  cat('Formula:',deparse(object$formula),'\n')
  cat('Type: ',.get_model_type(object),'\n')
  cat('\n')


  cat('Coefficients:','\n')
  cat('Mean Model: \n');print(t(.get_beta(object)$mu),...); cat('\n\n')
  cat('Within-group (log) SD: \n'); print(t(.get_gamma(object)$gamma),...); cat('\n')
  cat('Between-group (log) SD: \n'); print(t(.get_eta(object)$eta),...); cat('\n')
  cat('Random Effect Correlations: \n'); print(.get_omega(object)$omega,...); cat('\n')

  invisible(object)
}

#' Print method for ICCier summaries.
#'
#' @param x Output of \code{summary(ICCierObject)}.
#' @param matrix Logical. Print coefficients as table (FALSE) or as HLM-like matrix (TRUE).
#' @inheritParams print.ICCier
#'
#' @return Invisibly returns summary object.
#' @export
#' @keywords internal
print.summary.ICCier <- function(x,matrix=FALSE,...){
  object <- x

  dots <- list(...)
  digits <- dots$digits
  if(is.null(digits)){
    digits <- object$meta$digits
  }

  # Header
  cat('Formula:',deparse(object$formula),'\n')
  cat('Type: ',object$meta$type,'\n')
  cat('Number of observations:', object$meta$N,'\n')
  cat('Number of groups:',object$meta$K,'\n')
  cat('\n--------------------\n')
  .print_diagnostics(object$meta$diagnostics)
  cat('\n--------------------\n')

  # Formatting
  if(matrix){
    beta.sum <- .print_matrix(object,digits,'beta',...)
    gamma.sum <- .print_matrix(object,digits,'gamma',...)
    eta.sum <- .print_matrix(object,digits,'eta',...)
    icc.sum <- .print_matrix(object,digits,'icc',...)
  } else {
    beta.sum <- .print_table(object,digits,'beta',...)
    gamma.sum <- .print_table(object,digits,'gamma',...)
    eta.sum <- .print_table(object,digits,'eta',...)
    icc.sum <- .print_table(object,digits,'icc',...)
  }

  # Print
  cat('ICC Summary:','\n'); print(icc.sum,quote=FALSE)
  cat('\n--------------------\n')
  cat('Coefficients:','\n\n')
  cat('Mean Model: \n');print(beta.sum,quote=FALSE,...);cat('\n')
  cat('Within-group (log) SD: \n'); print(gamma.sum,quote=FALSE); cat('\n')
  cat('Between-group (log) SD: \n'); print(eta.sum,quote=FALSE); cat('\n')
  cat('Random Effect Correlations: \n'); print(round(object$estimate$omega,digits),...); cat('\n')

  invisible(object)
}

#' Generate table for summary
#'
#' Creates coefficient matrix ala lm and others.
#' @inheritParams .print_matrix
#' @keywords internal
.print_table <- function(object, digits, param, ...){
  prob <- object$prob
  probs <- c((1-prob)/2,1-(1-prob)/2)

  if(param == 'icc'){
    est <- c(object$estimate[['icc_mean']],object$estimate[['icc_sd']])
    ci.L <- c(object$ci$L[['icc_mean']],object$ci$L[['icc_sd']])
    ci.U <- c(object$ci$U[['icc_mean']],object$ci$U[['icc_sd']])
    out <- round(cbind(est,ci.L,ci.U),digits)
    rownames(out) <- c('Mean','SD')
    colnames(out) <- c('Estimate',paste0(probs*100,'%'))
    return(out)
  }
  # Extract
  est <- object$estimate[[param]]
  ci.L <- object$ci$L[[param]]
  ci.U <- object$ci$U[[param]]

  # Unroll
  names.l1 <- rownames(est)
  names.l2 <- colnames(est)
  nameMat <- t(sapply(names.l1,function(n){
                paste0(n,':',names.l2)
              }))
  names.vec <- as.vector(nameMat)
  names.vec <- gsub('(Intercept):(Intercept)','(Intercept)',names.vec,fixed=TRUE)
  if(param != 'eta'){
    names.vec <- gsub(':\\(Intercept\\)|\\(Intercept\\):','',names.vec)
  }

  out <- cbind(as.vector(est),as.vector(ci.L),as.vector(ci.U))
  out <- round(out,digits)
  colnames(out) <- c('Estimate',paste0(probs*100,'%'))
  rownames(out) <- names.vec
  if(param == 'eta'){
    out <- out[grepl('^Mean_',names.vec),,drop=FALSE]
  }
  return(out)

}

#' Generate character matrix for summary
#'
#' @param object summary.ICCier object
#' @param digits digits
#' @param param param (beta,gamma,eta,icc)
#' @param ... Arguments passed to format.
#'
#' @keywords internal
.print_matrix <- function(object, digits, param,...){
  # Extract
  if(param == 'icc'){
    est <-  matrix(unlist(object$estimate[c('icc_mean','icc_sd')],use.names = FALSE),ncol=1)
    ci.L <- matrix(unlist(object$ci$L[c('icc_mean','icc_sd')],use.names = FALSE),ncol=1)
    ci.U <- matrix(unlist(object$ci$U[c('icc_mean','icc_sd')],use.names = FALSE),ncol=1)
  } else {
    est <- object$estimate[[param]]
    ci.L <- object$ci$L[[param]]
    ci.U <- object$ci$U[[param]]
  }

  # Format
  est <-  format(round(est,digits),...)
  ci.L <- format(round(ci.L,digits),...)
  ci.U <- format(round(ci.U,digits),...)

  # Add brackets
  ci <- paste0('[',ci.L,' ',ci.U,']')
  ci <- matrix(ci,nrow=nrow(est))

  # Combine and rename
  out <- cbind(est, ci)
  if(param == 'icc'){
    colnames(out) <- c('',paste0(object$prob*100,'%'))
    rownames(out) <- c('Mean','SD')
  } else{
    colnames(out)[colnames(out) == ''] <- paste0(object$prob*100,'%')
  }
  names(dimnames(out)) <- names(dimnames(est))
  return(out)
}

.get_model_type <- function(object){
  type <- object$type
  typeString <- ifelse(type$conditional,paste0('Conditional (',ifelse(type$adjusted,'Adjusted','Unadjusted'),')'),'Unconditional')
  return(typeString)
}

.get_beta <- function(object,prob=.95,...){
  fnames <- .get_formula_names(object)

  mu <- .posterior_mean(object,'beta0')
  mu_group <- .posterior_mean(object, 'mu_group')
  mu.ci <- posterior_interval(object,prob=prob,pars='beta0')
  mu_group.ci <- posterior_interval(object,prob=prob,pars='mu_group')

  # Format as matrix
  K <- object$stan_data$K
  Q_l1 <- object$stan_data$Q_l1
  Q_l2 <- object$stan_data$Q_l2

  mu <- matrix(mu, Q_l2, Q_l1)
  mu_group <- matrix(mu_group,nrow=K)
  mu.L <- matrix(mu.ci[,1], Q_l2, Q_l1)
  mu.U <- matrix(mu.ci[,2], Q_l2, Q_l1)
  mu_group.L <- matrix(mu_group.ci[,1],nrow=K)
  mu_group.U <- matrix(mu_group.ci[,2],nrow=K)

  rownames(mu) <- rownames(mu.L) <- rownames(mu.U) <- fnames$l2.loc
  colnames(mu) <- colnames(mu.L) <- colnames(mu.U) <- fnames$l1.loc
  names(dimnames(mu)) <- names(dimnames(mu.L)) <- names(dimnames(mu.U)) <- c('Level 2','Level 1')

  colnames(mu_group) <- colnames(mu_group.L) <- colnames(mu_group.U) <- fnames$l1.loc

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
  names(dimnames(gamma)) <- names(dimnames(gamma.L)) <- names(dimnames(gamma.U)) <- c('Level 2','Level 1')

  colnames(gamma_group) <- colnames(gamma_group.L) <- colnames(gamma_group.U) <- fnames$l1

  out <- mget(c('gamma','gamma.L','gamma.U','gamma_group','gamma_group.L','gamma_group.U'))

  return(out)
}

.get_eta <- function(object,prob=.95,...){
  fnames <- .get_formula_names(object, prefix=TRUE)
  eta <- .posterior_mean(object, 'eta')
  eta.ci <- posterior_interval(object,prob=prob,pars = 'eta')

  P_l1 <- object$stan_data$P_l1
  P_l2 <- object$stan_data$P_l2
  R_l2 <- object$stan_data$R_l2
  Q_l1 <- object$stan_data$Q_l1
  eta <- matrix(eta, nrow = R_l2,ncol=P_l1 + Q_l1)
  eta.L <- matrix(eta.ci[,1],nrow=R_l2)
  eta.U <- matrix(eta.ci[,2],nrow=R_l2)

  rownames(eta) <- rownames(eta.L) <- rownames(eta.U) <- fnames$bet
  colnames(eta) <- colnames(eta.L) <- colnames(eta.U) <- c(fnames$l1.loc,fnames$l1)
  names(dimnames(eta)) <- names(dimnames(eta.L)) <- names(dimnames(eta.U)) <- c('Bet. Group','')

  out <- mget(c('eta','eta.L','eta.U'))

  return(out)
}

.get_omega <- function(object,prob=.95,...){
  fnames <- .get_formula_names(object,prefix=TRUE)
  omega <- .posterior_mean(object,'Omega')
  omega.ci <- posterior_interval(object,prob=prob,pars='Omega')

  D <- object$stan_data$P_l1 + object$stan_data$Q_l1
  omega <- matrix(omega,nrow=D)
  omega.L <- matrix(omega.ci[,1],D)
  omega.U <- matrix(omega.ci[,2],D)

  rownames(omega) <- colnames(omega) <- rownames(omega.L) <- colnames(omega.L) <-rownames(omega.U) <- colnames(omega.U) <- c(fnames$l1.loc,fnames$l1)

  out <- mget(c('omega','omega.L','omega.U'))

  return(out)

}

.get_icc_mean <- function(object,prob=.95,...){
  icc_mean <- .posterior_mean(object,'icc_mean')
  icc_mean.ci <- posterior_interval(object,prob=prob,pars='icc_mean')
  icc_mean.L <- icc_mean.ci[1]
  icc_mean.U <- icc_mean.ci[2]
  out <- mget(c('icc_mean','icc_mean.L','icc_mean.U'))
  return(out)

}
.get_icc_sd <- function(object, prob=.95,...){
  icc_sd <- .posterior_mean(object,'icc_sd')
  icc_sd.ci <- posterior_interval(object,prob=prob,pars='icc_sd')
  icc_sd.L <- icc_sd.ci[1]
  icc_sd.U <- icc_sd.ci[2]
  out <- mget(c('icc_sd','icc_sd.L','icc_sd.U'))
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

.get_formula_names <- function(object, prefix=FALSE){
  l1 <- colnames(object$stan_data$x_sca_l1)
  l2 <- colnames(object$stan_data$x_sca_l2)
  l1.loc <- colnames(object$stan_data$x_loc_l1)
  l2.loc <- colnames(object$stan_data$x_loc_l2)
  bet  <- colnames(object$stan_data$x_bet_l2)
  outcome <- colnames(model.part(object$formula,object$data,lhs=1))
  grouping <- colnames(model.part(object$formula,object$data,lhs=2))
  if(prefix){
    l1.loc <- paste0('Mean_',l1.loc)
    l2.loc <- paste0('Mean_',l2.loc)
  }
  return(mget(c('l1','l2','outcome','grouping','l1.loc','l2.loc','bet')))
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

#' ICCier method to predict new ICC values.
#'
#' @param object ICCier object
#' @param newdata Data to predict from. If NULL, calls \code{\link{fitted.ICCier}} on fit data.
#' @param draws Number of draws to use from posterior samples for prediction. Default: All of them.
#' @inheritParams fitted.ICCier
#' @param ... Not currently used.
#'
#' @return TBD
#' @export
#'
predict.ICCier <- function(object, newdata=NULL, draws=NULL,summary=TRUE,prob=.95,inc_group=TRUE, ...){
  fnames <- .get_formula_names(object)
  if(is.null(data)){
   return(fitted(object,summary,prob,inc_group))
  }
  total_iter <- (object$fit@sim$iter - object$fit@sim$warmup)*object$fit@sim$chains
  if(is.null(draws)){
    draws <- total_iter
  }
  if(draws > total_iter){
    stop('More draws specified than actually exist.')
  }

  if(!(fnames$grouping %in% colnames(newdata))){
    grouping_available <- FALSE
    message('No grouping variable. Using fixef only.')
    newdata[,fnames$grouping] <- NA
  } else{
    grouping_available <- TRUE
  }
  group_known <- newdata[,fnames$grouping] %in% object$data[,fnames$grouping]

  samps <- .extract_transform(object,draws)

}
# TODO: Needs to handle existing (known) groups as well as unknown.
# If known, pull from the group_map which integer they belong to.
# If unknown, assign a value, then randomly draw from RE distribution(s) instead on each iteration.
# If no groups specified, just use fixed effects, b/c there's no other information available.
# Remember: You only need to predict ICC = var(mean)/(var(mean) + var(within)) for each row, across samples.
# TODO: Actually, you need random_z for known groups. random includes the effect of eta already; so to ...
# predict new values, you need var(within) ~~ x_sca_l1%*%t(predicted_gamma + u_gammas), and u_gammas ...
# ~ MVN(0, Sigma(eta,x_sca_l2)). B/c x_sca_l2 can differ, we need _z separated out, predict new RE SDs, then ...
# create new Sigma.

#' Extracts random effect samples.
#'
#' gamma_group and beta_group are returned.
#' However, for prediction, we need the group-specific random effects.
#' This function takes the object and returns the computed REs.
#'
#' group_gamma = x_sca_l2%*%gamma + group_gamma_random
#'
#' @param object ICCier object.
#'
#' @keywords internal
#' @return Samples of RE matrices: Arrays with dimensions KxP_l1xdraws.
.get_random_effect_samples <- function(object,draws=NULL){
  if(is.null(draws)){
    draws <- (object$fit@sim$iter - object$fit@sim$warmup)*object$fit@sim$chains
  }
  K <- object$stan_data$K
  P_l1 <- object$stan_data$P_l1
  P_l2 <- object$stan_data$P_l2
  X <- object$stan_data$x_sca_l2

  samps.mu <- as.matrix(object$fit,pars='beta0')[1:draws,,drop=FALSE]
  samps.gamma <- array(t(as.matrix(object$fit,pars='gamma')[1:draws,]),dim=c(P_l2,P_l1,draws))
  samps.mu_group <- as.matrix(object$fit,pars='mu_group')[1:draws,,drop=FALSE]
  samps.mu_group <- array(t(samps.mu_group),dim=c(K,1,draws))
  samps.gamma_group <- as.matrix(object$fit,pars='gamma_group')[1:draws,,drop=FALSE]
  samps.gamma_group <- array(t(samps.gamma_group),dim=c(K,P_l1,draws))

  group_mu_random <- sapply(1:draws,FUN=function(x){
    samps.mu_group[,,x] - samps.mu[x,1]
  },simplify='array')
  group_mu_random <- array(group_mu_random,dim=c(K,1,draws))
  group_gamma_random <- sapply(1:draws,FUN=function(x){
    samps.gamma_group[,,x] - X%*%samps.gamma[,,x]
  },simplify = 'array')
  group_gamma_random <- array(group_gamma_random,dim=c(K,P_l1,draws))

  colnames(group_mu_random) <- '(Intercept.L)'
  colnames(group_gamma_random) <- colnames(X)

  list(mu_random = group_mu_random, gamma_random = group_gamma_random)
}
# TODO: Will need to compute the actual random effects, not just gamma_group. Can always take samps and subtract off rather than reworking everything else.

#' Extract standardized RE
#'
#' log(sds) = x_sca_l2%*%eta
#' RE = RE_z%\*%t(diag(sds)%\*%L_cor)
#' RE_z = RE %*% solve(t(diag(sds)%*%L_cor))
#'
#' @param object ICCier object
#' @inheritParams predict.ICCier
#'
#' @return Kx(P_l1 + 1)xS array of standardized REs.
#' @keywords internal
#'
.get_random_effect_z_samples <- function(object,draws=NULL){
  if(is.null(draws)){
    draws <- (object$fit@sim$iter - object$fit@sim$warmup)*object$fit@sim$chains
  }
  rand_samps <- .get_random_effect_samples(object,draws)
  samps <- sapply(1:draws,function(x){cbind(rand_samps$mu_random[,,x],rand_samps$gamma_random[,,x])},simplify='array')
  omega <- array(t(as.matrix(object$fit,pars='Omega')[1:draws,]),dim=c(object$stan_data$P_l1 + 1, object$stan_data$P_l1 + 1,draws))
  eta <- array(t(as.matrix(object$fit,pars='eta')[1:draws,]),dim=c(object$stan_data$P_l2, object$stan_data$P_l1 + 1, draws))

  L_omega <- array(apply(omega,3,function(x){t(chol(x))}),dim=dim(omega))
  sds <- sapply(1:draws,FUN = function(x){
    (exp(object$stan_data$x_sca_l2 %*% eta[,,x]))
  }, simplify='array')

  sapply(1:draws,FUN=function(x){
    t(sapply(1:dim(sds)[1], FUN=function(y){
      samps[y,,x] %*% solve(t(diag(sds[y,,x]) %*% L_omega[,,x]))
    },simplify='matrix')
    )
  },simplify='array')

}

#' Extracts samples, turns them into formed matrices for prediction
#'
#' Extracts \code{draws} of the generative parameters.
#' For matrices, converts to array; e.g., gamma is a [P_l2,P_l1,draws] array.
#' This will ease the computation of predictions, because it can be applied over draws.
#'
#' @param object ICCier object.
#' @param draws Number of draws
#'
#' @keywords internal
#'
#' @return List of arrays.
.extract_transform <- function(object,draws){
  K <- object$stan_data$K
  P_l2 <- object$stan_data$P_l2
  P_l1 <- object$stan_data$P_l1
  # Only these are needed for predicting ICC
  samps <- as.matrix(object$fit, pars = c('gamma','eta','gamma_group','Omega'))[1:draws,]
  gamma.cols <- grep('gamma.*',colnames(samps),value = TRUE)
  eta.cols <- grep('eta.*',colnames(samps),value = TRUE)
  gamma_group.cols <- grep('gamma_group.*',colnames(samps),value = TRUE)
  Omega.cols <- grep('Omega.*',colnames(samps),value = TRUE)

  gamma <- array(t(samps[,gamma.cols]),dim=c(P_l2,P_l1,draws))
  eta <- array(t(samps[,eta.cols]),dim=c(P_l2,P_l1 + 1,draws))
  Omega <- array(t(samps[,Omega.cols]),dim=c(P_l1 + 1,P_l1 + 1,draws))
  gamma_group <- array(t(samps[,gamma_group.cols]),dim=c(K,P_l1,draws))
  gamma_random <- .get_random_effect_samples(object,draws)$gamma_random
  return(mget(c('gamma','eta','gamma_random','Omega')))
}

#' ICCier method to extract ICC values.
#'
#' @param object ICCier object
#' @param summary Logical. Whether to return summary (mean, intervals) of ICCs (TRUE), or posterior samples (FALSE)
#' @inheritParams posterior_interval.ICCier
#' @param inc_group Logical. Whether to include the grouping variable with the estimates.
#' @param ... Not currently used.
#'
#' @return If \code{summary=TRUE}, then the mean and \code{prob}% intervals are returned for each observation in the model frame.
#' If \code{summary=FALSE}, an S by N matrix containing the S posterior samples for N observations.
#' @export
#'
fitted.ICCier <- function(object, summary=TRUE, prob=.95,inc_group=TRUE){
  if(summary) {
    out <- as.data.frame(cbind(mean=matrix(.posterior_mean(object,pars='icc'),ncol=1),posterior_interval(object,prob=prob,pars='icc')))
    colnames(out)[1] <- 'mean'
    if(inc_group){
      out[,.get_formula_names(object)$grouping] <- object$group_map[,1]
    }
  } else {
    out <- as.matrix(object$fit,pars='icc')
  }

  return(out)

}

#' ICCier method to predict new ICC values.
#'
#' @param object ICCier object
#' @param data Data to predict from.
#' @param ...
#'
#' @return TBD
#' @keywords internal
#'
posterior_predict.ICCier <- function(object, data, ...){
 ## TODO
}

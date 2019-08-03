#' ICCier method to predict new ICC values.
#'
#' Predicts new ICC values.
#'
#' If the grouping variable column is not included in \code{newdata}, then \code{predict} only predicts
#' from the fixed values.
#' This is \emph{not recommended}.
#' If the grouping variable column is present in \code{newdata}, but a value is either \code{NA} or not previously fit,
#' then \code{predict} marginalizes over the entire random effects distribution.
#' Finally, if the grouping variable is present and known, then \code{predict} marginalizes over the
#' uncertainty in that person's (standardized, orthogonalized) random effects.
#'
#' We recommend including the grouping column, even if the grouping identifier for the predicted case
#' is unknown. This ensures that uncertainty in the predicted ICC includes the uncertainty of that
#' case's random effects.
#'
#' @param object ICCier object.
#' @param newdata Data to predict from. If NULL, calls \code{\link{fitted.ICCier}} on fit data.
#' @param draws Number of draws to use from posterior samples for prediction. Default: All of them.
#' @inheritParams fitted.ICCier
#' @param ... Not currently used.
#' @importFrom mvtnorm rmvnorm
#'
#' @inherit fitted.ICCier return
#' @export
#'
predict.ICCier <- function(object, newdata=NULL, draws=NULL,summary=TRUE,prob=.95,inc_group=TRUE,occasion=NULL, ...){
  fnames <- .get_formula_names(object)
  magic_NA <- 'NA_ICCier'
  magic_ignore <- 'Ignore_ICCier'
  if(is.null(data)){
   return(fitted(object,summary,prob,inc_group))
  }
  total_iter <- nsamples(object)
  if(is.null(draws)){
    draws <- total_iter
  }
  if(draws > total_iter){
    stop('More draws specified than actually exist.')
  }

  if(!(fnames$grouping %in% colnames(newdata))){
    grouping_available <- FALSE
    message('No grouping variable. Using fixef only.')
    newdata[,fnames$grouping] <- magic_ignore
  } else{
    grouping_available <- TRUE
    group_known <- newdata[,fnames$grouping] %in% object$group_map$group_L2[,fnames$grouping]
    newdata[!group_known,fnames$grouping] <- magic_NA
  }

  samps <- .extract_transform(object,draws)
  dat <- .parse_formula(object$formula,newdata,predict=TRUE)
  out <- sapply(1:nrow(dat$stan_data$x_sca_l1),FUN = function(i){ # Each row
    sapply(1:draws, function(s){ # Each posterior sample
      sds <- as.vector(exp(dat$stan_data$x_sca_l2[i,,drop=FALSE] %*% samps$eta[,,s]))
      sds.diag <- diag(sds,length(sds),length(sds))
      var.mu <- sds[1]^2
      gamma_group <- dat$stan_data$x_sca_l2[i,,drop=FALSE] %*% samps$gamma[,,s]
      if(grouping_available){
        if(group_known[i]){
          # Get group_numeric for this row from previous fit.
          group_numeric_i <- object$group_map$group_L2[object$group_map$group_L2[,1] == dat$group_map$group_L1[i,1], 'group_numeric']
          # Add fixed to RE_z(diag(sd)L).
          gamma_group <- gamma_group + (samps$random_z[,,s][group_numeric_i,,drop=FALSE] %*% t(sds.diag%*%t(chol(samps$Omega[,,s]))))[-1]
        } else {
          gamma_group <- gamma_group + mvtnorm::rmvnorm(1,sigma=sds.diag%*%samps$Omega[,,s]%*%sds.diag)[-1]
        }
      }
      shat <- exp(dat$stan_data$x_sca_l1[i,,drop=FALSE] %*% t(gamma_group))
      if(!is.null(occasion)){ # Composite/avg score ICC
        icc <- var.mu / (var.mu + shat^2/occasion[i])
      } else { # Raw score ICC
        icc <- var.mu / (var.mu + shat^2)
      }
      icc
    })
  })
  colnames(out) <- paste0('icc[',1:nrow(newdata),']')
  if(summary){
    L <- (1-prob)/2
    U <- 1 - L
    out <- as.data.frame(cbind(mean=colMeans(out),t(apply(out,2,function(x){
      quantile(x,probs=c(L,U))
    }))))
    if(inc_group & grouping_available){
      out[,fnames$grouping] <- newdata[,fnames$grouping]
      out[out[,fnames$grouping] == magic_NA,fnames$grouping] <- NA
    }
  }
  return(out)

}

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
    draws <- nsamples(object)
  }
  K <- object$stan_data$K
  P_l1 <- object$stan_data$P_l1
  P_l2 <- object$stan_data$P_l2
  X <- object$stan_data$x_sca_l2

  # Grab (formatted) matrices
  samps.mu <- as.matrix(object$fit,pars='beta0')[1:draws,,drop=FALSE]
  samps.gamma <- array(t(as.matrix(object$fit,pars='gamma')[1:draws,]),dim=c(P_l2,P_l1,draws))
  samps.mu_group <- as.matrix(object$fit,pars='mu_group')[1:draws,,drop=FALSE]
  samps.mu_group <- array(t(samps.mu_group),dim=c(K,1,draws))
  samps.gamma_group <- as.matrix(object$fit,pars='gamma_group')[1:draws,,drop=FALSE]
  samps.gamma_group <- array(t(samps.gamma_group),dim=c(K,P_l1,draws))

  # Solve for RE
  group_mu_random <- sapply(1:draws,FUN=function(x){
    samps.mu_group[,,x] - samps.mu[x,1]
  },simplify='array')
  group_mu_random <- array(group_mu_random,dim=c(K,1,draws))

  group_gamma_random <- sapply(1:draws,FUN=function(x){
    samps.gamma_group[,,x] - X%*%samps.gamma[,,x]
  },simplify = 'array')
  group_gamma_random <- array(group_gamma_random,dim=c(K,P_l1,draws))

  colnames(group_mu_random) <- 'Mean'
  colnames(group_gamma_random) <- colnames(object$stan_data$x_sca_l1)

  list(mu_random = group_mu_random, gamma_random = group_gamma_random)
}

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
    draws <- nsamples(object)
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
  samps <- as.matrix(object$fit, pars = c('gamma','eta','Omega'))[1:draws,]
  gamma.cols <- grep('gamma.*',colnames(samps),value = TRUE)
  eta.cols <- grep('eta.*',colnames(samps),value = TRUE)
  gamma_group.cols <- grep('gamma_group.*',colnames(samps),value = TRUE)
  Omega.cols <- grep('Omega.*',colnames(samps),value = TRUE)

  gamma <- array(t(samps[,gamma.cols]),dim=c(P_l2,P_l1,draws))
  eta <- array(t(samps[,eta.cols]),dim=c(P_l2,P_l1 + 1,draws))
  Omega <- array(t(samps[,Omega.cols]),dim=c(P_l1 + 1,P_l1 + 1,draws))
  random_z <- .get_random_effect_z_samples(object,draws)
  return(mget(c('gamma','eta','random_z','Omega')))
}

#' ICCier method to extract ICC values.
#'
#' @param object ICCier object
#' @param summary Logical. Whether to return summary (mean, intervals) of ICCs (TRUE), or posterior samples (FALSE)
#' @inheritParams posterior_interval.ICCier
#' @param inc_group Logical. Whether to include the grouping variable with the estimates.
#' @param occasion Default: NULL. Vector representing the occasion number. One value per row in newdata. For estimating composite score reliability. If unspecified, set to NULL (default), and raw score ICCs are estimated.
#' @param ... Not currently used.
#'
#' @return If \code{summary=TRUE}, then the mean and \code{prob}\% intervals are returned for each observation in the model frame.
#' If \code{summary=FALSE}, an S by N matrix containing the S posterior samples for N observations.
#' @export
#'
fitted.ICCier <- function(object, summary=TRUE, prob=.95,inc_group=TRUE,occasion=NULL){
  if(summary) {
    out <- as.data.frame(cbind(mean=matrix(.posterior_mean(object,pars='icc'),ncol=1),posterior_interval(object,prob=prob,pars='icc')))
    colnames(out)[1] <- 'mean'
    if(inc_group){
      out[,.get_formula_names(object)$grouping] <- object$group_map$group_L1[,1]
    }
  } else {
    out <- as.matrix(object$fit,pars='icc')
  }

  if(!is.null(occasion)){
    out <- predict.ICCier(object,object$data,summary=summary,prob=prob,inc_group=inc_group,occasion=occasion)
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


#' ICCier method for extracting group-specific values
#'
#' Takes ICCier object and returns group-specific values.
#'
#' @param object ICCier object
#' @inheritParams fitted.ICCier
#'
#' @return
#' @export
#'
#' @examples
coef.ICCier <- function(object,summary = TRUE,prob = .95){

}

#' ICCier method for extracting random effect values
#'
#' Takes ICCier object and returns random effects.
#'
#' @param object ICCier object
#' @inheritParams fitted.ICCier
#'
#' @return
#' @export
#'
#' @examples
ranef.ICCier <- function(object,summary = TRUE, prob = .95){

}

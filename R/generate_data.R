#' Generate data for testing.
#'
#' Generates data for use in testing.
#' Takes a formula for the model spec (See \code{\link{ICCier}}), returns list of params, metadata, and stan data.
#'
#' @param formula Formula. See \code{\link{ICCier}}.
#' @param n Number of observations per group.
#' @param K Number of groups.
#' @param beta Location model fixed coefficients. Q_L2xQ_L1.
#' @param gamma Scale model fixed coefficients. P_L2xP_l1.
#' @param eta Between-group coefficients. R_L2x(Q_l1 + P_l1).
#' @param cor_structure 3-length vector describing correlation matrix of REs.
#'
#' @return ICCier_dataget object. List of params, meta, and data.
#' @import Formula
#' @keywords internal
datagen.formula <- function(formula,n=NA,K=NA,beta,gamma,eta,cor_structure,...){
  dots <- list(...)
  N <- n*K
  formula <- as.Formula(formula)
  Q_l1 <- ncol(beta)
  P_l1 <- ncol(gamma)
  Q_l2 <- nrow(beta)
  P_l2 <- nrow(gamma)
  R_l2 <- nrow(eta)

  # Get names
  outcome <- all.vars(formula(formula,lhs=1,rhs=0))
  group <- all.vars(formula(formula,lhs=2,rhs=0))
  L2 <- all.vars(formula(formula,lhs=0,rhs=c(2,3,5)))
  L1 <- all.vars(formula(formula,lhs=0,rhs=c(1,4)))

  # Branch: If data is supplied, use that to create model matrices.
  if(!is.null(dots$data)){
    ds <- dots$data
    ds[,outcome] <- 1
    pf <- .parse_formula(formula,ds)
    data <- pf$stan_data
    group_map <- pf$group_map
    N <- nrow(ds)
    K <- length(unique(ds[,group]))
  } else {
    # Create L2 model frame
    x_l2 <- data.frame(group = 1:K)
    colnames(x_l2)[1] <- group
    if(length(L2) > 0){
      x_l2[,L2] <- mvtnorm::rmvnorm(K,sigma=diag(1,length(L2)))
    }

    # Create L1 model frame
    x_l1 <- data.frame(group=rep(1:K,each=n))
    colnames(x_l1)[1] <- group
    if(length(L1) > 0){
      x_l1[,L1] <- mvtnorm::rmvnorm(N,sigma=diag(1,length(L1)))
    }

    x_sca_l1 <- model.matrix(formula(formula,rhs=1,lhs=0),x_l1)
    x_sca_l2 <- model.matrix(formula(formula,rhs=2,lhs=0),x_l2)
    x_bet_l2 <- model.matrix(formula(formula,rhs=3,lhs=0),x_l2)
    x_loc_l1 <- model.matrix(formula(formula,rhs=4,lhs=0),x_l1)
    x_loc_l2 <- model.matrix(formula(formula,rhs=5,lhs=0),x_l2)

    group_map <- list(group_L1=model.frame(formula,x_l1,lhs=2,rhs=0),group_L2=model.frame(formula,x_l2,lhs=2,rhs=0))
    data <- mget(c('x_sca_l1','x_sca_l2','x_bet_l2','x_loc_l1','x_loc_l2','N','n','K','Q_l1','Q_l2','P_l1','P_l2','R_l2'))

  }


  # Package up spec
  params <- list(beta=beta,gamma=gamma,eta=eta,Omega=.omegagen(cor_structure,P_l1,Q_l1),betaGamma_random_cor_L = t(chol(.omegagen(cor_structure,P_l1,Q_l1))))
  meta <- list(n=n,K=K,N=N,Q_l1=Q_l1,P_l1=P_l1,Q_l2=Q_l2,P_l2=P_l2,R_l2=R_l2,group_map=group_map,outcome = outcome,formula=formula)
  data$group = group_map$group_L1[[1]]

  spec <- list(params=params,meta=meta,data=data)
  return(.datagen(spec))

}

#' Generate Omega matrix
#'
#' Generates Omega matrix from a cor_structure specification
#'
#' @param cor_structure 3-length Vector of (loc-loc, loc-sca, sca-sca) correlations.
#' @param P_l1 Number of random scale effects
#' @param Q_l1 Number of random loc effects
#'
#' @return Correlation matrix
#' @keywords internal
.omegagen <- function(cor_structure,P_l1,Q_l1){
  Omega <- matrix(0,P_l1 + Q_l1, P_l1 + Q_l1)
  diag(Omega) <- 1
  Omega[1:Q_l1,1:Q_l1][lower.tri(Omega[1:Q_l1,1:Q_l1])] <- cor_structure[1]
  Omega[(Q_l1 + 1):nrow(Omega),(Q_l1 + 1):ncol(Omega)][lower.tri(Omega[(Q_l1 + 1):nrow(Omega),(Q_l1 + 1):ncol(Omega)])] <- cor_structure[3]
  Omega[(Q_l1 + 1):nrow(Omega), 1:Q_l1] <- cor_structure[2]
  Omega <- Matrix::forceSymmetric(Omega,'L')
  Omega <- as.matrix(Omega)
  return(Omega)
}

#' Datagen helper.
#'
#' The spec consists of everything but the random effects and outcome.
#' Specifically, the spec should be a list (params, meta, data).
#' Params should contain the beta, gamma, eta, Omega, betaGamma_random_cor_L matrices.
#' Meta should contain all dimensions (e.g., n, K, N, P_l1, etc), the group_map, outcome name, and formula.
#' This is used for both this function, and dataframe converters, to fill in variable names.
#' Data should contain fixed design matrices and other stan data. \code{.datagen} is responsible for generating RE_z,RE,group-specific, and y.
#'
#' @param spec List containing params, meta, and data. See details.
#'
#' @return Completed datagen list, containing params, meta, and data.
#'
#' @keywords internal
.datagen <- function(spec){

  # Generate REs
  betaGamma_random_z <- mvtnorm::rmvnorm(spec$meta$K,sigma=diag(1,spec$meta$P_l1 + spec$meta$Q_l1,spec$meta$P_l1 + spec$meta$Q_l1))
  betaGamma_random <- matrix(0,nrow(betaGamma_random_z),ncol(betaGamma_random_z))
  betaGamma_random_sigma <- exp(spec$data$x_bet_l2 %*% spec$params$eta)
  for(i in 1:spec$meta$K){
    betaGamma_random[i,] <- betaGamma_random_z[i,]%*%t(diag(betaGamma_random_sigma[i,],length(betaGamma_random_sigma[i,]),length(betaGamma_random_sigma[i,])) %*% spec$params$betaGamma_random_cor_L)
  }

  # Generate group-specific
  betaGamma_group <- cbind(spec$data$x_loc_l2 %*% spec$params$beta, spec$data$x_sca_l2 %*% spec$params$gamma) + betaGamma_random
  beta_group <- betaGamma_group[,1:spec$meta$Q_l1,drop=FALSE]
  gamma_group <- betaGamma_group[,(spec$meta$Q_l1 + 1):ncol(betaGamma_group),drop=FALSE]

  # Generate outcome
  yhat <- rowSums(spec$data$x_loc_l1 * beta_group[spec$data$group,,drop=FALSE])
  shat <- exp(rowSums(spec$data$x_sca_l1 * gamma_group[spec$data$group,,drop=FALSE]))
  y <- yhat + rnorm(spec$meta$N,0,shat)
  icc <- betaGamma_random_sigma[,1]^2 / (betaGamma_random_sigma[,1]^2 + shat^2)

  spec$params$betaGamma_random_z <- betaGamma_random_z
  spec$params$betaGamma_random <- betaGamma_random
  spec$params$betaGamma_random_sigma <- betaGamma_random_sigma
  spec$params$beta_group <- beta_group
  spec$params$gamma_group <- gamma_group
  spec$data$y <- y
  spec$params$icc <- icc

  class(spec) <- 'ICCier_datagen'

  return(spec)

}

#' Generate data for testing. (Deprecated)
#'
#' Generates data for use in testing. (Deprecated: Use datagen.formula() instead.)
#'
#' @param n Number of observations per person.
#' @param K  Number of persons.
#' @param beta 1x1 matrix representing the global mean of the outcome.
#' @param gamma QxQ_random matrix representing the effects of Q level 2 variables on the Q_random level 1 scale parameters.
#' @param eta Qx(Q_random+1) matrix representing the effects of Q level 2 variables on the mean variance and Q_random level 1 scale variances.
#' @param cor_structure 3-length vector describing the correlation matrix of random effects.
#'
#' @importFrom mvtnorm rmvnorm
#' @importFrom Matrix forceSymmetric
#'
#' @return ...
#'
#' @keywords internal
datagen <- function(n,K,beta,gamma,eta,cor_structure){
  N <- n*K
  Q <- nrow(beta)
  P <- nrow(gamma)
  Q_random <- ncol(beta)
  P_random <- ncol(gamma)
  # P_random <- length(beta_sigma)
  # Q_random <- length(gamma_sigma)
  group <- rep(1:K,each=n)

  if(ncol(beta) != Q_random | ncol(gamma) != P_random) {
    stop('Mismatch between number of beta/gamma cols and implied number of random effects in beta_sigma/gamma_sigma\n.')
  }

  ## Create cor matrix
  Omega <- matrix(0,P_random + Q_random, P_random + Q_random)
  diag(Omega) <- 1
  Omega[1:Q_random,1:Q_random][lower.tri(Omega[1:Q_random,1:Q_random])] <- cor_structure[1]
  Omega[(Q_random + 1):nrow(Omega),(Q_random + 1):ncol(Omega)][lower.tri(Omega[(Q_random + 1):nrow(Omega),(Q_random + 1):ncol(Omega)])] <- cor_structure[3]
  Omega[(Q_random + 1):nrow(Omega), 1:Q_random] <- cor_structure[2]
  Omega <- forceSymmetric(Omega,'L')
  Omega <- as.matrix(Omega)
  betaGamma_random_cor_L <- t(chol(Omega)) # Confusingly, this is now L tri, not U tri like usual. The hell, R?

  ## Generate L1 effects
  X_loc_l2 <- mvtnorm::rmvnorm(K,sigma=diag(1,Q,Q))
  X_sca_l2 <- mvtnorm::rmvnorm(K,sigma=diag(1,P,P))
  X_loc_l2[,1] <- 1
  X_sca_l2[,1] <- 1

  betaGamma_random_z <- mvtnorm::rmvnorm(K,sigma=diag(1,Q_random + P_random, Q_random + P_random))
  betaGamma_random <- matrix(0,K,P_random + Q_random)
  betaGamma_random_sigma = exp(X_sca_l2 %*% eta)
  for(i in 1:nrow(betaGamma_random)){
    betaGamma_random[i,] = betaGamma_random_z[i,]%*%t(diag(betaGamma_random_sigma[i,]) %*% betaGamma_random_cor_L)
  }
  betaGamma <- cbind(X_loc_l2%*%beta, X_sca_l2%*%gamma) + betaGamma_random
  beta_group <- betaGamma[,1:Q_random,drop=FALSE]
  gamma_group <- betaGamma[,(Q_random + 1):ncol(betaGamma),drop=FALSE]


  ## Generate L1 observations
  X_loc_l1 <- mvtnorm::rmvnorm(N,sigma=diag(1,nrow=Q_random,ncol=Q_random))
  X_sca_l1 <- mvtnorm::rmvnorm(N,sigma=diag(1,nrow=P_random,ncol=P_random))
  X_sca_l1[,1] <- 1
  X_loc_l1[,1] <- 1
  # X_loc_l1 <- cbind(1,mvtnorm::rmvnorm(N,sigma=diag(1,nrow=P_random-1,ncol=P_random-1)))
  # X_sca_l1 <- cbind(1,mvtnorm::rmvnorm(N,sigma=diag(1,nrow=Q_random-1,ncol=Q_random-1)))

  y <- rowSums(X_loc_l1 * beta_group[group,,drop=FALSE]) + rnorm(N,0,exp(rowSums(X_sca_l1 * gamma_group[group,,drop=FALSE])))
  icc <- betaGamma_random_sigma[group,1]^2 / (betaGamma_random_sigma[group,1]^2 + exp(rowSums(X_sca_l1 * gamma_group[group,,drop=FALSE]))^2)

  ## Package up
  params <- list(beta=beta,gamma=gamma,eta=eta,Omega=Omega,betaGamma_random_cor_L=betaGamma_random_cor_L,beta_group=beta_group,gamma_group=gamma_group,betaGamma_random=betaGamma_random,betaGamma_random_sigma = betaGamma_random_sigma,icc=icc)
  meta <- list(n=n,K=K,N=N,P=P,Q=Q,P_random=P_random,Q_random=Q_random)
  data <- list(n=n,K=K,N=N,P_l2=P,P_l1=P_random,Q_l1 = Q_random, Q_l2 = Q,y=y,x_loc_l2=X_loc_l2,x_sca_l2=X_sca_l2,x_loc_l1=X_loc_l1,x_sca_l1=X_sca_l1,group=group)
  list(params=params,meta=meta,data=data)

}

#' Generate data frame for testing. (Deprecated)
#'
#' Uses same as \code{\link{datagen}}, but outputs a data frame.
#'
#' @param n n
#' @param K K
#' @param beta beta
#' @param gamma gamma
#' @param eta eta
#' @param cor_structure cor_structure
#'
#' @return ...
#'
#' @keywords internal
generate_df <- function(n,K,beta,gamma,eta,cor_structure){
  d <- datagen(n,K,beta,gamma,eta,cor_structure)
  ds <- convert_datagen(d)
  return(ds)
}

#' Convert datagen to data frame. (Deprecated)
#'
#' @param d datagen output.
#'
#' @return data.frame.
#' @keywords internal
convert_datagen <- function(d){
  ds <- data.frame(y=d$data$y,group=d$data$group)
  if(d$meta$P_random > 1){
    ds[,paste0('sca_l1.',1:(d$meta$P_random - 1))] <- d$data$x_sca_l1[,2:d$meta$P_random]
  }
  if(d$meta$Q_random > 1){
    ds[,paste0('loc_l1.',1:(d$meta$Q_random - 1))] <- d$data$x_loc_l1[,2:d$meta$Q_random]
  }
  if(d$meta$P > 1){
    ds[,paste0('sca_l2.',1:(d$meta$P - 1))] <- d$data$x_sca_l2[d$data$group,2:d$meta$P]
  }
  if(d$meta$Q > 1){
    ds[,paste0('loc_l2.',1:(d$meta$Q - 1))] <- d$data$x_loc_l2[d$data$group,2:d$meta$Q]
  }
  return(ds)
}

#' Convert datagen to data.frame
#'
#' @param d datagen output.
#'
#' @return data.frame.
#' @keywords internal
as.data.frame.ICCier_datagen <- function(d){
  ds <- data.frame(y=d$data$y,group=d$data$group)
  colnames(ds) <- c(d$meta$outcome,colnames(d$meta$group_map$group_L1)[1])

  names_l1 <- list(x_sca_l1=colnames(d$data$x_sca_l1),x_loc_l1=colnames(d$data$x_loc_l1))
  names_l2 <- list(x_sca_l2=colnames(d$data$x_sca_l2),x_loc_l2=colnames(d$data$x_loc_l2),x_bet_l2=colnames(d$data$x_bet_l2))

  ds[,names_l1[['x_sca_l1']]] <- d$data$x_sca_l1
  ds[,names_l1[['x_loc_l1']]] <- d$data$x_loc_l1
  ds[,names_l2[['x_sca_l2']]] <- d$data$x_sca_l2[d$data$group,]
  ds[,names_l2[['x_loc_l2']]] <- d$data$x_loc_l2[d$data$group,]
  ds[,names_l2[['x_bet_l2']]] <- d$data$x_bet_l2[d$data$group,]
  ds[names(ds) == '(Intercept)'] <- NULL

  return(ds)

}

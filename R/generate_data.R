#' Generate data for testing.
#'
#' Generates data for use in testing.
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

#' Generate data frame for testing.
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

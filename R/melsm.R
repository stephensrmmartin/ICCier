
#' Run MELSM model.
#'
#' Runs the MELSM model and returns an ICCier object.
#'
#' \code{ICCier} uses the mixed effects location scale model to estimate an intercept-only
#' location model with random effect of person.
#' The within-person variances are log-linearly modelled from a set of observation-level and
#' person-level predictors, with coefficients \eqn{\gamma}.
#' The between-person variances are also log-linearly modelled from a set of person-level predictors,
#' with coefficients \eqn{\eta}.
#'
#' The formula syntax is as follows:
#'
#' \code{outcome | personID ~ Level_1_formula | Level_2_formula}
#'
#' The model is implemented in a 'maximal' manner, meaning \emph{all} level 1 effects are assumed to randomly vary and correlate.
#' Moreover, \emph{all} level 2 variables predict each level 1 parameter.
#' This means you \emph{should not include cross-level interaction terms}, because they are implicit in the model formulation.
#' Level 1 or Level 2 interaction terms may be included.
#'
#' For now, the Level 2 formula predicts both the level 1 scale parameters, as well as the level 2 random effect variances.
#' If you wish to only have person-specific IICs, use \code{1} as the level 1 formula (intercept-only).
#'
#' @param formula Formula representing the model. See details.
#' @param data Data frame containing all variables.
#' @param ... Arguments passed to \code{\link[rstan]{sampling}}.
#'
#' @return ICCier object. List containing the model formula, data, stan_data, model fit, and mapping between original ID and numeric ID.
#' @export
#'
ICCier <- function(formula, data, ...){
  d <- .parse_formula(formula, data)
  stanOut <- rstan::sampling(stanmodels$melsmICC,data=d$stan_data,pars=c('beta0','gamma','eta','mu_group','gamma_group','icc','log_lik','Omega'),...)
  out <- list(formula=Formula(formula), data=d$model.frame, stan_data = d$stan_data,fit=stanOut, group_map = d$group_map)
  class(out) <- c('ICCier')
  return(out)
}

#' Parses formula using Formula.
#'
#' @import Formula
#'
#' @param formula formula
#' @param data data
#'
#' @return List containing stan data, group mapping, and the model frame.
#' @keywords internal
#'
.parse_formula <- function(formula, data){
  n_orig <- nrow(data)

  f <- Formula::Formula(formula)
  length.f <- length(f)
  if(length.f[1] != 2){
    stop('Both the outcome and person-level indicator variables must be specified.')
  }
  if(length.f[2] != 2){
    stop('Both the level 1 and level 2 formulas must be specified.')
  }

  mf <- model.frame(f, data,na.action='na.omit')
  if(nrow(mf) - nrow(n_orig) > 0){
    message(paste0('Dropping ',nrow(mf) - n_orig,' incomplete cases.'))
  }

  group <- model.frame(f,mf,lhs=2,rhs=0,drop.unused.levels = TRUE)
  group$group_numeric <- as.numeric(as.factor(group[,1]))

  mf[,2] <- group$group_numeric

  N <- nrow(mf)
  K <- length(unique(group$group_numeric))

  x_sca_l1 <- model.matrix(f, mf, rhs=1)
  P_l1 <- ncol(x_sca_l1)

  x_sca_l2.mf <- model.frame(f,mf,lhs=2,rhs=2)
  x_sca_l2 <- as.data.frame(do.call(rbind,lapply(split(x_sca_l2.mf,f=group$group_numeric),function(x){x[1,]})))
  x_sca_l2 <- x_sca_l2[order(x_sca_l2[,1]),-1,drop=FALSE]
  x_sca_l2 <- model.matrix(f,x_sca_l2,rhs=2)
  P_l2 <- ncol(x_sca_l2)

  y <- mf[,1]
  stan_data <- mget(c('N','K','P_l1','P_l2','x_sca_l1','x_sca_l2','y'))
  stan_data$group <- group$group_numeric
  return(list(stan_data=stan_data,group_map = group, model.frame = mf))
}

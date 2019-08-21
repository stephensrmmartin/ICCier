
#' Run MELSM model.
#'
#' Runs the MELSM model and returns an ICCier object.
#'
#' \code{ICCier} uses the mixed effects location scale model to estimate an unconditional
#' (intercept-only) or conditional location model with random effect of person.
#' The within-person variances (i.e., residual, or error variances) are log-linearly modelled from a set of observation-level and
#' person-level predictors, with coefficients \eqn{\gamma}.
#' The between-person variances are also log-linearly modelled from a set of person-level predictors,
#' with coefficients \eqn{\eta}.
#'
#' The formula syntax is as follows:
#'
#' For an unconditional model:
#' \code{outcome | personID ~ Level_1_formula | Level_2_formula}
#'
#' For a conditional model:
#' \code{outcome | personID ~ Level_1_formula | Level_2_formula | Level_1_conditional | Level_2_conditional}
#'
#' The model is implemented in a 'maximal' manner, meaning \emph{all} level 1 effects are assumed to randomly vary and correlate.
#' Moreover, \emph{all} level 2 variables predict each level 1 parameter.
#' This means you \emph{should not include cross-level interaction terms}, because they are implicit in the model formulation.
#' Level 1 or Level 2 interaction terms may be included.
#'
#' For now, the Level 2 formula predicts both the level 1 scale parameters, as well as the level 2 random effect variances.
#' If you wish to only have person-specific ICCs, use \code{1} as the level 1 formula (intercept-only).
#'
#' When using a \emph{conditional} model, the default ICCs are "unadjusted".
#' In our case, the "unadjusted" ICC is still the random intercept variance divided by
#' the random intercept variance and error variance.
#' This makes sense, if the location model(s) are meant to be controlling variables.
#' If you wish to have the so-called \emph{adjusted} ICC, use \code{adjusted = TRUE}.
#' The \emph{adjusted} ICC instead uses the expected variance due to \emph{all random factors},
#' divided by itself and the error variance.
#' The adjusted ICC is therefore the proportion of random variance due to the random factors.
#' This makes sense if the goal is to examine individual differences, and therefore examine the
#' random effects themselves.
#'
#' @param formula Formula representing the model. See details.
#' @param data Data frame containing all variables.
#' @param ... Arguments passed to \code{\link[rstan]{sampling}}.
#' By default, \code{sampling} is called with \code{chains=4,iter=2000,adapt_delta=.95,init=0}.
#' If \code{options('mc.cores')} is defined, then \code{cores=getOption("mc.cores")}; otherwise all detected cores are used.
#'
#' @return ICCier object. List containing the model formula, data, stan_data, model fit, and mapping between original ID and numeric ID.
#' @importFrom parallel detectCores
#' @export
#'
ICCier <- function(formula, data, ...){
  # Sane defaults
  dots <- list(...)
  if(is.null(dots$control)){
    dots$control <- list(adapt_delta=.95)
  }
  if(is.null(dots$control$adapt_delta)){
    dots$control$adapt_delta = .95
  }
  if(is.null(dots$init)){
    dots$init <- 0
  }
  if(is.null(dots$cores)){
    cores <- getOption('mc.cores')
    if(is.null(cores)){
      dots$cores <- parallel::detectCores()
    }
  }
  if(is.null(dots$chains)){
    dots$chains <- 4
  }

  d <- .parse_formula(formula, data)
  if(d$conditional){
    model <- stanmodels$melsmCondICC
  } else {
    model <- stanmodels$melsmICC
  }
  if(!is.null(dots$adjusted)){
    adjusted <- dots$adjusted
    if(!is.logical(adjusted)) {
      stop('"adjusted" must be TRUE/FALSE.')
    }
    if(!d$conditional){
      warning("'adjusted' is not applicable when using an unconditional model. Ignoring argument.")
      adjusted <- FALSE
    }
    d$stan_data$adjust_icc <- as.numeric(adjusted)
    dots$adjusted <- NULL
  } else {
    adjusted <- FALSE
    d$stan_data$adjust_icc <- as.numeric(adjusted)
  }

  args <- c(list(object=model, data=d$stan_data,
            pars = c('beta0','gamma','eta','mu_group','gamma_group','icc','log_lik','Omega','icc_mean','icc_sd')),
            dots)
  stanOut <- do.call('sampling',args=args)

  out <- list(formula=Formula(formula), data=d$model.frame, stan_data = d$stan_data,fit=stanOut, group_map = d$group_map,type=list(conditional=d$conditional,adjusted=adjusted))

  diagnostics <- .get_diagnostics(out)
  out$diagnostics <- diagnostics

  class(out) <- c('ICCier')
  return(out)
}

#' Parses formula using Formula.
#'
#' @import Formula
#'
#' @param formula formula
#' @param data data
#' @param predict Logical. Default: FALSE. If True, y is ignored, and the data requires grouping variable and predictors.
#'
#' @return List containing stan data, group mapping, and the model frame.
#' @keywords internal
#'
.parse_formula <- function(formula, data,predict=FALSE){
  n_orig <- nrow(data)

  f <- Formula::Formula(formula)
  length.f <- length(f)
  conditional <- length.f[2] > 2

  if(length.f[1] != 2){
    stop('Both the outcome and person-level indicator variables must be specified.')
  }
  if(length.f[2] < 2){
    stop('Both the level 1 and level 2 formulas must be specified.')
  }
  if(length.f[2] < 4 & conditional){
    stop('If specifying a location model, both level 1 and level 2 formulas must be specified.')
  }
  fnames <- attr(terms(f),'term.labels')

  if(predict){
    # No outcome
    mf <- model.frame(f,data,na.action='na.omit',lhs=2)
  } else {
    mf <- model.frame(f, data,na.action='na.omit')
  }
  if(n_orig - nrow(mf) > 0){
    message(paste0('Dropping ',n_orig - nrow(mf) ,' incomplete cases.'))
  }

  group_L1 <- model.frame(f,mf,lhs=2,rhs=0,drop.unused.levels = TRUE)
  group_L1$group_numeric <- as.numeric(as.factor(group_L1[,1]))
  group_L2 <- as.data.frame(do.call(rbind, lapply(split(group_L1, f=group_L1$group_numeric), FUN=function(x){x[1,]})))
  group <- list(group_L1 = group_L1, group_L2= group_L2)

  mf[,fnames[2]] <- group$group_L1$group_numeric

  N <- nrow(mf)
  K <- length(unique(group$group_L1$group_numeric))

  x_sca_l1 <- model.matrix(f, mf, rhs=1)
  P_l1 <- ncol(x_sca_l1)

  if(conditional){
    f.location <- Formula(formula(f,rhs=c(3,4),lhs=0))
  } else{
    f.location <- Formula(~1|1)
  }
  x_loc_l1 <- model.matrix(f.location,mf,rhs=1)
  Q_l1 <- ncol(x_loc_l1)

  if(predict){
    # Leave matrix as-is for prediction.
    x_sca_l2 <- model.matrix(f,mf,rhs=2)
    x_loc_l2 <- model.matrix(f.location,mf,rhs=2)
    group$group_L2 <- group$group_L1
  } else {
    x_l2.mf <- model.frame(f,mf,lhs=2)
    x_l2 <- as.data.frame(do.call(rbind,lapply(split(x_l2.mf,f=group$group_L1$group_numeric),function(x){x[1,]})))
    x_l2 <- x_l2[order(x_l2[,1]),-1,drop=FALSE]

    x_sca_l2 <- model.matrix(f,x_l2,rhs=2)
    x_loc_l2 <- model.matrix(f.location,x_l2,rhs=2)
  }
  P_l2 <- ncol(x_sca_l2)
  Q_l2 <- ncol(x_loc_l2)

  if(predict){
    y <- NA
  } else {
    y <- mf[,fnames[1]]
  }
  stan_data <- mget(c('N','K','P_l1','P_l2','x_sca_l1','x_sca_l2','y','Q_l1','Q_l2','x_loc_l1','x_loc_l2'))
  stan_data$group <- group$group_L1$group_numeric
  mf[,fnames[2]] <- group$group_L1[,fnames[2]]
  return(list(stan_data=stan_data,group_map = group, model.frame = mf,conditional=conditional))
}

.get_diagnostics <- function(object){
  rhats <- rstan::summary(object$fit,pars=c('beta0','gamma','eta','mu_group','gamma_group'))$summary[,'Rhat']

  n_effs <- rstan::summary(object$fit,pars=c('beta0','gamma','eta','mu_group','gamma_group'))$summary[,'n_eff']

  div <- rstan::get_num_divergent(object$fit)

  tree.max <- rstan::get_num_max_treedepth(object$fit)

  bfmi <- rstan::get_bfmi(object$fit)

  return(mget(c('rhats','n_effs','div','tree.max','bfmi')))

}

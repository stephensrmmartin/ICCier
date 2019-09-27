#' Fit ICCier model.
#'
#' Runs the ICCier model and returns an ICCier object.
#'
#' The ICC is computed from a between-group variance and a within-group variance.
#' Unlike the traditional ICC, \code{ICCier} allows the within-group variance to vary across groups.
#' Moreover, \code{ICCier} can model both the between-group variance and within-group variance.
#' This effectively provides insight into predictors of reliability.
#'
#' \code{ICCier} uses the mixed effects location scale model (MELSM) to model the variance components required for the ICC.
#' The model specification is detailed in the \code{Model Specification} section below.
#' By default, unconditional ICC(1) values are estimated.
#' To obtain ICC(2) values, use the \code{fitted} function with the "occasion" argument specified.
#'
#' Conditional ICCs are estimated if a mean model is specified.
#' By default, the conditional ICC is "unadjusted".
#' To obtain "adjusted" ICCs, specify \code{adjusted=TRUE}.
#' See the vignette for more information.
#' In our case, the "unadjusted" ICC is still the random intercept variance divided by
#' the random intercept variance and error variance.
#' This makes sense, if the location model(s) are meant to be controlling variables.
#' The \emph{adjusted} ICC instead uses the expected variance due to \emph{all random factors},
#' divided by itself and the error variance.
#' The adjusted ICC is therefore the proportion of random variance due to the random factors.
#' This makes sense if the goal is to examine individual differences, and therefore examine the
#' random effects themselves.
#'
#' Note that \code{ICCier} estimates a \emph{maximal} model, meaning \emph{all} within-group predictors randomly vary and correlate.
#' Moreover, \emph{all} between-group predictors predict within-group coefficients.
#' This means you \emph{should not include cross-level interaction terms}, because they are implicit in the model formulation.
#' Within-group or between-group interaction terms may be included.
#'
#' \code{ICCier} uses the mixed effects location scale model (MELSM) to estimate an unconditional
#' (intercept-only) or conditional location model with random effect of person.
#' The within-person variances (i.e., residual, or error variances) are log-linearly modelled from a set of observation-level and
#' person-level predictors, with coefficients \eqn{\gamma}.
#' The between-person variances are also log-linearly modelled from a set of person-level predictors,
#' with coefficients \eqn{\eta}.
#'
#' \subsection{Model specification}{
#' The \code{ICCier} model imposes a model on the between-group SDs, within-group SDs, and the mean if a conditional ICC is desired.
#' These models can be specified in two ways.
#' One way is through a multiple-formula syntax, and another through a single-formula syntax.
#' These are described in turn.
#'
#' \subsection{Multiple-formula Syntax}{
#'     The multiple-formula syntax at minimum requires \code{x} (outcome) and \code{group} (the grouping variable).
#'     These are specified as raw variable names, not as character strings (i.e., \code{x=recall}, \emph{not} \code{x='recall'}).
#'     If only these two are included, then the between and within-group variances are modeled with intercepts only.
#'     The formulas are exemplified as follows:
#'     \describe{
#'       \item{\code{between ~ Sex + Age}}{The Between-group SD is modeled from between-group variables Sex and Age (and an intercept)}
#'       \item{\code{within ~ trial + condition | Sex + Age}}{The Within-group SD is modeled from within-group variables trial and condition (and an intercept), and between-group predictors Sex and Age.}
#'       \item{\code{mean ~ trial + condition | 1}}{Conditional Model. The mean is modeled from within-group variables trial and condition (and an intercept), and no between-group covariates.}
#'     }
#'     Any of these can be excluded or included.
#'     If they are excluded, they default to intercept-only models.
#' }
#' \subsection{Single-formula Syntax}{
#'     The single formula syntax can take two forms.
#'     One form is shorter and describes an unconditional model.
#'     The second, longer form describes a conditional model.
#'     In total, the full formula consists of:
#'     \enumerate{
#'       \item Outcome (e.g., recall)
#'       \item Grouping variable (e.g., subject)
#'       \item Within-group predictors of within-group variance (e.g., trial and condition)
#'       \item Between-group predictors of within-group variance (e.g., Sex and Age)
#'       \item Between-group predictors of between-group variance (e.g., Sex and Age)
#'       \item Within-group predictors of the mean. For conditional models. (e.g., trial and condition)
#'       \item Between-group predictors of the mean. For conditional models. (e.g., intercept only)
#'     }
#'     For an unconditional model:
#'
#'     \code{outcome | group ~ WG_Level_1_formula | WG_Level_2_formula | BG_Level_2_formula}
#'
#'     For a conditional model:
#'
#'     \code{outcome | group ~ WG_Level_1_formula | WG_Level_2_formula | BG_Level_2_formula | Mean_Level_1 | Mean_Level_2}
#' }
#' }
#'
#' @param x The outcome variable. Raw variable name (not a string).
#' @param group The grouping variable. Raw variable name (not a string).
#' @param formula Formula representing the model. See details.
#' @param data Data frame containing all variables.
#' @param ... Multiple formulas for ICCier model (See Details). Arguments passed to \code{\link[rstan]{sampling}}.
#' By default, \code{sampling} is called with \code{chains=4,iter=2000,adapt_delta=.95,init=0}.
#' If \code{options('mc.cores')} is defined, then \code{cores=getOption("mc.cores")}; otherwise all detected cores are used.
#'
#' @return ICCier object. List containing the model formula, data, stan_data, model fit, and mapping between original ID and numeric ID.
#' @importFrom parallel detectCores
#' @export
#'
ICCier <- function(x,...){
  if(missing(x) | class(substitute(x)) == 'call'){
    x <- formula()
  } else {
    x <- substitute(x)
  }
  UseMethod('ICCier',x)
}

#' Default ICCier Method
#'
#' @export
#' @describeIn ICCier Multiple formula method.
ICCier.default <- function(x, group, data, ...){
  x <- substitute(x)
  group <- substitute(group)
  # if(is.null(x) | is.character(x)){
  if(!is.name(x)){
    stop('x must be specified, as an unquoted variable name.')
  }
  # if(is.null(group) | is.character(group)){
  if(!is.name(group)){
    stop('Group must be specified, as an unquoted variable name.')
  }
  if(!is.data.frame(data)){
    stop('Data must be specified.')
  }

  formList <- list(outcome=x,group=group)

  dots <- list(...)
  which.form <- sapply(dots,is.formula)
  if(any(which.form)){
    forms <- dots[which.form]
    dots[which.form] <- NULL
    formList <- c(formList,forms)
  }

  form <- .formulate(formList)
  do.call(ICCier,c(list(formula=form,data=data),dots))

}

#' Formula ICCier Method.
#'
#' @export
#' @describeIn ICCier Single formula method.
ICCier.formula <- function(formula, data, ...){
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

#' Parses formula list into standard Formula
#'
#' Takes list of formulas, converts into the standard formula.
#' E.g., \code{list(outcome = y, group = subjectID, between ~ age, within ~ day|1, mean ~ 1|1)}
#'
#' @param formList List containing between, within, location, and group formulas.
#'
#' @return Formula.
#' @keywords internal
#'
.formulate <- function(formList){
  if(length(formList) < 2){
    stop('Formula list must have at least outcome (x) and group specified.')
  }

  forms <- list(outcome = formList[['outcome']],group = formList[['group']], between = Formula(between ~ 1), within = Formula(within ~ 1|1), mean = Formula(mean ~ 1|1))

  which.forms <- sapply(formList,is.formula)
  if(any(which.forms)){
    formList[which.forms] <- sapply(formList[which.forms],Formula::as.Formula)
    formType <- sapply(formList[which.forms],function(x){all.vars(formula(x,lhs=1,rhs=0))})
    if(any(!(formType %in% c('between','within','mean')))){
      stop('Only between, within, and mean formulas can be specified.')
    }
    forms[formType] <- formList[which.forms]
  }

  conditional <- length(all.vars(formula(forms[['mean']],lhs=0))) > 0

  # Piece together
  outcome <- deparse(forms[['outcome']])
  group <- deparse(forms[['group']])
  within <- Reduce(paste,deparse(forms[['within']][[3]]))
  between <- Reduce(paste,deparse(forms[['between']][[3]]))
  mean <- Reduce(paste,deparse(forms[['mean']][[3]]))



  fc <- paste(outcome,'|',group,'~',within,'|',between)
  if(conditional){
    fc <- paste(fc,'|',mean)
  }
  form <- Formula::as.Formula(fc)

  return(form)
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
  conditional <- length.f[2] > 3

  if(length.f[1] != 2){
    stop('Both the outcome and grouping variables must be specified.')
  }
  if(length.f[2] < 3){
    # stop('Both the level 1 and level 2 formulas must be specified.')
    stop('The level1 and level2 within-group variance model, and the between-group variance model must be specified.')
  }
  if(length.f[2] < 5 & conditional){
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
    f.location <- Formula(formula(f,rhs=c(4,5),lhs=0))
  } else{
    f.location <- Formula(~1|1)
  }
  x_loc_l1 <- model.matrix(f.location,mf,rhs=1)
  Q_l1 <- ncol(x_loc_l1)

  if(predict){
    # Leave matrix as-is for prediction.
    x_sca_l2 <- model.matrix(f,mf,rhs=2)
    x_loc_l2 <- model.matrix(f.location,mf,rhs=2)
    x_bet_l2 <- model.matrix(f,mf,rhs=3)
    group$group_L2 <- group$group_L1
  } else {
    x_l2.mf <- model.frame(f,mf,lhs=2)
    x_l2 <- as.data.frame(do.call(rbind,lapply(split(x_l2.mf,f=group$group_L1$group_numeric),function(x){x[1,]})))
    x_l2 <- x_l2[order(x_l2[,1]),-1,drop=FALSE]

    x_sca_l2 <- model.matrix(f,x_l2,rhs=2)
    x_loc_l2 <- model.matrix(f.location,x_l2,rhs=2)
    x_bet_l2 <- model.matrix(f,x_l2,rhs=3)
  }
  P_l2 <- ncol(x_sca_l2)
  Q_l2 <- ncol(x_loc_l2)
  R_l2 <- ncol(x_bet_l2)

  if(predict){
    y <- NA
  } else {
    y <- mf[,fnames[1]]
  }
  stan_data <- mget(c('N','K','P_l1','P_l2','x_sca_l1','x_sca_l2','y','Q_l1','Q_l2','x_loc_l1','x_loc_l2','R_l2','x_bet_l2'))
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

#' Tests if object is formula
#'
#' @param x Object to test
#'
#' @return Logical. TRUE if formula
#' @keywords internal
is.formula <- function(x){
  inherits(x,'formula')
}

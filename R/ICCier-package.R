#' The 'ICCier' package.
#'
#' @description ICCier is a package to estimate intraclass coefficients (ICCs) across people, for each person, or even for each individual observation within persons.
#' It does so by using the mixed effects location scale model (MELSM) to model both the within person variance and the between person variance as a function of observation and person-level
#' predictors.
#'
#' @docType package
#' @name ICCier-package
#' @useDynLib ICCier, .registration = TRUE
#' @import methods
#' @import Rcpp
#' @import rstantools
#' @importFrom rstan sampling
#'
#' @references
#' Stan Development Team (2018). RStan: the R interface to Stan. R package version 2.18.2. http://mc-stan.org
#'
NULL

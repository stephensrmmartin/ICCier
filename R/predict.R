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

  fnames <- .get_formula_names(object$formula)
  if(!(fnames$grouping %in% colnames(newdata))){
    grouping_available <- FALSE
    message('No grouping variable. Using fixef only.')
  } else{
    grouping_available <- TRUE
  }

}
# TODO: Needs to handle existing (known) groups as well as unknown.
# If known, pull from the group_map which integer they belong to.
# If unknown, assign a value, then randomly draw from RE distribution(s) instead on each iteration.
# If no groups specified, just use fixed effects, b/c there's no other information available.
# Remember: You only need to predict ICC = var(mean)/(var(mean) + var(within)) for each row, across samples.


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

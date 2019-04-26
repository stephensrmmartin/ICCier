#' ICCier method to predict new ICC values.
#'
#' @param object ICCier object
#' @param data Data to predict from. If NULL, calls \code{\link{fitted.ICCier}} on fit data.
#' @param ... Not currently used.
#'
#' @return TBD
#' @export
#'
predict.ICCier <- function(object, data, ...){
 ## TODO
}

#' ICCier method to extract ICC values.
#'
#' @param object ICCier object
#' @param summary Logical. Whether to return summary (mean, intervals) of ICCs (TRUE), or posterior samples (FALSE)
#' @inheritParams posterior_interval.ICCier
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

% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/predict.R
\name{.extract_transform}
\alias{.extract_transform}
\title{Extracts samples, turns them into formed matrices for prediction}
\usage{
.extract_transform(object, draws)
}
\arguments{
\item{object}{ICCier object.}

\item{draws}{Number of draws}
}
\value{
List of arrays.
}
\description{
Extracts \code{draws} of the generative parameters.
For matrices, converts to array; e.g., gamma is a [P_l2,P_l1,draws] array.
This will ease the computation of predictions, because it can be applied over draws.
}
\keyword{internal}

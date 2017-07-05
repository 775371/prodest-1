\name{FSres}
\alias{FSres,prod-method}
\alias{FSres}
\docType{methods}
%%%%%%%%%%%%%%%% TITLE %%%%%%%%%%%%%%%%%%%
\title{Generate the vector of the first stage residuals}
%%%%%%%%%%%%%%%% DESCRIPTION %%%%%%%%%%%%%%%%%%%
\description{
  This method provides the way to estimate the first stage residuals from a \link{prod} \code{S4} object - estimates from \link{prodestOP}, \link{prodestLP}, \link{prodestACF}, \link{prodestWRDG} and \link{prodestWRDG_GMM} - defined in the \code{prodest} package
}
%%%%%%%%%%%%%%%% USAGE %%%%%%%%%%%%%%%%%%%
\usage{
  FSres(object)
}
%%%%%%%%%%%%%%%% ARGUMENTS %%%%%%%%%%%%%%%%%%%
\arguments{
  %
  \item{object}{
    object of class \link{prod}.
  }
  %
}
%%%%%%%%%%%%%%%% DETAILS %%%%%%%%%%%%%%%%%%%
\details{
  \code{FSres} accepts an \code{S4} \link{prod} object and returns the vector of firs stage residuals.
}
%%%%%%%%%%%%%%%% AUTHOR %%%%%%%%%%%%%%%%%%%
\author{Gabriele Rovigatti}

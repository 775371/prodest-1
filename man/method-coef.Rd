\name{coef}
\alias{coef,prod-method}
\alias{coef}
\docType{methods}
%%%%%%%%%%%%%%%% TITLE %%%%%%%%%%%%%%%%%%%
\title{Print the estimated parameters}
%%%%%%%%%%%%%%%% DESCRIPTION %%%%%%%%%%%%%%%%%%%
\description{
  This method provides the way to extract and print the estimated parameters from a \link{prod} \code{S4} object - estimates from \link{prodestOP}, \link{prodestLP}, \link{prodestACF}, \link{prodestWRDG} and \link{prodestWRDG_GMM} - defined in the \code{prodest} package
}
%%%%%%%%%%%%%%%% USAGE %%%%%%%%%%%%%%%%%%%
\usage{
  coef(object,...)
}
%%%%%%%%%%%%%%%% ARGUMENTS %%%%%%%%%%%%%%%%%%%
\arguments{
  %
  \item{object}{
    object of class \link{prod}.
  }
  %
  \item{...}{
    Additional arguments.
  }
  %
}
%%%%%%%%%%%%%%%% DETAILS %%%%%%%%%%%%%%%%%%%
\details{
  \code{coef} accepts an \code{S4} \link{prod} object and prints the vector of estimated parameters.
}
%%%%%%%%%%%%%%%% AUTHOR %%%%%%%%%%%%%%%%%%%
\author{Gabriele Rovigatti}

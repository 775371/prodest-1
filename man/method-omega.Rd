\name{omega}
\alias{omega,prod-method}
\alias{omega}
\docType{methods}
%%%%%%%%%%%%%%%% TITLE %%%%%%%%%%%%%%%%%%%
\title{Generate the omega estimates}
%%%%%%%%%%%%%%%% DESCRIPTION %%%%%%%%%%%%%%%%%%%
\description{
  This method provides the way to estimate the omega residuals from a \link{prod} \code{S4} object - estimates from \link{prodestOP}, \link{prodestLP}, \link{prodestACF}, \link{prodestWRDG} and \link{prodestWRDG_GMM} - defined in the \code{prodest} package
}
%%%%%%%%%%%%%%%% USAGE %%%%%%%%%%%%%%%%%%%
\usage{
  omega(object)
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
  \code{omega} accepts an \code{S4} \link{prod} object and returns a vector of omega estimates.
}
%%%%%%%%%%%%%%%% VALUE %%%%%%%%%%%%%%%%%%%
\value{
  \itemize{
    \item A vector of productivity estimates - omega.
  }
}
%%%%%%%%%%%%%%%% AUTHOR %%%%%%%%%%%%%%%%%%%
\author{Gabriele Rovigatti}

\name{withinvar}
\alias{withinvar}
%%%%%%%%%%%%%%%% TITLE %%%%%%%%%%%%%%%%%%%
\title{Generate the variance of the demeaned variable}
%%%%%%%%%%%%%%%% DESCRIPTION %%%%%%%%%%%%%%%%%%%
\description{
  \code{withinvar()} subtracts the mean of a vector from the vector itself, and then returns its variance.
}
%%%%%%%%% USAGE %%%%%%%%%%%
  \usage{
    withinvar(inmat)
  }
%%%%%%%%% ARGUMENTS %%%%%%%%%%%
  \arguments{
    %
    \item{inmat}{
      Vector of log(value added output).
    }
    %
  }
%%%%%%%%% DETAILS %%%%%%%%%
  \details{
    \code{withinvar()} accepts a vector as input, then subtracts the mean from it and returns the variance.
  }
%%%%%%%%%%%%%%%% AUTHOR %%%%%%%%%%%%%%%%%%%
\author{Gabriele Rovigatti}

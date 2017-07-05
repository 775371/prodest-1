\name{finalACF}
\alias{finalACF}
%%%%%%%%%%%%%%%% TITLE %%%%%%%%%%%%%%%%%%%
\title{ACF estimation routine}
%%%%%%%%%%%%%%%% DESCRIPTION %%%%%%%%%%%%%%%%%%%
\description{
  \code{finalACF} is the function linking the data cleaning part of the routine with the final function to be boostrapped.
}
%%%%%%%%% USAGE %%%%%%%%%%%
\usage{
  finalACF(ind, data, fnum, snum, cnum, opt, theta0, boot = FALSE)
}
%%%%%%%%% ARGUMENTS %%%%%%%%%%%
\arguments{
  %
  \item{ind}{
    Vector of indices to reshuffle the data.
  }
  %
  \item{data}{
    \code{data.frame} with the data to perform the estimation on.
  }
  %
  \item{fnum}{
    Number of free variables.
  }
  %
  \item{snum}{
    Number of state variables.
  }
  %
  \item{cnum}{
    Number of control variables.
  }
  %
  \item{opt}{
    String with the optimizer.
  }
  %
  \item{theta0}{
    Vector of starting points.
  }
  %
  \item{boot}{
    Binary indicator for the estimation routine being the baseline estimation (\code{boot = FALSE}, the default) or a bootstrap repetition.
  }
  %
}
%%%%%%%%% DETAILS %%%%%%%%%
\details{
  \code{finalACF()} accepts at least 7 inputs: a vector of reshuffled indices, the \code{data.frame} with the data, the number of free, state and control variables, the starting points and the optimizer. It collects the results of \code{gACF()} function - baseline and bootstrapped - calculates the standard errors and stores all in a \code{prod} object.
}
%%%%%%%%%%%%%%%% AUTHOR %%%%%%%%%%%%%%%%%%%
\author{Gabriele Rovigatti}

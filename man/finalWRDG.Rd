\name{finalWRDG}
\alias{finalWRDG}
%%%%%%%%%%%%%%%% TITLE %%%%%%%%%%%%%%%%%%%
\title{Wooldridge - GMM estimation}
%%%%%%%%%%%%%%%% DESCRIPTION %%%%%%%%%%%%%%%%%%%
\description{
  \code{finalWRDG} is the function linking the data cleaning part of the routine with the final function to be boostrapped.
}
%%%%%%%%% USAGE %%%%%%%%%%%
  \usage{
    finalWRDG(ind, data, idvar, timevar, theta0, opt, W)
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
    \item{idvar}{
      Vector of panel identifiers.
    }
    %
    \item{timevar}{
      Vector of time variable.
    }
    %
    \item{theta0}{
      Vector of starting points.
    }
    %
    \item{opt}{
      String with the optimizer.
    }
    %
    \item{W}{
      Weighting matrix.
    }
    %
  }
%%%%%%%%% DETAILS %%%%%%%%%
  \details{
    \code{finalWRDG()} accepts 7 inputs: a vector of reshuffled indices, the \code{data.frame} with the data, the vector of panel and time identifiers, the starting points, the optimizer and the weighting matrix. It collects the results of \code{gWRDG()} function - baseline and bootstrapped - calculates the standard errors and stores all in a \code{prod} object.
  }
%%%%%%%%%%%%%%%% AUTHOR %%%%%%%%%%%%%%%%%%%
\author{Gabriele Rovigatti}

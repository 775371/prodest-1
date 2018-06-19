\name{finalOPLP}
\alias{finalOPLP}
%%%%%%%%%%%%%%%% TITLE %%%%%%%%%%%%%%%%%%%
\title{OP and LP estimation routine}
%%%%%%%%%%%%%%%% DESCRIPTION %%%%%%%%%%%%%%%%%%%
\description{
  \code{finalOPLP} is the function linking the data cleaning part of the routine with the final function to be bootstrapped.
}
%%%%%%%%% USAGE %%%%%%%%%%%
\usage{
  finalOPLP(ind, data, fnum, snum, cnum, opt, theta0, boot, tol, att)
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
    Binary indicator for the estimation routine being the baseline estimation or a bootstrap repetition.
  }
  %
  \item{tol}{
    Optimization tolerance set.
  }
  %
  \item{att}{
    Indicator for attrition in the data - i.e., if  firms exit the market.
  }
  %
}
%%%%%%%%% DETAILS %%%%%%%%%
\details{
  \code{finalOPLP()} accepts at 9 inputs: a vector of reshuffled indices, the \code{data.frame} with the data, the number of free, state and control variables, the starting points, the optimizer, an indicator for bootstrapped repetitions and the optimization tolerance. It collects the results of \code{gACF()} function - baseline and bootstrapped - calculates the standard errors and stores all in a \code{prod} object.
}
%%%%%%%%%%%%%%%% AUTHOR %%%%%%%%%%%%%%%%%%%
\author{Gabriele Rovigatti}

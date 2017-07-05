\name{gWRDG}
\alias{gWRDG}
%%%%%%%%%%%%%%%% TITLE %%%%%%%%%%%%%%%%%%%
\title{Wooldridge - GMM estimation}
%%%%%%%%%%%%%%%% DESCRIPTION %%%%%%%%%%%%%%%%%%%
\description{
  \code{gWRDG} returns the parameters estimates of WRDG models. It is part of the \code{prodestWRDG()} routine.
}
%%%%%%%%% USAGE %%%%%%%%%%%
  \usage{
    gWRDG(theta, Y, X1, X2, Z1, Z2, numR, W = NULL)
  }
%%%%%%%%% ARGUMENTS %%%%%%%%%%%
  \arguments{
    %
    \item{theta}{
      Vector of parameters to be estimated.
    }
    %
    \item{Y}{
      Vector of log(value added output).
    }
    %
    \item{X1}{
      Matrix of regressors of the first equation.
    }
    %
    \item{X2}{
      Matrix of regressors of the second equation.
    }
    %
    \item{Z1}{
      Matrix of instruments of the first equation.
    }
    %
    \item{Z2}{
      Matrix of instruments of the second equation.
    }
    %
    \item{numR}{
      Number of state, free and control variables.
    }
    %
    \item{W}{
      Weighting matrix. By default, \code{W=NULL} in the first stage case and it is a diagonal matrix. In the second stage case, it must be provided (W*).
    }
    %
  }
%%%%%%%%% DETAILS %%%%%%%%%
  \details{
    \code{gWRDG()} estimates the parameters of the Wooldridge model. It accepts at least 7 inputs, generates and optimizes over the group of moment functions E(U_{it}Z^{k}_{it}).
  }
%%%%%%%%%%%%%%%% AUTHOR %%%%%%%%%%%%%%%%%%%
\author{Gabriele Rovigatti}

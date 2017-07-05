\name{gOPLP}
\alias{gOPLP}
%%%%%%%%%%%%%%%% TITLE %%%%%%%%%%%%%%%%%%%
\title{OP and LP Second Stage - GMM estimation}
%%%%%%%%%%%%%%%% DESCRIPTION %%%%%%%%%%%%%%%%%%%
\description{
  \code{gOPLP} returns the second stage parameters estimates of both OP and LP models. It is part of both \code{prodestOP()} and \code{prodestsLP()} routines.
}
%%%%%%%%% USAGE %%%%%%%%%%%
\usage{
  gOPLP(vtheta, mX, mlX, vphi, vlag.phi, vres, stol)
}
%%%%%%%%% ARGUMENTS %%%%%%%%%%%
\arguments{
  %
  \item{vtheta}{
    Vector of parameters to be estimated.
  }
  %
  \item{mX}{
    Matrix of regressors.
  }
  %
  \item{mlX}{
    matrix of lagged regressors.
  }
  %
  \item{vphi}{
    Vector of fitted polynomial.
  }
  %
  \item{vlag.phi}{
    Lagged vector of fitted polynomial.
  }
  %
  \item{vres}{
    Vector of residuals of the free variables.
  }
  %
  \item{stol}{
    Number setting the tolerance of the routine.
  }
  %
}
%%%%%%%%% DETAILS %%%%%%%%%
\details{
  \code{gOPLP()} estimates the second stage of OP and LP routines. It accepts 7 inputs, generates and optimizes over the group of moment functions E(e_{it}X^{k}_{it}).
}
%%%%%%%%%%%%%%%% AUTHOR %%%%%%%%%%%%%%%%%%%
\author{Gabriele Rovigatti}

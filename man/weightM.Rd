\name{weightM}
\alias{weightM}
%%%%%%%%%%%%%%%% TITLE %%%%%%%%%%%%%%%%%%%
\title{Generate optimal GMM weighting matrix}
%%%%%%%%%%%%%%%% DESCRIPTION %%%%%%%%%%%%%%%%%%%
\description{
  In a Woodrdige estimation setting, i.e., in a system GMM framework, this function returns the optimal weighting matrix or the variance-covariance matrix given 1st or 2nd stage estimation results.
}
%%%%%%%%% USAGE %%%%%%%%%%%
\usage{
  weightM(Y, X1, X2, Z1, Z2, betas, numR, SE = FALSE)
}
%%%%%%%%% ARGUMENTS %%%%%%%%%%%
\arguments{
  %
  \item{Y}{
    Vector of log(value added output).
  }
  %
  \item{X1}{
    Matrix of regressors for the first equation.
  }
  %
  \item{X2}{
    Matrix of regressors for the second equation.
  }
  %
  \item{Z1}{
    Matrix of instruments for the first equation.
  }
  %
  \item{Z2}{
    Matrix of instruments for the second equation.
  }
  %
  \item{betas}{
    Vector of first/second stage parameter estimates.
  }
  %
  \item{numR}{
    Number of state + number of free + number of control variables (i.e., number of constrainted parameters).
  }
  %
  \item{SE}{
    Binary indicator for first (\code{SE == FALSE}, the default) or second stage.
  }
  %
}
%%%%%%%%% DETAILS %%%%%%%%%
\details{
  \code{weightM()} accepts at least 7 inputs: Y, X1, X2, Z1, Z2, betas and numR. With these, computes the optimal weighting matrix in a system GMM framework, i.e. W* = sigma*Z'Z. If it is called during the first stage, it returns W*, otherwise will return an estimate of the parameters' standard errors, i.e., the square root of the diagonal of the variance-covariance matrix: 1/N( (X'Z) W* (Z'X) )^{-1}.
}
%%%%%%%%%%%%%%%% AUTHOR %%%%%%%%%%%%%%%%%%%
\author{Gabriele Rovigatti}

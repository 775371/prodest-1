\name{block.boot.resample}
\alias{block.boot.resample}
%%%%%%%%%%%%%%%% TITLE %%%%%%%%%%%%%%%%%%%
\title{Cluster Bootstrap Resampling}
%%%%%%%%%%%%%%%% DESCRIPTION %%%%%%%%%%%%%%%%%%%
\description{
  Function to generate R vectors of resampled IDs. It works reshuffling the row number of the original data - which is stored in the input \code{idvar} along with the relative IDs. The output is a list (N_{i}x1xR), where N_{i} is a random number depending on the reshuffle.
}
%%%%%%%%% USAGE %%%%%%%%%%%
\usage{
  block.boot.resample(idvar, R)
}
%%%%%%%%% ARGUMENTS %%%%%%%%%%%
\arguments{
  %
  \item{idvar}{
    Vector of IDs to be resampled.
  }
  %
  \item{R}{
    Number of samples to be computed.
  }
  %
}
%%%%%%%%% DETAILS %%%%%%%%%
\details{
  \code{block.boot.resample()} accepts two inputs: a vector of IDs - i.e., the vector of panel identifier - and the number of resamplings. For each resampling, it reshuffles the IDs and outputs a vector whose row number is newly-created 'bootstrap' ID, while the value of each cell is the relative row to be reshuffled. This way, each individual can be sampled multiple times, keeping all her number of observations, without generating duplicates.
}
%%%%%%%%%%%%%%%% AUTHOR %%%%%%%%%%%%%%%%%%%
\author{Gabriele Rovigatti}

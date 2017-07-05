\name{lagPanel}
\alias{lagPanel}
%%%%%%%%%%%%%%%% TITLE %%%%%%%%%%%%%%%%%%%
\title{Generate lagged input variables}
%%%%%%%%%%%%%%%% DESCRIPTION %%%%%%%%%%%%%%%%%%%
\description{
  Function to generate lagged variables in a panel.
}
%%%%%%%%% USAGE %%%%%%%%%%%
  \usage{
    lagPanel(idvar, timevar, value)
  }
%%%%%%%%% ARGUMENTS %%%%%%%%%%%
  \arguments{
    %
    \item{idvar}{
      vector of panel identifiers.
    }
    %
    \item{timevar}{
      vector of time identifiers.
    }
    %
     \item{value}{
      variable vector to be lagged.
    }
    %
  }
%%%%%%%%% DETAILS %%%%%%%%%
  \details{
    \code{lagPanel()} accepts three inputs (the ID, the time and the variable to be lagged) and returns the vector of lagged variable. Lagged inputs with no correspondence - i.e., X_{-1} - are returned as NA.
  }
%%%%%%%%%%%%%%%% AUTHOR %%%%%%%%%%%%%%%%%%%
\author{Gabriele Rovigatti}

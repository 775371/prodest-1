\name{printProd}
\alias{printProd}
%%%%%%%%%% TITLE %%%%%%%%%%%%%%%%
\title{
  Print output - prod objects
}
%%%%%%%%%% DESCRIPTION %%%%%%%%%%%%%%%%
\description{
  The \code{printProd()} function accepts a \code{list} of \code{prod} class objects and returns a screen printed tabular in lateX format of the results.
}
%%%%%%%%%% USAGE %%%%%%%%%%%%%%%%
\usage{
  printProd(mods, modnames = NULL, parnames = NULL, outfile = NULL,
            ptime = FALSE, nboot = FALSE, screen = FALSE)
}
%%%%%%%%%% ARGUMENTS %%%%%%%%%%%%%%%%
\arguments{
  \item{mods}{
    a \code{list} of \code{prod} objects.}
  %
  \item{modnames}{
    an optional vector of model names. By default, model names are the \code{@ModelMethod} values in \code{prod} objects.}
  %
  \item{parnames}{
    an optional vector of parameter names. By default, parameter names are the \code{names()} vector of \code{@Estimatespars} in \code{prod} objects.}
  %
  \item{outfile}{
    optional string with the path and directory to store a text file (.txt, .tex, etc. depending on the specified extension) with the tabular. By default \code{outfile = NULL}.}
  %
  \item{ptime}{
    add a row showing the computational time. By default \code{ptime = FALSE}.}
  %
  \item{nboot}{
    add a row showing the number of bootstrap repetitions. By default \code{nboot = FALSE}.}
  %
  \item{screen}{
    print the table on-screen without teX format. By default \code{screen = FALSE}.}
  %
}
%%%%%%%%% VALUE %%%%%%%%%%%
  \value{
The output of the function \code{printProd} is either a screen printed tabular in lateX format of \code{prod} object results or a text file tabular in lateX format of \code{prod} object results.
  }
%%%%%%%%% AUTHOR %%%%%%%%%%%
  \author{
    Gabriele Rovigatti
  }
%%%%%%%%% EXAMPLES %%%%%%%%%%%
  \examples{
    data("chilean")

    # run various models
    WRDGfit <- prodestWRDG_GMM(chilean$Y, fX = cbind(chilean$fX1, chilean$fX2),
                          chilean$sX, chilean$pX, chilean$idvar, chilean$timevar)
    OPfit <- prodestOP(chilean$Y, fX = cbind(chilean$fX1, chilean$fX2), chilean$sX,
                        chilean$pX, chilean$idvar, chilean$timevar)

    # show the output in latex - tabular format
    printProd(list(OPfit, WRDGfit), modnames = c('Olley-Pakes', 'Wooldridge'),
              parnames = c('bunsk', 'bsk', 'bk'))
    # show the output on-screen - no teX format
    printProd(list(OPfit, WRDGfit), modnames = c('Olley-Pakes', 'Wooldridge'),
              parnames = c('bunsk', 'bsk', 'bk'), screen = TRUE)
  }

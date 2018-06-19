############ ROBINSON/WOOLDRIDGE - IV REGRESSION ###############

# function to estimate Wooldridge #
prodestROB <- function(Y, fX, sX, pX, idvar, timevar, cX = NULL){
  Start = Sys.time() # start tracking time
  Y <- checkM(Y) # change all input to matrix
  fX <- checkM(fX)
  sX <- checkM(sX)
  pX <- checkM(pX)
  idvar <- checkM(idvar)
  timevar <- checkM(timevar)
  opt = 'optim' # this is going to be changed once the routine will work with DEoptim and solnp
  snum <- ncol(sX) # find the number of input variables
  fnum <- ncol(fX)
  if (!is.null(cX)) {cX <- checkM(cX); cnum <- ncol(cX)} else {cnum <- 0} # if is there any control, take it into account, else fix the number of controls to 0
  lag.fX = fX # generate fX lags
  for (i in 1:fnum) {
    lag.fX[, i] = lagPanel(fX[, i], idvar = idvar, timevar = timevar)
  }
  polyframe <- data.frame(sX,pX) # vars to be used in polynomial approximation
  regvars <- cbind(model.matrix( ~.^2-1, data = polyframe), sX^2, pX^2) # generate a polynomial of the desired level
  lagregvars <- regvars
  for (i in 1:dim(regvars)[2]) {
    lagregvars[, i] <- lagPanel(idvar = idvar, timevar = timevar, regvars[ ,i])
  }
  data <- model.frame(Y ~ fX + sX + lag.fX + regvars + lagregvars + idvar + timevar) # data.frame of usable observations --> regvars

  iv.out <- ivreg(Y ~ fX + sX + lagregvars | lag.fX + sX + lagregvars, data = data)
  res.names <- c(colnames(fX, do.NULL = FALSE, prefix = 'fX'),
                 colnames(sX, do.NULL = FALSE, prefix = 'sX') ) # generate the list of names for results

  betapar <- iv.out$coefficients[2: (snum + fnum + cnum + 1)]
  betase <- coef(summary(iv.out))[, 2][2: (snum + fnum + cnum + 1)]
  names(betapar) <- res.names # change results' names
  names(betase) <- res.names # change results' names
  elapsed.time = Sys.time() - Start # total running time
  out <- new("prod",
             Model = list(method = 'ROB-IV', FSbetas = NULL, boot.repetitions = NA, elapsed.time = elapsed.time, theta0 = NA,
                          opt = NA, opt.outcome = NULL, nCores = 1),
             Data = list(Y = Y, free = fX, state = sX, proxy = pX, control = cX, idvar = idvar, timevar = timevar,
                         FSresiduals = NULL),
             Estimates = list(pars = betapar, std.errors = betase))
  return(out)
}
## end of prodest Robinson-Wooldridge ##

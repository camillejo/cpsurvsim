#' Inverse CDF for the exponential distribution
#'
#' \code{exp_icdf} returns a value from the inverse CDF of the
#' exponential function.
#'
#' This function uses the exponential function in the form
#' \deqn{f(t)=\theta e^{-\theta t}}
#' to get the inverse CDF
#' \deqn{F^{-1}(u)=(-log(1-u))/\theta.} It can be
#' implemented directly and is also called by the functions
#' \code{\link{exp_memsim}} and \code{\link{exp_cdfsim}}.
#'
#' @param u Numerical value(s) to be converted to exponential variable(s)
#' @param theta Scale parameter \eqn{\theta}
#'
#' @return If inputs are numeric, output is a value or vector of values
#' from the inverse CDF of the exponential function.
#'
#' @examples
#' simdta <- exp_icdf(u = runif(500), theta = 0.05)
#'
#' @export
exp_icdf <- function(u, theta) {
  if(is.numeric(theta) == F | Hmisc::all.is.numeric(u, "test") == F |
     all(u >= 0) == F | theta <= 0) {
    stop("All input values must be numeric and >= 0.")
  }
  return(-log(1 - u) / theta)
}

#' Memoryless simulation for the exponential distribution
#'
#' \code{exp_memsim} returns a dataset simulated from the exponential
#' hazard distribution with K change-points.
#'
#' This function simulates data between change-points from independent
#' exponential distributions using the inverse exponential CDF implemented in
#' the function \code{exp_icdf}.
#'
#' @param theta Scale parameter \eqn{\theta}
#' @param n Sample size
#' @param endtime Maximum study time, point at which all participants
#' are censored
#' @param tau Change-point(s) \eqn{\tau}
#'
#' @return Dataset with n participants with a survival time
#' and censoring indicator.
#'
#' @examples
#' nochangepoint <- exp_memsim(theta = 0.05, n = 10, endtime = 20)
#' onechangepoint <- exp_memsim(theta = c(0.05, 0.01), n = 10,
#'   endtime = 20, tau = 10)
#' twochangepoints <- exp_memsim(theta = c(0.05, 0.01, 0.05),
#'   n = 10, endtime = 20, tau = c(8, 12))
#'
#' @export

exp_memsim <- function(theta, n, endtime, tau = NA) {
  # controls for incorrect input
  if(Hmisc::all.is.numeric(theta, "test") == F |
     is.numeric(n) == F | is.numeric(endtime)==F |
     all(theta > 0) == F | n < 1 | endtime <= 0) {
    stop("All input values must be numeric with value >= 0.")
  }
  n <- as.integer(n)
  id <- NULL
  time <- NULL

  #no change-point
  if(is.na(tau[1]) == T) {
    simdta <- data.frame(id = c(1:n))
    x <- stats::runif(n)
    dt <- cpsurvsim::exp_icdf(u = x, theta = theta)
    simdta$time <- ifelse(dt >= endtime, endtime, dt)
    simdta$censor <- ifelse(simdta$time == endtime, 0, 1)
    dta <- data.frame(time = simdta$time, censor = simdta$censor)
    return(dta)
  }
  #at least one change-point
  if(is.na(tau[1]) == F) {
    if(Hmisc::all.is.numeric(tau, "test") == F | all(tau > 0) == F) {
      stop("Tau must be numeric and > 0.")
    }
    if(endtime < tau[length(tau)]) {
      print("Warning: Change-points occur after endtime.")
    }
    if(length(theta) != (length(tau)+1)) {
      stop("Length of theta and tau not compatible.")
    }
    alltime <- c(0, tau, endtime)
    taudiff <- alltime[2:length(alltime)] - alltime[1:(length(alltime) - 1)]
    s <- n
    nphases <- length(tau) + 1
    phasedta <-list()
    #phase 1
    phasedta[[1]] <- data.frame(id = c(1:n))
    x1 <- stats::runif(s)
    dt <- cpsurvsim::exp_icdf(u = x1, theta = theta[1])
    phasedta[[1]]$time <- ifelse(dt >= taudiff[1], taudiff[1], dt)
    s <- sum(dt >= taudiff[1])
    #other phases
    for(i in 2:nphases) {
      phasedta[[i]] <- subset(phasedta[[i - 1]],
                              phasedta[[i - 1]]$time >= taudiff[i - 1],
                              select = id)
      x <- stats::runif(s)
      p <- cpsurvsim::exp_icdf(u = x, theta = theta[i])
      phasedta[[i]]$time <- ifelse(p >= taudiff[i], taudiff[i], p)
      s <- sum(phasedta[[i]]$time >= taudiff[i])
      colnames(phasedta[[i]]) <-c ("id", paste0("time", i))
    }
    #combine
    simdta <- plyr::join_all(phasedta, by = 'id', type = 'full')
    simdta$survtime <- rowSums(simdta[, -1], na.rm = TRUE)
    simdta$censor <- ifelse(simdta$survtime == endtime, 0, 1)
    dta <- data.frame(time = simdta$survtime, censor = simdta$censor)
    return(dta)
  }
}

#' Inverse CDF simulation for the exponential distribution
#'
#' This function simulates data for the exponential change-point hazard function
#' by simulating values of the exponential function and substituting them into
#' the inverse hazard function. This function allows for up to four
#' change-points.
#'
#' @param theta Scale parameter \eqn{\theta}
#' @param n Sample size
#' @param endtime Maximum study time, point at which all participants
#' are censored
#' @param tau Change-point(s) \eqn{\tau}
#'
#' @return Dataset with n participants with a survival time
#' and censoring indicator.
#'
#' @examples
#' nochangepoint <- exp_cdfsim(theta = 0.05, n = 10, endtime = 20)
#' onechangepoint <- exp_cdfsim(theta = c(0.05, 0.01), n = 10,
#'   endtime = 20, tau = 10)
#' twochangepoints <- exp_cdfsim(theta = c(0.05, 0.01, 0.05),
#'   n = 10, endtime = 20, tau = c(8, 12))
#'
#' @export

exp_cdfsim <- function(theta, n, endtime, tau = NA) {
  # controls for incorrect input
  if(Hmisc::all.is.numeric(theta, "test") == F |
     is.numeric(n) == F | is.numeric(endtime)==F |
     all(theta > 0) == F | n < 1 | endtime <= 0) {
    stop("All input values must be numeric and >= 0.")
  }
  if(length(tau) > 4){
    stop("This function only allows for up to 4 change-points.")
  }
  n <- as.integer(n)
  x <- stats::rexp(n)
  if(is.na(tau[1]) == T) {
    t <- x / theta
  }
  if(is.na(tau[1]) == F) {
    if(Hmisc::all.is.numeric(tau, "test") == F | all(tau > 0) == F) {
      stop("Tau must be numeric and > 0.")
    }
    if(endtime < tau[length(tau)]) {
      print("Warning: Change-points occur after endtime.")
    }
    if(length(theta) != (length(tau)+1)) {
      stop("Length of theta and tau not compatible.")
    }
  }
  if(length(tau) == 1) {
    first <- theta[1] * tau
    cdfcp1 <- function(v) {
      ifelse(v < first, v / theta[1] ,((v - first) / theta[2]) + tau)
    }
    t <- cdfcp1(x)
  }
  if(length(tau) == 2) {
    first <- theta[1] * tau[1]
    second <- first + theta[2] * (tau[2] - tau[1])
    cdfcp2 <- function(v) {
      ifelse(v < first, v / theta[1],
             ifelse(v < second, ((v - first) / theta[2]) + tau[1],
                                ((v - second) / theta[3]) + tau[2]))
    }
    t <- cdfcp2(x)
  }
  if(length(tau) == 3) {
    first <- theta[1] * tau[1]
    second <- first + theta[2] * (tau[2] - tau[1])
    third <- second + theta[3] * (tau[3] - tau[2])
    cdfcp3 <- function(v) {
      ifelse(v < first, v / theta[1],
             ifelse(v < second, ((v - first) / theta[2]) + tau[1],
                    ifelse(v < third, ((v - second) / theta[3]) + tau[2],
                                      ((v - third) / theta[4]) + tau[3])))
    }
    t <- cdfcp3(x)
  }
  if(length(tau) == 4) {
    first <- theta[1] * tau[1]
    second <- first + theta[2] * (tau[2] - tau[1])
    third <- second + theta[3] * (tau[3] - tau[2])
    fourth <- third + theta[4] * (tau[4] - tau[3])
    cdfcp4 <- function(v) {
      ifelse(v < first, v / theta[1],
             ifelse(v < second, ((v - first) / theta[2]) + tau[1],
                    ifelse(v < third, ((v - second) / theta[3]) + tau[2],
                           ifelse(v < fourth, ((v - third) / theta[4]) + tau[3],
                                              ((v - fourth) / theta[5]) + tau[4]))))
    }
    t <- cdfcp4(x)
  }
  C <- rep(endtime, length(x)) #all censored at endtime
  time <- pmin(t, C)  #observed time is min of censored and true
  censor <- as.numeric(time != endtime) #if not endtime then dropout
  dta <- data.frame(time = time, censor = censor)
  return(dta)
}

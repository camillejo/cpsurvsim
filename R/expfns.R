#' Inverse CDF for the exponential distribution
#'
#' \code{exp_icdf} simulates values from the inverse CDF of the
#' exponential distribution.
#'
#' This function uses the exponential distribution of the form
#' \deqn{f(t)=\theta exp(-\theta t)}
#' to get the inverse CDF
#' \deqn{F^(-1)(u)=(-log(1-u))/\theta} where \eqn{u}
#' is a uniform random variable. It can be
#' implemented directly and is also called by the function
#' \code{\link{exp_memsim}}.
#'
#' @param n Number of output exponential values
#' @param theta Scale parameter \eqn{\theta}
#'
#' @return Output is a value or a vector of values
#' from the exponential distribution.
#'
#' @examples
#' simdta <- exp_icdf(n = 10, theta = 0.05)
#'
#' @export
exp_icdf <- function(n, theta) {
  if(is.numeric(n) == FALSE | n <= 0) {
    stop("n must be an integer greater than 0.")
  }
  if(is.numeric(theta) == FALSE | theta <= 0) {
    stop("Theta must be numeric and > 0.")
  }
  x <- stats::rexp(n)
  return(x / theta)
}

#' Memoryless simulation for the exponential change-point hazard distribution
#'
#' \code{exp_memsim} simulates time-to-event data from the exponential change-point
#' hazard distribution by implementing the memoryless method.
#'
#' This function simulates time-to-event data between \eqn{K} change-points from
#' independent exponential distributions using the inverse CDF implemented
#' in \code{exp_icdf}. This method applies Type I right censoring at the endtime
#' specified by the user.
#'
#' @param n Sample size
#' @param endtime Maximum study time, point at which all participants
#' are censored
#' @param theta Scale parameter \eqn{\theta}
#' @param tau Change-point(s) \eqn{\tau}
#'
#' @return Dataset with n participants including a survival time
#' and censoring indicator (0 = censored, 1 = event).
#'
#' @examples
#' nochangepoint <- exp_memsim( n = 10, endtime = 20, theta = 0.05)
#' onechangepoint <- exp_memsim(n = 10, endtime = 20,
#'   theta = c(0.05, 0.01), tau = 10)
#' twochangepoints <- exp_memsim(n = 10, endtime = 20,
#'   theta = c(0.05, 0.01, 0.05), tau = c(8, 12))
#'
#' # Pay attention to how you parameterize your model!
#' # This simulates a decreasing hazard
#' set.seed(1245)
#' decreasingHazard <- exp_memsim(n = 10, endtime = 20,
#'   theta = c(0.5, 0.2, 0.01), tau = c(8, 12))
#' # This tries to fit an increasing hazard, resulting in biased estimates
#' cp2.nll <- function(par, tau = tau, dta = dta){
#'   theta1 <- par[1]
#'   theta2 <- par[2]
#'   theta3 <- par[3]
#'   ll <- log(theta1) * sum(dta$time < tau[1])+
#'         log(theta2) * sum((tau[1] <= dta$time) * (dta$time < tau[2])) +
#'         log(theta3) * sum((dta$time >= tau[2]) * dta$censor) -
#'         theta1 * sum(dta$time * (dta$time < tau[1]) +
#'           tau[1] * (dta$time >= tau[1])) -
#'         theta2 * sum((dta$time - tau[1]) * (dta$time >= tau[1]) *
#'           (dta$time < tau[2]) + (tau[2] - tau[1]) * (dta$time >= tau[2])) -
#'         theta3 * sum((dta$time - tau[2]) * (dta$time >= tau[2]))
#'   return(-ll)
#' }
#' optim(par = c(0.001, 0.1, 0.5), fn = cp2.nll,
#'       tau = c(8, 12), dta = decreasingHazard)
#'
#' @export

exp_memsim <- function(n, endtime, theta, tau = NA) {
  # controls for incorrect input
  if(is.numeric(n) == FALSE | n <= 0) {
    stop("n must be an integer greater than 0.")
  }
  if(Hmisc::all.is.numeric(theta, "test") == FALSE | all(theta > 0) == FALSE |
     is.numeric(endtime) == FALSE | endtime <= 0) {
    stop("Endtime and theta must be numeric and > 0.")
  }
  n <- as.integer(n)
  id <- NULL
  time <- NULL

  #no change-point
  if(is.na(tau[1]) == TRUE) {
    simdta <- data.frame(id = c(1:n))
    dt <- cpsurvsim::exp_icdf(n = n, theta = theta)
    simdta$time <- ifelse(dt >= endtime, endtime, dt)
    simdta$censor <- ifelse(simdta$time == endtime, 0, 1)
    dta <- data.frame(time = simdta$time, censor = simdta$censor)
    return(dta)
  }
  #at least one change-point
  if(is.na(tau[1]) == FALSE) {
    if(Hmisc::all.is.numeric(tau, "test") == FALSE | all(tau > 0) == FALSE) {
      stop("Tau must be numeric and > 0.")
    }
    if(endtime < tau[length(tau)]) {
      warning("Warning: Change-points occur after endtime.")
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
    dt <- cpsurvsim::exp_icdf(n = s, theta = theta[1])
    phasedta[[1]]$time <- ifelse(dt >= taudiff[1], taudiff[1], dt)
    s <- sum(dt >= taudiff[1])
    #other phases
    for(i in 2:nphases) {
      if(s == 0){
        break
      }
      phasedta[[i]] <- subset(phasedta[[i - 1]],
                              phasedta[[i - 1]]$time >= taudiff[i - 1],
                              select = id)
      p <- cpsurvsim::exp_icdf(n = s, theta = theta[i])
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

#' Inverse CDF simulation for the exponential change-point hazard distribution
#'
#' \code{exp_cdfsim} simulates time-to-event data from the exponential change-point
#' hazard distribution by implementing the inverse CDF method.
#'
#' This function simulates data for the exponential change-point hazard
#' distribution with \eqn{K} change-points by simulating values of the exponential
#' distribution and substituting them into the inverse hazard function. This
#' method applies Type I right censoring at the endtime specified by the user.
#' This function allows for up to four change-points.
#'
#' @param n Sample size
#' @param endtime Maximum study time, point at which all participants
#' are censored
#' @param theta Scale parameter \eqn{\theta}
#' @param tau Change-point(s) \eqn{\tau}
#'
#' @return Dataset with n participants including a survival time
#' and censoring indicator (0 = censored, 1 = event).
#'
#' @examples
#' nochangepoint <- exp_cdfsim(n = 10, endtime = 20, theta = 0.05)
#' onechangepoint <- exp_cdfsim(n = 10, endtime = 20,
#'   theta = c(0.05, 0.01), tau = 10)
#' twochangepoints <- exp_cdfsim(n = 10, endtime = 20,
#'   theta = c(0.05, 0.01, 0.05), tau = c(8, 12))
#'
#' # Pay attention to how you parameterize your model!
#' # This simulates a decreasing hazard
#' set.seed(7830)
#' decreasingHazard <- exp_cdfsim(n = 10, endtime = 20,
#'   theta = c(0.5, 0.2, 0.01), tau = c(8, 12))
#' # This tries to fit an increasing hazard, resulting in biased estimates
#' cp2.nll <- function(par, tau = tau, dta = dta){
#'   theta1 <- par[1]
#'   theta2 <- par[2]
#'   theta3 <- par[3]
#'   ll <- log(theta1) * sum(dta$time < tau[1])+
#'         log(theta2) * sum((tau[1] <= dta$time) * (dta$time < tau[2])) +
#'         log(theta3) * sum((dta$time >= tau[2]) * dta$censor) -
#'         theta1 * sum(dta$time * (dta$time < tau[1]) +
#'           tau[1] * (dta$time >= tau[1])) -
#'         theta2 * sum((dta$time - tau[1]) * (dta$time >= tau[1]) *
#'           (dta$time < tau[2]) + (tau[2] - tau[1]) * (dta$time >= tau[2])) -
#'         theta3 * sum((dta$time - tau[2]) * (dta$time >= tau[2]))
#'   return(-ll)
#' }
#' optim(par = c(0.001, 0.1, 0.5), fn = cp2.nll,
#'       tau = c(8, 12), dta = decreasingHazard)
#' @export

exp_cdfsim <- function(n, endtime, theta, tau = NA) {
  # controls for incorrect input
  if(is.numeric(n) == FALSE | n <= 0) {
    stop("n must be an integer greater than 0.")
  }
  if(Hmisc::all.is.numeric(theta, "test") == FALSE | all(theta > 0) == FALSE |
     is.numeric(endtime) == FALSE | endtime <= 0) {
    stop("Endtime and theta must be numeric and > 0.")
  }
  if(length(tau) > 4){
    stop("This function only allows for up to 4 change-points.")
  }
  n <- as.integer(n)
  x <- stats::rexp(n)
  if(is.na(tau[1]) == TRUE) {
    t <- x / theta
  }
  if(is.na(tau[1]) == FALSE) {
    if(Hmisc::all.is.numeric(tau, "test") == FALSE | all(tau > 0) == FALSE) {
      stop("Tau must be numeric and > 0.")
    }
    if(endtime < tau[length(tau)]) {
      warning("Warning: Change-points occur after endtime.")
    }
    if(length(theta) != (length(tau)+1)) {
      stop("Length of theta and tau not compatible.")
    }
  }
  if(length(tau) == 1 & is.na(tau[1]) == FALSE) {
    first <- theta[1] * tau
    cdfcp1 <- function(v) {
      ifelse(v < first, v / theta[1] ,((v - first) / theta[2]) + tau)
    }
    t <- cdfcp1(x)
  }
  if(length(tau) == 2 & is.na(tau[1]) == FALSE) {
    first <- theta[1] * tau[1]
    second <- first + theta[2] * (tau[2] - tau[1])
    cdfcp2 <- function(v) {
      ifelse(v < first, v / theta[1],
             ifelse(v < second, ((v - first) / theta[2]) + tau[1],
                                ((v - second) / theta[3]) + tau[2]))
    }
    t <- cdfcp2(x)
  }
  if(length(tau) == 3 & is.na(tau[1]) == FALSE) {
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
  if(length(tau) == 4 & is.na(tau[1]) == FALSE) {
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
  endtime <- as.numeric(endtime)
  C <- rep(endtime, length(x)) #all censored at endtime
  time <- pmin(t, C)  #observed time is min of censored and true
  censor <- as.numeric(time != endtime) #if not endtime then dropout
  dta <- data.frame(time = time, censor = censor)
  return(dta)
}

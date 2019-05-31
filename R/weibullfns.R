#' Inverse CDF for the Weibull distribution
#'
#' \code{weib_icdf} returns a value from the inverse CDF of the
#' Weibull distribution.
#'
#' This function uses the Weibull density of the form
#' \deqn{f(t)=\theta t^(\gamma - 1)exp(-\theta/\gamma t^(\gamma))}
#' to get the inverse CDF
#' \deqn{F^(-1)(u)=(-\gamma/\theta log(1-u))^(1/\gamma).} It can be
#' implemented directly and is also called by the functions
#' \code{\link{weib_memsim}} and \code{\link{weib_cdfsim}}.
#'
#' @param u Numerical value(s) to be converted to Weibull variable(s)
#' @param theta Scale parameter \eqn{\theta}
#' @param gamma Shape parameter \eqn{\gamma}
#'
#' @return Output is a value or vector of values
#' from the inverse CDF of the Weibull distribution.
#'
#' @examples
#' simdta <- weib_icdf(u = runif(10), theta = 0.05, gamma = 2)
#'
#' @export

weib_icdf <- function(u, theta, gamma) {
  if(is.numeric(theta) == FALSE | is.numeric(gamma) == FALSE |
     Hmisc::all.is.numeric(u, "test") == FALSE | theta <= 0 |
     gamma <= 0 | all(u >= 0) == FALSE) {
    stop("All input values must be numeric and >= 0.")
  }
  return(((-gamma / theta) * log(1 - u)) ^ (1 / gamma))
}


#' Memoryless simulation for the Weibull change-point hazard distribution
#'
#' \code{weib_memsim} simulates time-to-event data from the Weibull change-point
#' hazard distribution by implementing the memoryless method.
#'
#' This function simulates time-to-event data between \eqn{K} change-points \eqn{\tau}
#' from independent Weibull distributions using the inverse Weibull CDF
#' implemented in \code{\link{weib_icdf}}. This method applies Type I right
#' censoring at the endtime specified by the user. \eqn{\gamma} is
#' held constant.
#'
#' @param theta Scale parameter \eqn{\theta}
#' @param gamma Shape parameter \eqn{\gamma}
#' @param n Sample size
#' @param endtime Maximum study time, point at which all participants
#' are censored
#' @param tau Change-point(s) \eqn{\tau}
#'
#' @return Dataset with n participants including a survival time
#' and censoring indicator (0 = censored, 1 = event).
#'
#' @examples
#' nochangepoint <- weib_memsim(theta = 0.05, gamma = 2, n = 10,
#'   endtime = 20)
#' onechangepoint <- weib_memsim(theta = c(0.05, 0.01), gamma = 2,
#'   n = 10, endtime = 20, tau = 10)
#' twochangepoints <- weib_memsim(theta = c(0.05, 0.01, 0.05),
#'   gamma = 2, n = 10, endtime = 20, tau = c(8, 12))
#'
#' @export

weib_memsim <- function(theta, gamma, n, endtime, tau = NA) {
  # controls for incorrect input
  if(Hmisc::all.is.numeric(theta, "test") == FALSE | is.numeric(gamma) == FALSE |
     is.numeric(n) == FALSE | is.numeric(endtime)==FALSE |
     all(theta > 0) == FALSE | gamma <= 0 | n < 1 | endtime <= 0) {
    stop("All input values must be numeric and >= 0.")
  }
  n <- as.integer(n)
  id <- NULL
  time <- NULL

  #no change-point
  if(is.na(tau[1]) == TRUE) {
    simdta <- data.frame(id = c(1:n))
    x <- stats::runif(n)
    dt <- cpsurvsim::weib_icdf(u = x, theta = theta, gamma = gamma)
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
    x1 <- stats::runif(s)
    dt <- cpsurvsim::weib_icdf(u = x1, theta = theta[1], gamma = gamma)
    phasedta[[1]]$time <- ifelse(dt >= taudiff[1], taudiff[1], dt)
    s <- sum(dt >= taudiff[1])
    #other phases
    for(i in 2:nphases) {
      phasedta[[i]] <- subset(phasedta[[i - 1]],
                              phasedta[[i - 1]]$time >= taudiff[i - 1],
                              select = id)
      x <- stats::runif(s)
      p <- cpsurvsim::weib_icdf(u = x, theta = theta[i], gamma = gamma)
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

#' Inverse CDF simulation for the Weibull change-point hazard distribution
#'
#' \code{weib_cdfsim} simulates time-to-event data from the Weibull change-point
#' hazard distribution by implementing the inverse CDF method.
#'
#' This function simulates data from the Weibull change-point hazard distribution
#' with \eqn{K} change-points by simulating values of the exponential distribution and
#' substituting them into the inverse hazard function. This method applies Type I
#' right censoring at the endtime specified by the user. This function allows for
#' up to four change-points and \eqn{\gamma} is held constant.
#'
#' @param theta Scale parameter \eqn{\theta}
#' @param gamma Shape parameter \eqn{\gamma}
#' @param n Sample size
#' @param endtime Maximum study time, point at which all participants
#' are censored
#' @param tau Change-point(s) \eqn{\tau}
#'
#' @return Dataset with n participants including a survival time
#' and censoring indicator (0 = censored, 1 = event).
#'
#' @examples
#' nochangepoint <- weib_cdfsim(theta = 0.5, gamma = 2, n = 10,
#'   endtime = 20)
#' onechangepoint <- weib_cdfsim(theta = c(0.05, 0.01), gamma = 2,
#'   n = 10, endtime = 20, tau = 10)
#' twochangepoints <- weib_cdfsim(theta = c(0.05, 0.01, 0.05),
#'   gamma = 2, n = 10, endtime = 20, tau = c(8, 12))
#'
#' @export

weib_cdfsim <- function(theta, gamma, n, endtime, tau = NA) {
  # controls for incorrect input
  if(Hmisc::all.is.numeric(theta, "test") == FALSE | is.numeric(gamma) == FALSE |
     is.numeric(n) == FALSE | is.numeric(endtime)==FALSE |
     all(theta > 0) == FALSE | gamma <= 0 | n < 1 | endtime <= 0) {
    stop("All input values must be numeric and > 0.")
  }
  if(length(tau) > 4){
    stop("This function only allows for up to 4 change-points.")
  }
  n <- as.integer(n)
  x <- stats::rexp(n)
  if(is.na(tau[1]) == TRUE) {
    t <- cpsurvsim::weib_icdf(x, theta = theta, gamma = gamma)
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
  if(length(tau) == 1) {
    first <- (theta / gamma) * tau ^ gamma
    cdfcp1 <- function(v) {
      ifelse(v < first, ((gamma / theta[1]) * v) ^ (1 / gamma),
             ((gamma / theta[2]) * (v - first) + tau ^ gamma) ^ (1 / gamma))
    }
    t <- cdfcp1(x)
  }

  if(length(tau) == 2) {
    first <- (theta[1] / gamma) * tau[1] ^ gamma #set interval 1
    second <- first + (theta[2] / gamma) * (tau[2] ^ gamma - tau[1] ^ gamma) #set interval 2
    cdfcp2 <- function(v) {
      ifelse(v < first, ((gamma / theta[1]) * v) ^ (1 / gamma),
             ifelse(v < second,
                    ((gamma / theta[2]) * (v - first) + tau[1] ^ gamma) ^ (1 / gamma),
                    ((gamma / theta[3]) * (v - second) + tau[2] ^ gamma) ^ (1 / gamma)))
    }
    t <- cdfcp2(x)
  }

  if(length(tau) == 3) {
    first <- (theta[1] / gamma) * tau[1] ^ gamma
    second <- first + (theta[2] / gamma) * (tau[2] ^ gamma - tau[1] ^ gamma)
    third <- second + (theta[3] / gamma) * (tau[3] ^ gamma - tau[2] ^ gamma)
    cdfcp3 <- function(v) {
      ifelse(v < first, ((gamma / theta[1]) * v) ^ (1 / gamma),
             ifelse(v < second,
                    ((gamma / theta[2]) * (v - first) + tau[1] ^ gamma) ^ (1 / gamma),
                    ifelse(v < third,
                           ((gamma / theta[3]) * (v - second) + tau[2] ^ gamma) ^ (1 / gamma),
                           ((gamma / theta[4]) * (v - third) + tau[3] ^ gamma) ^ (1 / gamma))))
    }
    t <- cdfcp3(x)
  }

  if(length(tau) == 4) {
    first <- (theta[1] / gamma) * tau[1] ^ gamma
    second <- first + (theta[2] / gamma) * (tau[2] ^ gamma - tau[1] ^ gamma)
    third <- second + (theta[3] / gamma) * (tau[3] ^ gamma - tau[2] ^ gamma)
    fourth <- third + (theta[4] / gamma) * (tau[4] ^ gamma - tau[3] ^ gamma)
    cdfcp4 <- function(v) {
      ifelse(v < first, ((gamma / theta[1]) * v) ^ (1 / gamma),
             ifelse(v < second,
                    ((gamma / theta[2]) * (v - first) + tau[1] ^ gamma) ^ (1 / gamma),
                    ifelse(v < third,
                           ((gamma / theta[3]) * (v - second) + tau[2] ^ gamma) ^ (1 / gamma),
                           ifelse(v < fourth,
                                  ((gamma / theta[4]) * (v-third) + tau[3] ^ gamma) ^ (1 / gamma),
                                  ((gamma / theta[5]) * (v-fourth) + tau[4] ^ gamma) ^ (1 / gamma)))))
    }
    t <- cdfcp4(x)
  }

  C <- rep(endtime, length(x)) #all censored at endtime
  time <- pmin(t, C)  #observed time is min of censored and true
  censor <- as.numeric(time != endtime) #if not endtime then dropout
  dta <- data.frame(time = time, censor = censor)
  return(dta)
}

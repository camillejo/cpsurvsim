#' Inverse CDF value generation for the Weibull distribution
#'
#' \code{weib_icdf} returns a value from the Weibull distribution by
#' using the inverse CDF.
#'
#' This function uses the Weibull density of the form
#' \deqn{f(t)=\theta t^(\gamma - 1)exp(-\theta/\gamma t^(\gamma))}
#' to get the inverse CDF
#' \deqn{F^(-1)(u)=(-\gamma/\theta log(1-u))^(1/\gamma)} where \eqn{u}
#' is a uniform random variable. It can be implemented directly and is
#' also called by the function \code{\link{weib_memsim}}.
#'
#' @param n Number of output Weibull values
#' @param theta Scale parameter \eqn{\theta}
#' @param gamma Shape parameter \eqn{\gamma}
#'
#' @return Output is a value or vector of values
#' from the Weibull distribution.
#'
#' @examples
#' simdta <- weib_icdf(n = 10, theta = 0.05, gamma = 2)
#'
#' @export

weib_icdf <- function(n, gamma, theta) {
  if(is.numeric(n) == FALSE | n <= 0) {
    stop("n must be an integer greater than 0.")
  }
  if(is.numeric(theta) == FALSE | is.numeric(gamma) == FALSE |
     theta <= 0 | gamma <= 0) {
    stop("Theta and gamma must be numeric and > 0.")
  }
  x <- stats::rexp(n)
  return(((gamma / theta) * x) ^ (1 / gamma))
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
#' @param n Sample size
#' @param endtime Maximum study time, point at which all participants
#' are censored
#' @param gamma Shape parameter \eqn{\gamma}
#' @param theta Scale parameter \eqn{\theta}
#' @param tau Change-point(s) \eqn{\tau}
#'
#' @return Dataset with n participants including a survival time
#' and censoring indicator (0 = censored, 1 = event).
#'
#' @examples
#' nochangepoint <- weib_memsim(n = 10, endtime = 20, gamma = 2,
#'   theta = 0.05)
#' onechangepoint <- weib_memsim(n = 10, endtime = 20, gamma = 2,
#'   theta = c(0.05, 0.01), tau = 10)
#' twochangepoints <- weib_memsim(n = 10, endtime = 20, gamma = 2,
#'   theta = c(0.05, 0.01, 0.05), tau = c(8, 12))
#'
#' # Pay attention to how you parameterize your model!
#' # This simulates an increasing hazard
#' set.seed(5738)
#' increasingHazard <- weib_memsim(n = 100, endtime = 20, gamma = 2,
#'   theta = c(0.001, 0.005, 0.02), tau = c(8, 12))
#' # This tries to fit a decreasing hazard, resulting in biased estimates
#' cp2.nll <- function(par, tau = tau, gamma = gamma, dta = dta){
#'   theta1 <- par[1]
#'   theta2 <- par[2]
#'   theta3 <- par[3]
#'   ll <- (gamma - 1) * sum(dta$censor * log(dta$time)) +
#'     log(theta1) * sum((dta$time < tau[1])) +
#'     log(theta2) * sum((tau[1] <= dta$time) * (dta$time < tau[2])) +
#'     log(theta3) * sum((dta$time >= tau[2]) * dta$censor) -
#'     (theta1/gamma) * sum((dta$time^gamma) * (dta$time < tau[1]) +
#'       (tau[1]^gamma) * (dta$time >= tau[1])) -
#'     (theta2/gamma) * sum((dta$time^gamma - tau[1]^gamma) *
#'       (dta$time >= tau[1]) * (dta$time<tau[2]) +
#'       (tau[2]^gamma - tau[1]^gamma) * (dta$time >= tau[2])) -
#'     (theta3/gamma) * sum((dta$time^gamma - tau[2]^gamma) *
#'       (dta$time >= tau[2]))
#'   return(-ll)
#' }
#' optim(par = c(0.2, 0.02, 0.01), fn = cp2.nll,
#'   tau = c(8, 12), gamma = 2,
#'   dta = increasingHazard)
#'
#' @export

weib_memsim <- function(n, endtime, gamma, theta, tau = NA) {
  # controls for incorrect input
  if(is.numeric(n) == FALSE | n <= 0) {
    stop("n must be an integer greater than 0.")
  }
  if(Hmisc::all.is.numeric(theta, "test") == FALSE | is.numeric(gamma) == FALSE |
     all(theta > 0) == FALSE | gamma <= 0 | is.numeric(endtime) == FALSE |
     endtime <= 0) {
    stop("Endtime, theta and gamma must be numeric and > 0.")
  }
  n <- as.integer(n)
  id <- NULL
  time <- NULL

  #no change-point
  if(is.na(tau[1]) == TRUE) {
    simdta <- data.frame(id = c(1:n))
    dt <- cpsurvsim::weib_icdf(n = n, gamma = gamma, theta = theta)
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
    dt <- cpsurvsim::weib_icdf(n = as.numeric(s), gamma = gamma, theta = theta[1])
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
      p <- cpsurvsim::weib_icdf(n = as.numeric(s), gamma = gamma, theta = theta[i])
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
#' @param n Sample size
#' @param endtime Maximum study time, point at which all participants
#' are censored
#' @param gamma Shape parameter \eqn{\gamma}
#' @param theta Scale parameter \eqn{\theta}
#' @param tau Change-point(s) \eqn{\tau}
#'
#' @return Dataset with n participants including a survival time
#' and censoring indicator (0 = censored, 1 = event).
#'
#' @examples
#' nochangepoint <- weib_cdfsim(n = 10, endtime = 20, gamma = 2,
#'   theta = 0.5)
#' onechangepoint <- weib_cdfsim(n = 10, endtime = 20, gamma = 2,
#'   theta = c(0.05, 0.01), tau = 10)
#' twochangepoints <- weib_cdfsim(n = 10, endtime = 20, gamma = 2,
#'   theta = c(0.05, 0.01, 0.05), tau = c(8, 12))
#'
#'#' # Pay attention to how you parameterize your model!
#' # This simulates an increasing hazard
#' set.seed(9945)
#' increasingHazard <- weib_cdfsim(n = 100, endtime = 20, gamma = 2,
#'   theta = c(0.001, 0.005, 0.02), tau = c(8, 12))
#' # This tries to fit a decreasing hazard, resulting in biased estimates
#' cp2.nll <- function(par, tau = tau, gamma = gamma, dta = dta){
#'   theta1 <- par[1]
#'   theta2 <- par[2]
#'   theta3 <- par[3]
#'   ll <- (gamma - 1) * sum(dta$censor * log(dta$time)) +
#'     log(theta1) * sum((dta$time < tau[1])) +
#'     log(theta2) * sum((tau[1] <= dta$time) * (dta$time < tau[2])) +
#'     log(theta3) * sum((dta$time >= tau[2]) * dta$censor) -
#'     (theta1/gamma) * sum((dta$time^gamma) * (dta$time < tau[1]) +
#'       (tau[1]^gamma) * (dta$time >= tau[1])) -
#'     (theta2/gamma) * sum((dta$time^gamma - tau[1]^gamma) *
#'       (dta$time >= tau[1]) * (dta$time<tau[2]) +
#'       (tau[2]^gamma - tau[1]^gamma) * (dta$time >= tau[2])) -
#'     (theta3/gamma) * sum((dta$time^gamma - tau[2]^gamma) *
#'       (dta$time >= tau[2]))
#'   return(-ll)
#' }
#' optim(par = c(0.2, 0.02, 0.01), fn = cp2.nll,
#'   tau = c(8, 12), gamma = 2,
#'   dta = increasingHazard)
#'
#' @export

weib_cdfsim <- function(n, endtime, gamma, theta, tau = NA) {
  # controls for incorrect input
  if(is.numeric(n) == FALSE | n <= 0) {
    stop("n must be an integer greater than 0.")
  }
  if(Hmisc::all.is.numeric(theta, "test") == FALSE | is.numeric(gamma) == FALSE |
     all(theta > 0) == FALSE | gamma <= 0 | is.numeric(endtime) == FALSE |
     endtime <= 0) {
    stop("Endtime, theta and gamma must be numeric and > 0.")
  }
  if(length(tau) > 4){
    stop("This function only allows for up to 4 change-points.")
  }
  n <- as.numeric(n)
  x <- stats::rexp(n)
  if(is.na(tau[1]) == TRUE) {
    t <- ((gamma / theta) * x) ^ (1 / gamma)
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
    first <- (theta / gamma) * tau ^ gamma
    cdfcp1 <- function(v) {
      ifelse(v < first, ((gamma / theta[1]) * v) ^ (1 / gamma),
             ((gamma / theta[2]) * (v - first) + tau ^ gamma) ^ (1 / gamma))
    }
    t <- cdfcp1(x)
  }

  if(length(tau) == 2 & is.na(tau[1]) == FALSE) {
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

  if(length(tau) == 3 & is.na(tau[1]) == FALSE) {
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

  if(length(tau) == 4 & is.na(tau[1]) == FALSE) {
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
  endtime <- as.numeric(endtime)
  C <- rep(endtime, length(x)) #all censored at endtime
  time <- pmin(t, C)  #observed time is min of censored and true
  censor <- as.numeric(time != endtime) #if not endtime then dropout
  dta <- data.frame(time = time, censor = censor)
  return(dta)
}

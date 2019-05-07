#' cpsurvsim: Simulating Survival Data from Change-Point Hazard Distributions
#'
#' The cpsurvsim package simulates time-to-event data
#' with type I right censoring using two methods: the inverse CDF
#' method and a memoryless method (for more information on simulation
#' methods, see the vignette). We include two parametric
#' distributions: exponential and Weibull.
#'
#' @section cpsurvsim functions:
#' For the exponential distribution, the \code{\link{exp_icdf}}
#' function simulates values from the inverse exponential distribution.
#' \code{\link{exp_cdfsim}} and \code{\link{exp_memsim}} return
#' time-to-event datasets simulated using the inverse CDF and memoryless
#' methods respectively.
#'
#' For the Weibull distribution, the \code{\link{weib_icdf}} function
#' simulates values from the inverse Weibull distribution.
#' \code{\link{weib_cdfsim}} and \code{\link{weib_memsim}} return
#' time-to-event datasets simulated using the inverse CDF and memoryless
#' methods respectively.
#'
#' @docType package
#' @name cpsurvsim
#'
#' @importFrom Hmisc all.is.numeric
#' @importFrom plyr join_all
#' @importFrom stats rexp runif
#' @import knitr
#' @import rmarkdown
NULL

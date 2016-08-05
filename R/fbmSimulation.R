#' Simulate a discrete fractional Brownian motion
#'
#' \code{fbmSim} simulates a fractional Brownian motion using the Levinson's
#' method.
#'
#' @param n	Length of the time series.
#' @param H Hurst exponent.
#' @return A numeric vector containing the simulated time series.
#' @examples
#' # simulate a dfbm process with H=0.2 and plot it
#' z = fbmSim(100, 0.2)
#' ts.plot(z)
#' @note The code of the \code{fbmSim} is largelly based on the implementation
#'  provided in the Coeurjolly article (see references).
#' @references Jean-Francois, Coeurjolly. "Simulation and Identification of the
#' Fractional Brownian Motion: A bibliographical and comparative study." (2007).
#' @import Rcpp
#' @export
fbmSim = function(n, H) {
  .Call('fracdet_simulateFbmCpp', PACKAGE = 'fracdet', n, H)
}

#' Simulate a discrete fractional Gaussian noise
#' @param n	Length of the time series.
#' @param H Hurst exponent.
#' @return A numeric vector containing the simulated time series.
#' @seealso \code{\link[FGN]{SimulateFGN}}
#' @details The FFT is used during the simulation so it is most efficient if
#' n is a power of 2.
#' @note The \code{fgnSim} is just a wrapper of the \code{SimulateFGN} function
#' from the \code{FGN} package.
#' @author A.I. McLeod
#' @examples
#' #simulate a dfgn process with H=0.2 and plot it
#' z = fgnSim(100, 0.2)
#' ts.plot(z)
#' @importFrom FGN SimulateFGN
#' @export
fgnSim = function(n, H) {
  FGN::SimulateFGN(n,H)
}

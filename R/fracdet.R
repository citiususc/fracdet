# fracdet selective imports -----------------------------------------------
#' @importFrom graphics lines plot points
#' @importFrom stats coef lm mad nls predict sd var  vcov
NULL

# fracdet constructor -----------------------------------------------------
checkFracdetArgs = function(x, fbmPars) {
  if (!inherits(x, "wd")) {
    stop("x is not a \"wd\" object")
  }
  if (!inherits(fbmPars, "lm") && !inherits(fbmPars, "nls")) {
    stop("fbmPars shoulb be either a \"lm\" or a \"nls\" object")
  }
}

#' Fractal-deterministic model
#'
#' The \code{fracdet} class represents a wavelet transform that is assumed to
#' be modelled  by a simple fractal-deterministic model with long-memory
#' properties. The model consists on a  linear superposition of a fractional
#' brownian motion (fBm) \eqn{B} and a deterministic band-limited signal
#' \eqn{x}. That is, the observed signal is \eqn{Y = x + B}. The \code{fracdet}
#' receives as parameters the wavelet transform of \eqn{Y} and an R fitted model
#' object (\code{nls} or \code{lm}) containing the estimation of the parameters
#' characterizing the fBm (The Hurst exponent \eqn{H} and the "dispersion"
#' parameter \eqn{\sigma ^ 2}{sigma ^ 2}). The main method of the \code{fracdet} class
#' is \code{estimateDetSignal}, which permits the estimation of the deterministic
#' signal \eqn{x} using bayesian modelling techniques in the wavelet-domain.
#'
#' In addition to the specific methods for the \code{fracdet} class, all the
#' \code{wd} methods can be used with a \code{fracdet} object.
#'
#' @param x A \code{wd} object representing the wavelet transform of the observed
#' signal \eqn{Y}.
#' @param fbmPars Either a \code{nls} or a \code{lm} object (see
#' \code{\link[stats]{nls}} and \code{\link[stats]{nls}}) providing estimates
#' of the parameters characterizing the fBm part of the model (The Hurst
#' exponent \eqn{H} and the "dispersion" parameter \eqn{\sigma ^ 2}{sigma ^ 2}).
#' Thus, the names of the fitted variables are assumed to be \code{H} and
#' \code{sigma2}.
#' @return A S3 \code{fracdet} object that represents the observed \eqn{Y} signal
#' in the wavelet domain. The \code{fracdet} object will also contain the
#' estimates of the fBm parameters.
#' @note When computing the wavelet transform of the \eqn{Y} signal through the
#' \code{wd} method (see \code{\link[wavethresh]{wd}}) we recommend using the
#' "symmetric" boundary handling method (\code{bc == "symmetric"}) if
#' \eqn{H > 0.5}.
#' @examples
#' \dontrun{
#' ### deterministic signal estimation
#' # set parameters for the example
#' set.seed(3)
#' nlevels = 14
#' n = 2 ^ nlevels
#' xlim = c(200, 400)
#' H = 0.3
#' # simulate a simple fbm + sinusal signal
#' x = 1.25 * cospi(2 * 1:n / 10)
#' y = fbmSim(n = n, H = H) + x
#' plot(y, xlim = xlim,
#'      ylim = range(y[xlim[[1]]:xlim[[2]]]),
#'      type = "l")
#' lines(x, lty = 2, col = 2)
#' legend("topright", lty = 1:2, col = 1:2,
#'        legend = c("Y", "x"), bty = "n")
#' # compute the wavelet transform and the wavelet coefficients' variances
#' wy = wd(y, bc = "symmetric")
#' vpr = waveletVar(wy)
#' # do you note the increase in variance in level 11?
#' plot(vpr, xlim = c(4, nlevels - 1),
#'      ylim = range(vpr[5:nlevels]))
#' # estimate the fBm parameters avoiding levels 11 and 12 (with deterministic
#' # contributions). Level 12 is also avoid as a precaution
#' model = estimateFbmPars(vpr, use_resolution_levels = c(5:10, 13))
#'
#' # Create a fracdet object...
#' fd = fracdet(wy, model)
#' # ... and check the fbm parameters estimates
#' coef(fd)
#' # A plot of the fitted variances may be useful
#' plot(getWaveletVar(fd))
#' points(getFittedWaveletVar(fd),
#'       col = 2,
#'       pch = 2)
#' # we could also extract the nls-fit and use the predict function
#' # The nls-fit is performed in semilog-space. Thus, a transformation
#' # of the predicted values is required
#' nls_model = getWaveletVarModel(fd)
#' points(resolutionLevels(vpr),
#'        2 ^ predict(nls_model, newdata = data.frame(x = resolutionLevels(vpr))),
#'        col = 3,
#'        pch = 3)
#' # Estimate the deterministic signal (this may take a while) taking
#' # into account the deviations in level 11, 12. The estimateDetSignal uses
#' # this information using a Bayesian modelling approach
#' wx = estimateDetSignal(fd, estimate_from = 11:12)
#' x_est = wr(wx)
#' # compare the original signal and the estimation
#' old_par = par(mfrow = c(2,1))
#' plot(y, xlim = xlim,
#'      ylim = range(y[xlim[[1]]:xlim[[2]]]),
#'      type = "l")
#' lines(x, lty = 2, col = 2)
#' legend("topright", lty = 1:2, col = 1:2,
#'        legend = c("Y", "x"), bty = "n")
#' plot(x, type = "l", xlim = xlim,
#'      ylim = range(x) * c(1, 2))
#' lines(x_est, col = 2, lty = 2)
#' legend("topright", lty = 1:2, col = 1:2,
#'        legend = c("x", "x estimate"), bty = "n")
#' par(old_par)
#' }
#' @seealso \code{\link{waveletVar}}, \code{\link{getWaveletVar}},
#' \code{\link{getWaveletVarModel}}, \code{\link{getFittedWaveletVar}},
#' \code{\link{as.wd}}, \code{\link{estimateDetSignal}}
#' @export
#' @exportClass fracdet
#' @export
fracdet = function(x, fbmPars) {
  checkFracdetArgs(x, fbmPars)
  attr(x, "waveletVar") = waveletVar(x)
  attr(x, "fbmPars") = fbmPars
  class(x) = c("fracdet", class(x))
  x
}

# Basic functionality -----------------------------------------------------

#' Get the wavelet coefficients' variances
#'
#' @param x A \code{fracdet} object.
#' @return A \code{waveletVar} object representing the wavelet coefficients'
#' variances of the wavelet transform represented by the \code{fracdet} object.
#' @seealso \code{\link{fracdet}}
#' @export
getWaveletVar = function(x) {
  UseMethod("getWaveletVar", x)
}

#' @export
getWaveletVar.fracdet = function(x) {
  attr(x, "waveletVar")
}


#' Get fitted model of the fractional Brownian motion wavelet variances
#' @inheritParams getWaveletVar
#' @return Either a \code{nls} or a \code{lm} object representing the least-squared
#' estimates of the fractional Brownian motion parameters.
#' @seealso \code{\link{fracdet}}, \code{\link{coef.fracdet}}
#' @export
getWaveletVarModel = function(x){
  UseMethod("getWaveletVarModel", x)
}

#' @export
getWaveletVarModel.fracdet = function(x){
  attr(x, "fbmPars")
}

#' Extract fractional Brownian motion parameters
#' @param object A \code{fracdet} object.
#' @param ... Ignored (used for S3 generic/method consistency).
#' @return Named numeric vector with the least-squared
#' estimates of the fractional Brownian motion parameters
#' (\code{H} and \code{sigma2}).
#' @seealso \code{\link{fracdet}}
#' @export
coef.fracdet = function(object, ...) {
  coef(getWaveletVarModel.fracdet(object))
}

#' Estimates of the wavelet variances of the fractional Brownian motion
#'
#' \code{getFittedVar} estimates the wavelet coefficients' variances of the
#' fractional Brownian motion (fBm) part of the fractal-deterministic model.
#' The variances estimates are obtained considering the theoretical wavelet
#'  variance expression and the estimates of the fBm  parameters.
#'
#' @return A \code{waveletVar} representing the estimates of the wavelet
#' coefficients' variances of the fractional Brownian motion (fBm) part of the
#' fractal-deterministic model.
#' @seealso \code{\link{fracdet}}, \code{\link{theoreticalWaveletVar}}
#' @export
getFittedWaveletVar = function(x){
  UseMethod("getFittedWaveletVar", x)
}

#' @export
getFittedWaveletVar.fracdet = function(x){
  fbm_pars = coef(x)
  theoreticalWaveletVar(fbm_pars[["H"]], fbm_pars[["sigma2"]],
                        x$filter$family,
                        x$filter$filter.number,
                        length(getWaveletVar(x)))
}

#' @rdname getFittedWaveletVar
#' @param x,object A \code{fracdet} object.
#' @param ... Ignored (used for S3 generic/method consistency).
#' @export
fitted.fracdet = function(object, ...) {
  getFittedWaveletVar.fracdet(object)
}

#' Coerce to a wd object
#'
#' Coarce an object to a \code{wd} object if possible.
#' @param x Any R object.
#' @return A \code{wd} object (see \code{\link[wavethresh]{wd}}).
#' @export
as.wd = function(x) {
  UseMethod("as.wd", x)
}

#' @rdname as.wd
#' @export
as.wd.fracdet = function(x){
  class(x) = "wd"
  x
}

# methods of class wd -----------------------------------------------------
# Required due to an error in the "wd" class, that checks the object class by
# direct comparation instead of using inherits

accessC.fracdet = function(x, ...) { accessC.wd(as.wd(x), ...) }
accessD.fracdet = function(x, ...) { accessD.wd(as.wd(x), ...) }
convert.fracdet = function(x, ...) { convert.wd(as.wd(x), ...) }
draw.fracdet = function(x, ...) { draw.wd(as.wd(x), ...) }
image.fracdet = function(x, ...) { image.wd(as.wd(x), ...) }
IsEarly.fracdet = function(x) { IsEarly.wd(as.wd(x)) }
LocalSpec.fracdet = function(x, ...) { LocalSpec.wd(as.wd(x), ...) }
modernise.fracdet = function(x, ...) { modernise.wd(as.wd(x), ...) }
nullevels.fracdet = function(x, ...) { nullevels.wd(as.wd(x), ...) }
plot.fracdet = function(x, ...) { plot.wd(as.wd(x), ...) }
print.fracdet = function(x, ...) { print.wd(as.wd(x), ...) }
putC.fracdet = function(x, ...) { putC.wd(as.wd(x), ...) }
putD.fracdet = function(x, ...) { putD.wd(as.wd(x), ...) }
summary.fracdet = function(x, ...) { summary.wd(as.wd(x), ...) }
threshold.fracdet = function(x, ...) { threshold.wd(as.wd(x), ...) }
waveletVar.fracdet = function(x, ...) { waveletVar(as.wd(x), ...) }
wr.fracdet = function(x, ...) { wr.wd(as.wd(x), ...) }


# estimateDetSignal -----------------------------------------------------------

#' Estimate deterministic signal
#'
#' \code{estimateDetSignal} estimates the deterministic signal based on a
#' fractal-deterministic model that consists on a  linear superposition of a
#' fractional brownian motion (fBm) \eqn{B} and a deterministic band-limited
#' signal \eqn{x} (\eqn{Y = x + B}). The model is represented (in wavelet
#' domain) through a \code{fracdet} object which also contains estimates of the
#' fBm parameters (see \code{\link{fracdet}}). To obtain the estimate of the
#' deterministic component, all the knowledge we may have
#' about the signals is exploited using a bayesian framework. This knowledge
#' consist on:
#' \itemize{
#' \item The estimates of the fBm parameters (which completely characterize the
#' statistical properties of the fBm signal).
#' \item The well-known statistical properties of the fBm signals in wavelet
#' domain.
#' \item The deviations from the theoretical wavelet coefficients' variance
#' (as a function of the resolution level) permit the estimation of the energy
#' distribution accross the resolution levels of the deterministic signal.
#' }
#' Using all this information, a bayesian model of the distribution
#' of the deterministic and stochastic wavelet coefficients is built. The wavelet
#' transform of the deterministic signal is estimated by
#' calculating the decision rule that minimizes the posterior expected value of
#' a squared loss function.
#'
#' @inheritParams getWaveletVar
#' @param estimate_from Numeric vector specifying which resolution levels of the
#' wavelet transform should be used to estimate the deterministic signal.
#' @param df Degrees of freedom of the Student's t used to model the deterministic
#' wavelet coefficients with large energy. Default: {df = 5}.
#' @param nsim_bootstrap Number of bootstrap replicates used to estimate the
#' regression errors in each resolution level. Default: \code{nsim_boostrap = 1000}.
#' @param amplitude_fraction Fraction of the deterministic signal amplitude used
#' to compute a soft bound of what should be considered a negligible
#' deterministic wavelet coefficient. Default: \code{amplitude_fraction = 0.01}.
#' @param signal_amplitude An estimate of the deterministic signal amplitude. If
#' not provided, a proper estimate is computed.
#' @param gk_method An integer code specifying the integration rule to be used
#' in the computations. \code{gk_method} should be an integer between 1 and 6,
#' corresponding to the 15, 21, 31, 41, 51 and 61 point Gauss-Kronrod rules.
#' Default: \code{gk_method = 6}. See details.
#' @param rel_tol Desired relative error used in the numerical integration.
#' Defautl \code{rel_tol = 1e-7}. See details.
#' @param nsubintv  Maximum number of subintervals to be used in the numerical
#' integration. Default: \code{nsubintv = 1e5}. See details.
#' @param aggr_strategy Logical value. If \code{TRUE}, an aggresive strategy
#' for selecting the priors is used, leading to more aggresive estimations
#' of the deterministic signal. Default: \code{aggr_strategy = TRUE}.
#' @return A \code{wd} object representing the wavelet transform of the
#' deterministic signal estimation. To obtain the final estimate of \eqn{x},
#' the \code{wr} function (\code{\link[wavethresh]{wr}}) may be used.
#' @details
#' The decision rule involves the computation of several integrals. The numerical
#' integrals were computed using an adaptive integration scheme from the
#' \href{https://www.gnu.org/software/gsl/}{GSL - GNU Scientific Library}. For
#' further details about the numerical integration method, the interested reader
#'  is referred to the
#' \href{https://www.gnu.org/software/gsl/manual/html_node/QAG-adaptive-integration.html}{GSL documentation}.
#' @examples
#' \dontrun{
#' # set parameters for the example
#' set.seed(3)
#' nlevels = 14
#' n = 2 ^ nlevels
#' xlim = c(200, 400)
#' H = 0.3
#' # simulate a simple fbm + sinusal signal
#' x = 1.25 * cospi(2 * 1:n / 10)
#' y = fbmSim(n = n, H = H) + x
#' # compute the wavelet transform and the wavelet coefficients' variances
#' wy = wd(y, bc = "symmetric")
#' vpr = waveletVar(wy)
#' # do you note the increase in variance in level 11?
#' plot(vpr, xlim = c(4, nlevels - 1),
#'      ylim = range(vpr[5:nlevels]))
#' # estimate the fBm parameters avoiding levels 11 and 12 (with deterministic
#' # contributions). Level 12 is also avoid as a precaution
#' model = estimateFbmPars(vpr, use_resolution_levels = c(5:10, 13))
#' # Create a fracdet object...
#' fd = fracdet(wy, model)
#' # ... and check the fit
#' print(coef(fd))
#' # plot the experimental and the fitted wavelet's variances
#' plot(getWaveletVar(fd))
#' points(getFittedWaveletVar(fd),
#'        col = 2,
#'        pch = 2)
#' # Estimate the deterministic signal (this may take a while) taking
#' # into account the deviations in level 11, 12. The estimateDetSignal uses
#' # this information using a Bayesian modelling approach
#' wx = estimateDetSignal(fd, estimate_from = 11:12)
#' x_est = wr(wx)
#' # compare the original signal and the estimation
#' old_par = par(mfrow = c(2,1))
#' plot(y, xlim = xlim,
#'      ylim = range(y[xlim[[1]]:xlim[[2]]]),
#'      type = "l")
#' lines(x, lty = 2, col = 2)
#' legend("topright", lty = 1:2, col = 1:2,
#'        legend = c("Y", "x"), bty = "n")
#' plot(x, type = "l", xlim = xlim,
#'      ylim = range(x) * c(1, 2))
#' lines(x_est, col = 2, lty = 2)
#' legend("topright", lty = 1:2, col = 1:2,
#'        legend = c("x", "x estimate"), bty = "n")
#' par(old_par)
#' }
#' @seealso \code{\link[wavethresh]{wd}}, \code{\link[wavethresh]{wr}},
#' \code{\link{fracdet}}.
#' @export
estimateDetSignal = function(x, estimate_from,
                             df = 5, nsim_bootstrap = 1000,
                             amplitude_fraction = 0.01,
                             signal_amplitude = NULL,
                             gk_method = 6, rel_tol = 1e-7,
                             nsubintv = 1e5L,
                             aggr_strategy = TRUE){
  UseMethod("estimateDetSignal", x)
}

# gk_method = 1:6, corresponding to the 15, 21, 31, 41, 51 and 61 point Gauss-Kronrod rules
checkGKMethod = function(code){
  if (length(code) > 1) {
    stop("gk_method should have length 1")
  }
  if (!(code %in% 1:6)) {
    stop("Invalid gk_method: gk_method should be an integer from 1 to 6")
  }
}

# gk_method = 1:6, corresponding to the 15, 21, 31, 41, 51 and 61 point Gauss-Kronrod rules
#' @export
estimateDetSignal.fracdet = function(x, estimate_from,
                                     df = 5, nsim_bootstrap = 1000,
                                     amplitude_fraction = 0.01,
                                     signal_amplitude = NULL,
                                     gk_method = 6, rel_tol = 1e-7,
                                     nsubintv = 1e5L,
                                     aggr_strategy = TRUE) {
  checkGKMethod(gk_method)
  # extract variables for the sake of clarity
  vpr = getWaveletVar(x)
  nlevels = length(vpr)
  fitted_vpr = getFittedWaveletVar(x)
  fitted_model = getWaveletVarModel(x)
  family = x$filter$family
  filter_number = x$filter$filter.number
  fitted_vpr_std = bootstrapWaveletVarError(fitted_model, nsim_bootstrap,
                                            family, filter_number, nlevels)
  if (is.null(signal_amplitude)) {
    signal_amplitude = mad(diff(wr(x)))
  }
  wavelet_amplitudes = getWaveletAmplitudes(family, filter_number,nlevels)
  x = as.wd(x)
  # get the lengths of the wavelet coefficients in each wavelet resolution level
  details_lengths =  sapply(0:(nlevels - 1),
                            FUN = function(rl) {
                              length(accessD(x, level = rl, boundary = TRUE))
                            })

  # assume that the slowest trend is due to the fbm
  x = putC(x,
           level = 0,
           rep(0, details_lengths[[1]]),
           boundary = TRUE)

  # eliminate all the signal in those levels where the deterministic
  # contribution is negligible
  for (rl in setdiff(resolutionLevels(vpr), estimate_from)) {
    x = putD(x,
             level = rl,
             rep(0, details_lengths[[rl + 1]]),
             boundary = TRUE)
  }

  # process those resolution levels with deterministic contributions
  for (rl in estimate_from) {
    # check that the fitted variance <= experimental variance
    # in those resolution levels where we assume there is deterministic
    # contribution
    if (vpr[rl + 1] <= fitted_vpr[rl + 1]) {
      x = putD(x, level = rl,
               rep(0, details_lengths[[rl + 1]]),
               boundary = TRUE)
    } else {
      priors = getPriors(accessD(x, rl, boundary = TRUE),
                         vpr[[rl + 1]], fitted_vpr[[rl + 1]],
                         fitted_vpr_std[[rl + 1]],
                         family, filter_number, nlevels,
                         df, signal_amplitude, amplitude_fraction,
                         wavelet_amplitudes[[rl + 1]],
                         aggr_strategy)
      x = putD(x, level = rl, boundary = TRUE,
               bayesRule(accessD(x, rl, boundary = TRUE),
                         priors, gk_method, rel_tol, nsubintv)
               )
    }
  }
  x
}

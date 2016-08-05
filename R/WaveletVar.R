
# constructors of the WaveletVar class  -----------------------------------------------

# Build a WaveletVar object from a numeric vector (private function)
#
# @param vpr Numeric vector containing the  wavelet coefficients' variances on
# each wavelet resolution level (resolution levels are assumed to be arranged
# in ascending order).
# @param family Specifies the family of wavelets used in the computation of the
# wavelet transform. See \code{\link[wavethresh]{wd}}.
# @param filter_number An integer identifying the concrete filter used within the
# \code{family} of wavelets selected. \code{filter_number} is related with the
# smoothness of the wavelet. See \code{\link[wavethresh]{wd}}.
buildWaveletVar = function(vpr, family, filter_number) {
  names(vpr) = 0:(length(vpr) - 1)
  class(vpr) = c("WaveletVar", class(vpr))
  attr(vpr, "wt_info") = list(family = family,
                              filter_number = filter_number)
  vpr
}

#' Wavelet coefficients' variances.
#'
#' \code{WaveletVar} computes the wavelet coefficients' variances of an object
#' representing a wavelet transform in each of its resolution levels.
#'
#' In addition to the specific methods for the \code{WaveletVar} class, all the
#' methods for numeric vectors can be used with a \code{WaveletVar} object.
#'
#' @param x A \code{wd} object (see \code{\link[wavethresh]{wd}}). The wavelet
#' coefficients' variances are computed in each of the resolution levels
#' of \code{wd}.
#' @return A S3 \code{WaveletVar} object that stores the
#' Wavelet coefficients' variances depending on the resolution level.
#' @export
#' @exportClass Fracdet
#' @examples
#' use_H = 0.3
#' fbm = fbmSim(n = 2 ^ 10, H = use_H)
#' w_fbm = wd(fbm, bc = "symmetric")
#' vpr = WaveletVar(w_fbm)
#' plot(vpr)
#' # the theoreticalWaveletVar also returns a 'WaveletVar' object:
#' theo_vpr = theoreticalWaveletVar(H = use_H, sigma2 = 1,
#'                                  family = wtInfo(vpr)[["family"]],
#'                                  filter_number = wtInfo(vpr)[["filter_number"]],
#'                                  nlevels = length(vpr))
#' points(theo_vpr, col = "red")
#' @seealso \code{\link{wtInfo}}, \code{\link{resolutionLevels}},
#' \code{\link{theoreticalWaveletVar}}, \code{\link{Fracdet}},
#' \code{\link{estimatefBmPars}}
#' @import wavethresh
#' @export
WaveletVar = function(x){
  kMIN_POINTS = 5
  if (!inherits(x, "wd")) {
    stop("x is not an \"wd\" object")
  }
  x_nlevels = wavethresh::nlevelsWT(x)
  vpr = rep(0.0, x_nlevels)

  filter_len = length(
    wavethresh::filter.select(
      filter.number = x$filter$filter.number,
      family = x$filter$family
    )$H
  )
  # first index unaffected by the circularity assumption of the wavelet transform
  unaffected_index =  ceiling( (filter_len - 2) * (1 - 1 / 2 ^ (x_nlevels:1)) )

  # estimate the variance of the resolution level 0 ...
  vpr[[1]] = wavethresh::accessD(x, 0) ^ 2
  # ... and loop for the remainder levels
  resolution_levels = 0:(x_nlevels - 1)
  for (i in 2:x_nlevels) {
    available_points = (2 ^ resolution_levels[[i]] - 2 * unaffected_index[[i]])
    if (available_points >  kMIN_POINTS) {
      segment_index = unaffected_index[[i]]:(2 ^ resolution_levels[[i]] -
                                               unaffected_index[[i]])
      vpr[[i]] = var(wavethresh::accessD(x,
                                         resolution_levels[[i]])[segment_index])

    }else {
      # not enough points for a variance-estimation free from the circularity
      # assumption: we use all points of the level
      vpr[[i]] = var(wavethresh::accessD(x, resolution_levels[[i]]))
    }
  }

  buildWaveletVar(vpr,
                  family = x$filter$family,
                  filter_number = x$filter$filter.number)
}

# Accesing attributes of the WaveletVar class ------------------------------------

#' Get the resolution levels at which the wavelet coefficients' variances was
#' computed.
#' @param x A \code{WaveletVar} object.
#' @return A numerical vector of the same length as \code{WaveletVar}.
#' @export
resolutionLevels = function(x){
  UseMethod("resolutionLevels")
}

#' @export
resolutionLevels.WaveletVar = function(x){
  0:(length(x) - 1)
}

#' Get the family and the filter-number used in the wavelet transform asociated
#' with the \code{WaveletVar} object.
#' @inheritParams resolutionLevels
#' @return An R list containing the family (\code{family} field) and the
#' filter-number (\code{filter_number}).
#' @export
wtInfo = function(x){
  UseMethod("wtInfo")
}

#' @export
wtInfo.WaveletVar = function(x){
  attr(x, "wt_info")
}

# Plotting functions for WaveletVar class ----------------------------------------
#' @export
plot.WaveletVar = function(x, main = "Variance per Resolution Level",
                           xlab = "Wavelet Resolution Level",
                           ylab = "Variance", log = "y", ...){
  plot(resolutionLevels(x), as.numeric(x), main = main, xlab = xlab, ylab = ylab,
       log = log, ...)
}

#' @export
points.WaveletVar = function(x, ...) {
  points(resolutionLevels(x), as.numeric(x), ...)
}

#' @export
lines.WaveletVar = function(x, ...){
  lines(resolutionLevels(x), as.numeric(x), ...)
}


# Theoretical expression of WaveletVar -----------------------------------------------
# private function
# get wavelet filters

getWaveletFilters = function(family, filter_number){
  h = wavethresh::filter.select(filter.number = filter_number,
                                family = family)$H
  len = length(h)
  # Definition of g[n]:
  # g[n] = (-1) ^ (1 - n) * h[1 - n]
  index = 1 - (len - 1):0
  g = rev(h) * (-1) ^ (1 - index)
  list(h = h, g = g)
}


#' Wavelet coefficients' variances of a fractional Brownian motion
#'
#' \code{theoreticalWaveletVar} calculates the theoretical wavelet coefficients'
#'  variance of a fractional Brownian motion (fBm) with known parameters.
#'
#' @param H The Hurst exponent of the fBm.
#' @param sigma2  Variance of the fractional Gaussian noise that results after
#' differentiating the fBm.
#' @param family Specifies the family of wavelets used in the computation of the
#' wavelet transform. See \code{\link[wavethresh]{wd}}.
#' @param filter_number An integer identifying the concrete filter used within the
#' \code{family} of wavelets selected. \code{filter_number} is related with the
#' smoothness of the wavelet. See \code{\link[wavethresh]{wd}}.
#' @param nlevels Number of resolution levels to compute. Thus, the resulting
#' wavelet coefficients' variances are the expected ones from a fBm of length
#' \eqn{2 ^ n}.
#' @return A \code{WaveletVar} object with the theoretical wavelet variance.
#' @examples
#' use_H = 0.3
#' fbm = fbmSim(n = 2 ^ 10, H = use_H)
#' w_fbm = wd(fbm, bc = "symmetric")
#' vpr = WaveletVar(w_fbm)
#' plot(vpr)
#' # the theoreticalWaveletVar also returns a 'WaveletVar' object:
#' theo_vpr = theoreticalWaveletVar(H = use_H, sigma2 = 1,
#'                                  family = wtInfo(vpr)[["family"]],
#'                                  filter_number = wtInfo(vpr)[["filter_number"]],
#'                                  nlevels = length(vpr))
#' points(theo_vpr, col = "red")
#' @seealso \code{\link{WaveletVar}}, \code{\link{estimatefBmPars}}
#' @export
#' @useDynLib fracdet
theoreticalWaveletVar = function(H, sigma2, family, filter_number, nlevels) {
  filters = getWaveletFilters(family = family,
                              filter_number = filter_number)
  buildWaveletVar(.Call('fracdet_theoreticalWaveletVarCpp', PACKAGE = 'fracdet',
                        filters$g, filters$h, H, sigma2, nlevels),
                  family = family,
                  filter_number = filter_number)
}





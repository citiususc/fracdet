
# constructors of the WaveletVar class  -----------------------------------------------

#' Build a WaveletVar object from a numeric vector (private function)
#' @param family Specifies the family of wavelets used in the computation of the
#' wavelet transform. See \code{\link[wavethresh]{wd}}.
#' @param filter_number An integer identifying the concrete filter used within the
#' \code{family} of wavelets selected. \code{filter_number} is related with the
#' smoothness of the wavelet. See \code{\link[wavethresh]{wd}}.
buildWaveletVar = function(vpr, family, filter_number) {
  names(vpr) = 0:(length(vpr) - 1)
  class(vpr) = c("WaveletVar", class(vpr))
  attr(vpr, "wt_info") = list(family = family,
                              filter_number = filter_number)
  vpr
}

#' Wavelet coefficients' variance depending on the resolution level
#'
#' @param x Either a numeric vector or a \code{wd} object (see
#' \code{\link[wavethresh]{wd}}). If \code{x} is a numeric vector, it should
#' contain the wavelet coefficients' variance on each wavelet resolution level
#' (resolution levels are assumed to be arranged in ascending order). If
#' \code{x} is a \code{wd} object, the wavelet coefficients' variance is computed
#' before constructing the \code{WaveletVar} object.
#' @param ... Additional arguments.
#' @return An S3 \code{WaveletVar} object that stores the
#' Wavelet coefficients' variance depending on the resolution level.
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
#' @import wavethresh
#' @rdname WaveletVar
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

#' @section Methods:
#' \code{resolutionLevels}: Get the resolution levels.
#' @rdname WaveletVar
#' @export
resolutionLevels = function(x){
  UseMethod("resolutionLevels")
}

#' @export
resolutionLevels.WaveletVar = function(x){
  0:(length(x) - 1)
}

#' @section Methods:
#' \code{wtInfo}: Get the family and the filter-number used in the
#'  wavelet transform asociated with the object. Returns a list with the proper
#'  wavelet information.
#' @rdname WaveletVar
#' @export
wtInfo = function(x){
  UseMethod("wtInfo")
}

#' @export
wtInfo.WaveletVar = function(x){
  attr(x, "wt_info")
}

# Plotting functions for WaveletVar class ----------------------------------------
#' @section Methods:
#' \code{plot}, \code{points} and \code{lines}: generic plotting utilities.
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





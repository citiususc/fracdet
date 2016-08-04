
# constructors of the WaveletVar class  -----------------------------------------------



buildWaveletVar = function(vpr, family, filter_number){
  names(vpr) = 0:(length(vpr) - 1)
  class(vpr) = c("WaveletVar", class(vpr))
  attr(vpr, "wt_info") = list(family = family,
                              filter_number = filter_number)
  vpr
}

#' constructor of the WaveletVar class from a wavethresh::wd object
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

#' @export
resolutionLevels = function(x){
  UseMethod("resolutionLevels")
}

#' @export
resolutionLevels.WaveletVar = function(x){
  0:(length(x) - 1)
}

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





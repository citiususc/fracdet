# fracdet constructor -----------------------------------------------------

#' @export
Fracdet = function(x, fbmPars) {
  UseMethod("Fracdet", fbmPars)
}

# private function
FracdetFromFit = function(x, fbmPars) {
  if (!inherits(x, "wd")) {
    stop("x is not a \"wd\" object")
  }
  attr(x, "WaveletVar") = WaveletVar.wd(x)
  attr(x, "fbmPars") = fbmPars
  class(x) = c("Fracdet", class(x))
  x
}

#' @export
Fracdet.lm = function(x, fbmPars) {
  FracdetFromFit(x, fbmPars)
}

#' @export
Fracdet.nls = function(x, fbmPars) {
  FracdetFromFit(x, fbmPars)
}

# Basic functionality -----------------------------------------------------

#' @export
getWaveletVar = function(x) {
  UseMethod("getWaveletVar", x)
}

#' @export
getWaveletVar.Fracdet = function(x) {
  attr(x, "WaveletVar")
}

#' @export
getFbmPars = function(x){
  UseMethod("getFbmPars", x)
}

#' @export
getFbmPars.Fracdet = function(x){
  attr(x, "fbmPars")
}

#' @export
coef.Fracdet = function(x) {
  coef(getFbmPars.Fracdet(x))
}

#' @export
getFittedWaveletVar = function(x){
  UseMethod("getFittedWaveletVar", x)
}

#' @export
getFittedWaveletVar.Fracdet = function(x){
  fbm_pars = coef(x)
  theoreticalWaveletVar(fbm_pars[["H"]], fbm_pars[["sigma2"]],
                        x$filter$family,
                        x$filter$filter.number,
                        length(getWaveletVar(x)))
}

#' @export
fitted.Fracdet = function(x) {
  getFittedWaveletVar.Fracdet(x)
}


#' @export
as.wd = function(x) {
  UseMethod("as.wd", x)
}

#' @export
as.wd.Fracdet = function(x){
  class(x) = "wd"
  x
}

# methods of class wd -----------------------------------------------------
# Required due to an error in the "wd" class, that checks the object class by
# direct comparation instead of using inherits

accessC.Fracdet = function(x, ...) { accessC.wd(as.wd(x), ...) }
accessD.Fracdet = function(x, ...) { accessD.wd(as.wd(x), ...) }
convert.Fracdet = function(x, ...) { convert.wd(as.wd(x), ...) }
draw.Fracdet = function(x, ...) { draw.wd(as.wd(x), ...) }
image.Fracdet = function(x, ...) { image.wd(as.wd(x), ...) }
IsEarly.Fracdet = function(x) { IsEarly.wd(as.wd(x)) }
LocalSpec.Fracdet = function(x, ...) { LocalSpec.wd(as.wd(x), ...) }
modernise.Fracdet = function(x, ...) { modernise.wd(as.wd(x), ...) }
nullevels.Fracdet = function(x, ...) { nullevels.wd(as.wd(x), ...) }
plot.Fracdet = function(x, ...) { plot.wd(as.wd(x), ...) }
print.Fracdet = function(x, ...) { print.wd(as.wd(x), ...) }
putC.Fracdet = function(x, ...) { putC.wd(as.wd(x), ...) }
putD.Fracdet = function(x, ...) { putD.wd(as.wd(x), ...) }
summary.Fracdet = function(x, ...) { summary.wd(as.wd(x), ...) }
threshold.Fracdet = function(x, ...) { threshold.wd(as.wd(x), ...) }
WaveletVar.Fracdet = function(x, ...) { WaveletVar.wd(as.wd(x), ...) }
wr.Fracdet = function(x, ...) { wr.wd(as.wd(x), ...) }


# estimateDetSignal -----------------------------------------------------------

#' @export
estimateDetSignal = function(x, ...){
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
estimateDetSignal.Fracdet = function(x, estimate_from,
                                     n1 = 5, nsim_bootstrap = 1000,
                                     amplitude_percentage = 0.01,
                                     signal_amplitude = NULL,
                                     gk_method = 6, rel_tol = 1e-7,
                                     nsubintv = 1e5L) {
  checkGKMethod(gk_method)
  # extract variables for the sake of clarity
  vpr = getWaveletVar(x)
  nlevels = length(vpr)
  fitted_vpr = getFittedWaveletVar(x)
  fitted_model = getFbmPars(x)
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
                         n1, signal_amplitude, amplitude_percentage,
                         wavelet_amplitudes[[rl + 1]])
      x = putD(x, level = rl, boundary = TRUE,
               bayesRule(accessD(x, rl, boundary = TRUE),
                         priors, gk_method, rel_tol, nsubintv)
               )
    }
  }
  x
}

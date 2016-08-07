
# estimate start, upper and lower for nls ---------------------------------


getInitialParameters = function(nls_data, family, filter_number){
  if (nrow(nls_data) == 1) {
    start = list(H = 0.5, sigma2 = 2 ^ nls_data$y)
  } else {
    # y ~ x approximately follows a linear model
    lm_fit = lm(y ~ x, data = nls_data)
    initial_H = (-coef(lm_fit)[[2]] - 1) / 2
    # check that 0 < H < 1
    initial_H = ifelse(initial_H <= 0, 0.05, initial_H)
    initial_H = ifelse(initial_H >= 1, 0.95, initial_H)
    initial_sigma2 =
      2 ^ predict(lm_fit, newdata = data.frame(x = max(nls_data$x))) /
      theoreticalWaveletVar(initial_H, 1, family, filter_number, 1)
    start = list(H = initial_H, sigma2 = initial_sigma2)
  }
  start
}

getNlsLower = function(vpr, family, filter_number, initial_H){
  kDELTA = 1e-12
  kSHRINK = 1e-3
  minimum_sigma2 =
    min(vpr, na.rm = TRUE) / theoreticalWaveletVar(initial_H, 1, family, filter_number, 1)
  # divide the minimum estimated sigma2 to have some safety margin
  list(H = 0 + kDELTA, sigma2 = minimum_sigma2 * kSHRINK)
}

getNlsUpper = function(vpr, family, filter_number, initial_H){
  kDELTA = 1e-12
  kENLARGE = 1e3
  maximum_sigma2 =
    max(vpr, na.rm = TRUE) / theoreticalWaveletVar(initial_H, 1, family, filter_number, 1)
  list(H = 1 - kDELTA, sigma2 = maximum_sigma2 * kENLARGE)
}

# nls adapter -------------------------------------------------------------


nlsAdapter = function(...){
  args = list(...)
  # the data for the regression is passed as nls_data to avoid collisions just in
  # case the user provides data as part of ... in the estimateFbmPars function.
  # Rename it and eliminate nls_data
  args$data = args$nls_data
  args$nls_data = NULL
  # if weights is NULL, remove it from the list
  if (is.null(args$weights)) {
    args$weights = NULL
  }
  # call nls
  do.call(nls, args)
}


# estimateFbmPars ---------------------------------------------------------


#' Estimate fractional Brownian motion parameters from its wavelet coefficients' variances
#'
#' \code{estimateFbmPars} determines the fractional Brownian motion (fBm)
#' parameters by fitting the theoretical wavelet coefficients' variances
#' to given experimental variances. The estimation procedure is thus based
#' on nonlinear least-squares estimates and makes use of the \code{nls} function
#' \code{\link[stats]{nls}}.
#'
#' @param x A \code{waveletVar} object.
#' @param use_resolution_levels A numeric vector specifying which resolution levels
#' should be used in the fit.
#' @param start A named list or named numeric vector of starting estimates. The
#' names of the list should be \code{H} and \code{sigma2}.
#' @param lower,upper Vectors of lower and upper bounds. Note that the \code{H}
#' parameter should fulfil \eqn{0 < H < 1}, whereas that \code{sigma2} should be
#' \eqn{sigma2 > 0} If unspecified, proper bounds are computed.
#' @param use_weights Logical value indicating wheter the objective function is
#' weighted least squares or not. If \code{use_weights} is \code{TRUE}, proper
#' weights are automatically computed.
#' @param algorithm Character string specifying the algorithm to use (see
#' \code{\link[stats]{nls}}. Note that bounds can only be used with the "port"
#' algorithm (default choice).
#' @param ... Additional \code{nls} arguments (see \code{\link[stats]{nls}}).
#' @return A \code{nls} object representing the fitted model. See
#' \code{\link[stats]{nls}} for further details.
#' @examples
#' set.seed(10)
#' fbm = fbmSim(n = 2 ^ 13, H = 0.4)
#' vpr = waveletVar(wd(fbm, bc = "symmetric"))
#' plot(vpr)
#' # Estimate the fBm parameters using the largest resolution levels
#' # since the estimates of their variances are better
#' model = estimateFbmPars(vpr, use_resolution_levels = 5:12)
#' # The nls-fit is performed in semilog-space. Thus, a transformation
#' # of the predicted values is required
#' points(resolutionLevels(vpr),
#'        2 ^ predict(model, newdata = data.frame(x = 0:12)),
#'        col = 2,
#'        pch = 2)
#' # Since the estimates of the largest resolution levels are better we may
#' # use a weigthed regression scheme
#' wmodel = estimateFbmPars(vpr, use_resolution_levels = 5:12,
#'                          use_weights = TRUE)
#' points(resolutionLevels(vpr),
#'        2 ^ predict(wmodel, newdata = data.frame(x = 0:12)),
#'        col = 3,
#'        pch = 3)
#' legend("topright", pch = 1:3, col = 1:3, bty = "n",
#'        legend = c("Original data", "nls model", "weighted-nls model"))
#' @section Warning:
#' Do not use nls on artificial "zero-residual" data. See
#' \code{\link[stats]{nls}} for further details.
#' @seealso \code{\link{waveletVar}}, \code{\link{theoreticalWaveletVar}}
#' @export
estimateFbmPars = function(x, use_resolution_levels,
                           start = NULL, lower = NULL, upper = NULL,
                           use_weights = FALSE, algorithm = "port",
                           ...) {
  UseMethod("estimateFbmPars", x)
}



#' @export
estimateFbmPars.waveletVar <- function(x,
                                       use_resolution_levels,
                                       start = NULL, lower = NULL, upper = NULL,
                                       use_weights = FALSE, algorithm = "port",
                                       ...){
  # rename for clarity
  vpr = x
  # change resolution levels to vector indexes
  use_indx = use_resolution_levels + 1
  # get wavelet information
  filter_number = wtInfo(vpr)$filter_number
  family = wtInfo(vpr)$family
  nlevels = length(vpr)


  # prepare the data for the nls function:
  nls_data = data.frame(y = log2(vpr)[use_indx],
                        x = resolutionLevels(vpr)[use_indx])

  # select nls parameters if not provided
  if (is.null(start)) {
    start = getInitialParameters(nls_data, family, filter_number)
  }
  if (is.null(lower)) {
    lower = getNlsLower(vpr, family, filter_number, start[["H"]])
  }
  if (is.null(upper)) {
    upper = getNlsUpper(vpr, family, filter_number, start[["H"]])
  }

  # define function to fit
  functionToFit = function(x, H, sigma2) {
    theoreticalWaveletVar(H, sigma2, family, filter_number, nlevels)[x + 1]
  }

  if (use_weights) {
    weights = 2 ^ (nls_data$x) - 1
  } else {
    weights = NULL
  }

  fit = nlsAdapter(formula = y ~ log2(functionToFit(x, H, sigma2)),
                    nls_data = nls_data,
                    weights = weights,
                    algorithm = algorithm,
                    lower = lower,
                    upper = upper,
                    start = start,
                    ...)


  fit

}


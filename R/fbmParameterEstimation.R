
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


nls_adapter = function(...){
  args = list(...)
  # the data for the regression is passed as nls_data to avoid collisions just in
  # case the user provides data as part of ... in the estimatefBmPars function.
  # Rename it and eliminate nls_data
  args$data = args$nls_data
  args$nls_data = NULL
  # if weigths is NULL, remove it from the list
  if (is.null(args$weigths)) {
    args$weights = NULL
  }
  # call nls
  do.call(nls, args)
}


# estimatefBmPars ---------------------------------------------------------


#' @export
estimatefBmPars = function(x, ...) {
  UseMethod("estimatefBmPars", x)
}



#' @export
estimatefBmPars.WaveletVar <- function(vpr,
                                       use_resolution_levels,
                                       start = NULL, lower = NULL, upper = NULL,
                                       use_weights = FALSE, algorithm = "port",
                                       ...){
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
    weigths = 2 ^ (nls_data$x) - 1
  } else {
    weigths = NULL
  }

  fit = nls_adapter(formula = y ~ log2(functionToFit(x, H, sigma2)),
                    nls_data = nls_data,
                    weights = weigths,
                    algorithm = algorithm,
                    lower = lower,
                    upper = upper,
                    start = start,
                    ...)


  fit

}


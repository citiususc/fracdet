# Calculate WaveletVar error ----------------------------------------------


getInvalidSimIndex = function(pars_sim){
  which((pars_sim[,"H"] <= 0) | (pars_sim[,"H"] >= 1) |
          pars_sim[,"sigma2"] <= 0)
}

#' @importFrom MASS mvrnorm
bootstrapWaveletVarError = function(fitted_model, nsim,
                                    family, filter_number, nlevels) {

  fbm_coefs = coef(fitted_model)

  ## get variance-covariance matrix
  cov_matrix = vcov(fitted_model)

  # simulate nsim pars_sim
  pars_sim = MASS::mvrnorm(n = nsim,
                           mu = fbm_coefs,
                           Sigma = cov_matrix,
                           empirical = TRUE)
  # check for possible simulation errors
  invalid_index = getInvalidSimIndex(pars_sim)
  if (length(invalid_index) > 0) {
    pars_sim = pars_sim[-invalid_index,]
    # try to reach nsim valid simulations
    # nsim_additional is computed taking into account that the prob of a succesful
    # simulation is (nsim - length(invalid_index)) / nsim
    nsim_additional = nsim * length(invalid_index) / (nsim - length(invalid_index))
    pars_sim_additional = MASS::mvrnorm(n = max(round(nsim_additional * 1.01), 2),
                                        mu = fbm_coefs,
                                        Sigma = cov_matrix,
                                        empirical = TRUE)
    invalid_index = getInvalidSimIndex(pars_sim_additional)
    if (length(invalid_index) > 0) {
      pars_sim_additional = pars_sim_additional[-invalid_index,]
    }
    pars_sim = rbind(pars_sim, pars_sim_additional)
    if (nrow(pars_sim) == 0) {
      stop("missing boostrapped fbm parameters")
    }
  }
  sim_wavelet_var =
    apply(pars_sim, MARGIN = 1,
          function(fbm_pars) {
            theoreticalWaveletVar(fbm_pars[["H"]], fbm_pars[["sigma2"]],
                                  family, filter_number, nlevels)

          })

  apply(sim_wavelet_var, MARGIN = 1,
        FUN = sd)
}


# Calculate proper priors -------------------------------------------------


getWaveletAmplitudes =  function(family,
                                 filter_number,
                                 nlevels) {
  zwd = wd(rep(0, 2 ^ (nlevels)),
           family = family,
           filter.number = filter_number)
  amplitudes = rep(0, nlevels)

  # get amplitude at level 0
  amplitudes[[1]] = diff(range(wr(putD(zwd, level = 0, 1))))
  # get the amplitudes at the remaining resolution levels
  for (j in 1:(nlevels - 1)){
    put_values = rep(0, 2 ^ j)
    put_values[[ 2 ^ j / 2 ]] = 1
    amplitudes[[ j + 1 ]] = diff(range(wr(putD(zwd, level = j, put_values))))
  }
  amplitudes
}

getPriors = function(wavelet_details, vpr,
                     fitted_vpr, fitted_vpr_std,
                     family, filter_number, nlevels,
                     df_tstudent, signal_amplitude,
                     amplitude_percentage, wavelet_amplitude,
                     aggr_strategy){

  # Compute parameters for the Gamma distribution
  alpha = (fitted_vpr / fitted_vpr_std) ^ 2
  beta = alpha * fitted_vpr

  # Compute tau (dispersion of small deterministic coefficients)
  tau = amplitude_percentage * signal_amplitude / (3 * wavelet_amplitude)

  # Compute p (weight of the large deterministic coeff)
  ncoeffs = length(wavelet_details)
  noise_level = sqrt(2 * log(ncoeffs) * (fitted_vpr + tau ^ 2))

  p = sum(abs(wavelet_details) >= noise_level) / ncoeffs
  if (aggr_strategy || (p <= 1 / ncoeffs)) {
    # use two standard deviations to get a conservative estimation of noise level
    # and thus avoiding estimation p = 0
    conservative_noise_level = 2 * sqrt(fitted_vpr + tau ^ 2)
    p = max(1 / ncoeffs,
            sum(abs(wavelet_details) >= conservative_noise_level) / ncoeffs)

  }

  # Compute mu (dispersion of large deterministic coeff)
  mu2 = ((vpr - fitted_vpr - (1 - p) * tau ^ 2) / p) *
    (df_tstudent - 2) / df_tstudent
  mu = sqrt(max(0, mu2))

  list(alpha = alpha, beta = beta, tau = tau, mu = mu, p = p, n = df_tstudent)
}


# Bayes rule --------------------------------------------------------------
checkPriors = function(priors) {
  required_priors = c("alpha", "beta", "tau", "mu", "p", "n")
  if (!all(required_priors %in% names(priors))) {
    stop("Incorrect priors")
  }
}

bayesRule = function(input, priors, key = 6L, rel_tol = 1e-7, size = 1e5L){
  checkPriors(priors)
  .Call('fracdet_bayesRuleCpp', PACKAGE = 'fracdet',
        input, priors, key, 0.0, rel_tol, size)
}

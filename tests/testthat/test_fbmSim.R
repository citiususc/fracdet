skip_on_cran()

library('fracdet', quietly = TRUE)

test_that("fbmSim simulates the proper covariance structure", {
  # tolerance has to be quite high due to the difficulties in simulating 
  # the highest Hs
  kTOL = 2e-2
  n = 1000
  repl = 500
  max_lag = 5
  k = 0:(max_lag)
  H_to_test = c(0.1,0.3,0.5,0.6,0.7)
  for (H in H_to_test) {
    set.seed(10)
    fgn_acv = rowMeans(replicate(n = repl, 
                                 as.numeric(acf(diff(fbmSim(n,H)), 
                                                plot = FALSE,
                                                type = "covariance",
                                                lag.max = max_lag)$acf)))
    theo_acv = (abs(k - 1) ^ (2 * H) - 2 * abs(k) ^ (2 * H) + 
                  abs(k + 1) ^ (2 * H)) / 2
    
    expect_true(max(abs(fgn_acv - theo_acv)) < kTOL)
  }
})


  

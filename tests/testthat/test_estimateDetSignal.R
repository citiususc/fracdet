skip_on_cran()

library('fracdet', quietly = TRUE)

test_that("estimateDetSignal provides reasonable estimations", {

  kMAX_REL_MSE = 0.1
  nlevels = 13

  nbootstrap = 500
  n = 2 ^ nlevels
  snr = 2

  for (H in c(0.3, 0.5, 0.7)) {
    b = fbmSim(n = n, H = H)
    for (j in c(1,3)) {
      estimate_from = nlevels - j
      use_freq = (0.25 + 0.5) / 2 ^ j
      x = cospi(2 * use_freq * 1:n)
      snr = ifelse(j == 1, 2, 4)
      x =  snr * x / sd(x)
      n1 = sample(c(4, 5, 10), 1)
      amplitude_percentage = sample(c(0.01,0.05,0.1), 1)
      s = b + x
      plot(s, type= "l")
      ws = wd(s, bc = "symmetric")
      vpr = WaveletVar(ws)
      use_resolution_levels = setdiff(4:(nlevels - 1), c(estimate_from,estimate_from + 1))
      fit = fracdet::estimatefBmPars(vpr,
                                     use_resolution_levels)
      # print(coef(fit))
      fdmodel = Fracdet(ws, fit)
      # plot(vpr,main = paste(H, j, snr))
      # points(getFittedWaveletVar(fdmodel),
      #        col = 2, pch = 24, bg = 2)

      est = wr(estimateDetSignal(fdmodel, estimate_from = estimate_from))
      # plot(est, type = "l", xlim = c(0,100))
      # lines(x, col = 2)
      squared_error = (est - x)  ^ 2
      expect_true(mean(squared_error) / var(x) < kMAX_REL_MSE)
    }
  }
})




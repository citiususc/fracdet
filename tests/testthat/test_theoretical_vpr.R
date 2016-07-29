skip_on_cran()

library('fracdet')
data = readRDS("../testdata/vpr_data.RDS")

test_that("theoretical vpr works properly (intensive)", {
  for (i in seq_len(nrow(data))) {
    vpr = theoreticalWaveletVar(H = data[i,"H"],sigma2 = data[i,"sigma2"],
                                family = data[i,"family"],
                                filter_number = data[i,"filter_number"],
                                nlevels = data[i,"nlevels"])
    expect_equal(as.numeric(data[i,"vpr"][[1]]),
                 as.numeric(vpr))
    
  }
})


skip_on_cran()

library('fracdet')
data = readRDS("../testdata/vps_checks.RDS")

test_that("theoretical vpr works properly (intensive)", {
  for (i in seq_len(nrow(data))) {
    fd_vpr = theoreticalVpr(H = data[i,"H"],sigma2 = data[i,"sigma2"],
                            family = data[i,"family"],
                            filter_number = data[i,"filter.number"],
                            nlevels = data[i,"nlevels"])
    expect_equal(as.numeric(data[i,"vps"][[1]]),
                 as.numeric(fd_vpr))
  }
})


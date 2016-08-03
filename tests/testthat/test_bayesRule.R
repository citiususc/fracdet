skip_on_cran()

library('fracdet', quietly = TRUE)

test_that("bayes rule works as expected", {
  data = readRDS("../testdata/bayesRule.RDS")
  # tolerance has to be quite high due to the difficulties in simulating
  # the highest Hs
  result = apply(data,MARGIN = 1, function(x){
    bayesRule(x[["d"]], as.list(x[2:7]))
  })
  expect_equal(result, data$output, tolerance = 1e-5)
})




context("Test input to constructor")

test_that("Bad input results in failure",{

# example data from Phil
  data(dgp)
  gp13 <- 1*(gp==3)[gp %in% c(1, 3)]
  d13 <- d[gp %in% c(1, 3), gp %in% c(1, 3)]
  dummy <- 1:264
  # nothing supplied
  expect_error(gam(gp13 ~ s(dummy, bs="pco", k=15),
                   family="binomial"),
               "Please supply either a distance matrix or data and distance function!")
  # matrix and data but no function
  expect_error(gam(gp13 ~ s(dummy, bs="pco", k=15,
                            xt=list(D=d13, realdata=d13)),
                   family="binomial"),
               "Please supply either a distance matrix or data and distance function!")
  # matrix and data and function
  expect_error(gam(gp13 ~ s(dummy, bs="pco", k=15,
                            xt=list(D=d13, realdata=d13,dist_fn=dist)),
                   family="binomial"),
               "Please supply either a distance matrix or data and distance function!")
  # matrix and function
  expect_error(gam(gp13 ~ s(dummy, bs="pco", k=15,
                            xt=list(D=d13, dist_fn=dist)),
                   family="binomial"),
               "Please supply either a distance matrix or data and distance function!")


})

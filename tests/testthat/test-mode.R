

for(iter in 1:100) {
  x = runif(1, -100, 100)
  test_that("mode(fastest) deal with vector of constants - double", {
    expect_identical(vector_mode(rep(x, 100), min_value_count = 3, n_bins = 200, adjust = 1, kernel_width_in_sd = 3, na_mode = "check", bandwidth_method = "nrd_fastest"), x)
  })
  test_that("mode(fast) deal with vector of constants - double", {
    expect_identical(vector_mode(rep(x, 100), min_value_count = 3, n_bins = 200, adjust = 1, kernel_width_in_sd = 3, na_mode = "check", bandwidth_method = "nrd_fast"), x)
  })
  test_that("mode() deal with vector of constants - double", {
    expect_identical(vector_mode(rep(x, 100), min_value_count = 3, n_bins = 200, adjust = 1, kernel_width_in_sd = 3, na_mode = "check", bandwidth_method = "nrd"), x)
  })

  x = sample.int(200, size = 1) - 100
  test_that("mode(fastest) deal with vector of constants - integer", {
    expect_identical(vector_mode(rep(x, 100), min_value_count = 3, n_bins = 200, adjust = 1, kernel_width_in_sd = 3, na_mode = "check", bandwidth_method = "nrd_fastest"), x)
  })
  test_that("mode(fast) deal with vector of constants - integer", {
    expect_identical(vector_mode(rep(x, 100), min_value_count = 3, n_bins = 200, adjust = 1, kernel_width_in_sd = 3, na_mode = "check", bandwidth_method = "nrd_fast"), x)
  })
  test_that("mode() deal with vector of constants - integer", {
    expect_identical(vector_mode(rep(x, 100), min_value_count = 3, n_bins = 200, adjust = 1, kernel_width_in_sd = 3, na_mode = "check", bandwidth_method = "nrd"), x)
  })
}


# instead of comparing the x coordinate,
# which can differ between distributions (shapes) and density settings such as bandwidth,
# we here compare distance in height of the density.
# e.g. if X is a noisy normal distribution with a "flat top" then minor changes in KDE
# will affect x-coordinates, while the objective calling of the mode (top of distribution) might be very similar.


for(iter in 1:100) {
  for(withna in c(FALSE, TRUE)) {
    # generate test data; normal distribution with 10k data points
    mu = runif(1, min = -10, max = 10)
    x = rnorm(10000, mean = mu, sd = 1)
    if(withna) {
      x = c(NA, x, Inf, -Inf)
    }
    # save(x, file = "~/temp/x_testval_density_fail.RData")
    # load("~/temp/x_testval_density_fail.RData")

    # density function in base R
    d = stats::density(x, adjust = 0.9, bw = 'nrd', na.rm = TRUE)
    mode_baser = d$x[which.max(d$y)]

    mode_pairscale_fastest = vector_mode(x, adjust = 0.9, bandwidth_method = "nrd_fastest")
    mode_pairscale_fast = vector_mode(x, adjust = 0.9, bandwidth_method = "nrd_fast")
    mode_pairscale = vector_mode(x, adjust = 0.9, bandwidth_method = "nrd")

    density_xy_baser = list(x = mode_baser, y = max(d$y))
    density_xy_pairscale_fastest = approx(d$x, d$y, xout = mode_pairscale_fastest) # linear interpolation by x-coord
    density_xy_pairscale_fast = approx(d$x, d$y, xout = mode_pairscale_fast) # linear interpolation by x-coord
    density_xy_pairscale = approx(d$x, d$y, xout = mode_pairscale)

    # validate that there is less than 1% difference in density between where
    # base R functions calls the mode and where our Rcpp implementations do.
    y_diff_frac_fastest = abs(density_xy_baser$y - density_xy_pairscale_fastest$y) / density_xy_baser$y
    y_diff_frac_fast    = abs(density_xy_baser$y - density_xy_pairscale_fast$y) / density_xy_baser$y
    y_diff_frac         = abs(density_xy_baser$y - density_xy_pairscale$y) / density_xy_baser$y

    test_that(paste0("mode(fastest) random data iter=", iter), {
      expect(y_diff_frac_fastest < 0.01, failure_message = "y_diff_frac_fast < 0.01")
    })
    test_that(paste0("mode(fast) random data iter=", iter), {
      expect(y_diff_frac_fast < 0.01, failure_message = "y_diff_frac_fast < 0.01")
    })
    test_that(paste0("mode() random data iter=", iter), {
      expect(y_diff_frac < 0.01, failure_message = "y_diff_fract_robust < 0.01")
    })

  }

}

# # in case a test fails, use this code to plot the distribution and respective x-values
# plot(stats::density(x, adjust = 0.9, na.rm = TRUE))
# abline(v = c(mode_baser, mode_pairscale_fastest, mode_pairscale_fast, mode_pairscale), col = 1:4, lwd = 2, lty = c(1,2,2,2))
# c(mu, mode_baser, mode_pairscale_fastest, mode_pairscale_fast, mode_pairscale)






# r_density_mode = function(x) {
#   d = stats::density(x, na.rm = TRUE)
#   return(d$x[which.max(d$y)])
# }
#
# N = 1000
# for(iter in 1:50) {
#   for(xnoise in seq(from = 2, to = 4, by = 0.1)) {
#     x = c(stats::rnorm(N), stats::rnorm(floor(N * 0.1), mean = xnoise), stats::rnorm(N*0.2, mean=0, sd = 0.1))
#     # plot(density(x, adjust = 1))
#     x_sans_outliers = x[x > stats::quantile(x, probs = 0.025, na.rm = TRUE) & x < stats::quantile(x, probs = 0.9975, na.rm = TRUE)]
#     mode_ref = r_density_mode(x_sans_outliers)
#     mode = density_mode(x, minval = 1, n_points = 512, adjust = 1, max_dist = 3, na_mode = "check")
#     # relatively high tolerance because our implementation is expected to return "robust mode" so not exactly the same computation
#     expect_equal(all.equal(mode, mode_ref, tolerance = 0.1), TRUE)
#   }
# }


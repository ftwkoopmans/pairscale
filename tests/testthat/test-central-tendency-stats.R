
r_vector_madmean = function(x, threshold_std = 3) {
  x = x[is.finite(x)]
  if(length(x) <= 3) return (median(x))
  mean(x[ x >= median(x) - threshold_std * mad(x) & x <= median(x) + threshold_std * mad(x) ])
}



test_that("vector_trimmedmean() deal with Inf, test 1", {
  expect_identical(is.na(vector_trimmedmean(c(-Inf), min_value_count = 1, trim = 0.2, na_mode = "check")), TRUE)
})
test_that("vector_trimmedmean() deal with NA, test 1", {
  expect_identical(is.na(vector_trimmedmean(c(NA), min_value_count = 3, trim = 0.2, na_mode = "check")), TRUE)
})
test_that("vector_trimmedmean() deal with NA, test 2", {
  expect_identical(is.na(vector_trimmedmean(c(NA, NA), min_value_count = 3, trim = 0.2, na_mode = "check")), TRUE)
})
test_that("vector_trimmedmean() deal with NA, test 3", {
  expect_identical(is.na(vector_trimmedmean(c(NA, NA, 1), min_value_count = 3, trim = 0.2, na_mode = "check")), TRUE)
})
test_that("vector_trimmedmean() deal with NA, test 4", {
  expect_identical(is.na(vector_trimmedmean(c(NA, NA, 1, 2), min_value_count = 3, trim = 0.2, na_mode = "check")), TRUE)
})
test_that("vector_trimmedmean() deal with NA, test 5", {
  expect_identical(vector_trimmedmean(c(NA, NA, 1, 2, 3), min_value_count = 3, trim = 0.2, na_mode = "check"), 2)
})
test_that("vector_trimmedmean() deal with NA, test 6", {
  expect_identical(vector_trimmedmean(c(NA, NA, 2), min_value_count = 1, trim = 0.2, na_mode = "check"), 2)
})



test_that("vector_madmean() deal with Inf, test 1", {
  expect_identical(is.na(vector_madmean(c(-Inf), min_value_count = 1, na_mode = "check")), TRUE)
})
test_that("vector_madmean() deal with NA, test 1", {
  expect_identical(is.na(vector_madmean(c(NA), min_value_count = 3, na_mode = "check")), TRUE)
})
test_that("vector_madmean() deal with NA, test 2", {
  expect_identical(is.na(vector_madmean(c(NA, NA), min_value_count = 3, na_mode = "check")), TRUE)
})
test_that("vector_madmean() deal with NA, test 3", {
  expect_identical(is.na(vector_madmean(c(NA, NA, 1), min_value_count = 3, na_mode = "check")), TRUE)
})
test_that("vector_madmean() deal with NA, test 4", {
  expect_identical(is.na(vector_madmean(c(NA, NA, 1, 2), min_value_count = 3, na_mode = "check")), TRUE)
})
test_that("vector_madmean() deal with NA, test 5", {
  expect_identical(vector_madmean(c(NA, NA, 1, 2, 3), min_value_count = 3, na_mode = "check"), 2)
})
test_that("vector_madmean() deal with NA, test 6", {
  expect_identical(vector_madmean(c(NA, NA, 2), min_value_count = 1, na_mode = "check"), 2)
})



test_that("vector_median() deal with Inf, test 1", {
  expect_identical(is.na(vector_median(c(-Inf), min_value_count = 1)), TRUE)
})
test_that("vector_median() deal with NA, test 1", {
  expect_identical(is.na(vector_median(c(NA), min_value_count = 3)), TRUE)
})
test_that("vector_median() deal with NA, test 2", {
  expect_identical(is.na(vector_median(c(NA, NA), min_value_count = 3)), TRUE)
})
test_that("vector_median() deal with NA, test 3", {
  expect_identical(is.na(vector_median(c(NA, NA, 1), min_value_count = 3)), TRUE)
})
test_that("vector_median() deal with NA, test 4", {
  expect_identical(is.na(vector_median(c(NA, NA, 1, 2), min_value_count = 3)), TRUE)
})
test_that("vector_median() deal with NA, test 5", {
  expect_identical(vector_median(c(NA, NA, 1, 2, 3), min_value_count = 3), 2)
})
test_that("vector_median() deal with NA, test 6", {
  expect_identical(vector_median(c(NA, NA, 2), min_value_count = 1), 2)
})




test_that("vector_trimmedmean() deal with vector of constants", {
  expect_identical(vector_trimmedmean(rep(1,5), min_value_count = 3, trim = 0.2, na_mode = "check"), 1)
})
test_that("vector_trimmedmean() deal with vector of constants", {
  expect_identical(vector_madmean(rep(1,5), min_value_count = 3, na_mode = "check"), 1)
})
test_that("vector_trimmedmean() deal with vector of constants", {
  expect_identical(vector_median(rep(1,5), min_value_count = 3), 1)
})



# hardcoded testcase; x = c(0.4200354, 0.4798072, -1.2270806)

for(iter in 1:100) {
  # array with random values, N = 1..100
  x = rnorm(sample.int(100, 1))
  # try a range of trim settings, from 0 (no trim) to 0.45
  for(iter_trim in (0:9 / 20)) {
    test_that(paste0("vector_trimmedmean() random data + NA iter=", iter, " trim=", iter_trim), {
      expect_equal(vector_trimmedmean(c(x, NA), trim = iter_trim, na_mode = "check"), mean(x, trim = iter_trim), tolerance = 10^-6)
    })

    test_that(paste0("vector_trimmedmean() random data iter=", iter, " trim=", iter_trim), {
      expect_equal(vector_trimmedmean(x, trim = iter_trim, na_mode = "check"), mean(x, trim = iter_trim), tolerance = 10^-6)
    })

  }


  test_that(paste0("vector_madmean() random data iter=", iter), {
    expect_equal(vector_madmean(x, min_value_count = 1, threshold_std = 3, na_mode = "check"), r_vector_madmean(x, threshold_std = 3), tolerance = 10^-6)
  })

  test_that(paste0("vector_median() random data iter=", iter), {
    expect_equal(vector_median(x, min_value_count = 1), median(x, na.rm = TRUE), tolerance = 10^-6)
  })

}



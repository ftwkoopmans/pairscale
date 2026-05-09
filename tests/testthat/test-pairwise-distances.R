


# tricky matrix to normalize with "mode" approach, testing the density function
# in this example, the ratio between columns can be a constant value, or zero
x = cbind(1:5, 1:5 + 0.1, 1:5 + 0.1, 1:5 - 3, 1:5 - 3.3)
x_expect = rbind(
  c(0.0,-0.1,-0.1,3.0,3.3),
  c(0.1,0.0,0.0,3.1,3.4),
  c(0.1,0.0,0.0,3.1,3.4),
  c(-3.0,-3.1,-3.1,0.0,0.3),
  c(-3.3,-3.4,-3.4,-0.3,0.0)
)



for(colspec in list(1:2, 3:4, c(5,2))) {

  test_that(paste0("hardcoded data, subset of columns ", colspec[1],",",colspec[2]," ; pairdiff_median() vs expected values"), {
    expect_equal(pairdiff_median(x, cols = colspec, min_value_count = 3, na_mode = "present"), x_expect[colspec,colspec])
  })

  test_that(paste0("hardcoded data, subset of columns ", colspec[1],",",colspec[2]," ; pairdiff_mean() vs expected values"), {
    expect_equal(pairdiff_mean(x, cols = colspec, min_value_count = 3, na_mode = "present"), x_expect[colspec,colspec])
  })

  test_that(paste0("hardcoded data, subset of columns ", colspec[1],",",colspec[2]," ; pairdiff_trimmedmean() vs expected values"), {
    expect_equal(pairdiff_trimmedmean(x, cols = colspec, min_value_count = 3, na_mode = "present"), x_expect[colspec,colspec])
  })

  test_that(paste0("hardcoded data, subset of columns ", colspec[1],",",colspec[2]," ; pairdiff_madmean() vs expected values"), {
    expect_equal(pairdiff_madmean(x, cols = colspec, min_value_count = 3, threshold_std = 3, na_mode = "present"), x_expect[colspec,colspec])
  })


  test_that(paste0("hardcoded data, subset of columns ", colspec[1],",",colspec[2]," ; pairdiff_mode(fastest, na=present) vs expected values"), {
    expect_equal(pairdiff_mode(x, cols = colspec, min_value_count = 3, na_mode = "present", n_bin = 100, adjust = 1, bandwidth_method = "nrd_fastest"), x_expect[colspec,colspec])
  })

  test_that(paste0("hardcoded data, subset of columns ", colspec[1],",",colspec[2]," ; pairdiff_mode(fastest, na=unchecked) vs expected values"), {
    expect_equal(pairdiff_mode(x, cols = colspec, min_value_count = 3, na_mode = "unchecked", n_bin = 100, adjust = 1, bandwidth_method = "nrd_fastest"), x_expect[colspec,colspec])
  })


  test_that(paste0("hardcoded data, subset of columns ", colspec[1],",",colspec[2]," ; pairdiff_mode(fast, na=present) vs expected values"), {
    expect_equal(pairdiff_mode(x, cols = colspec, min_value_count = 3, na_mode = "present", n_bin = 100, adjust = 1, bandwidth_method = "nrd_fast"), x_expect[colspec,colspec])
  })

  test_that(paste0("hardcoded data, subset of columns ", colspec[1],",",colspec[2]," ; pairdiff_mode(fast, na=unchecked) vs expected values"), {
    expect_equal(pairdiff_mode(x, cols = colspec, min_value_count = 3, na_mode = "unchecked", n_bin = 100, adjust = 1, bandwidth_method = "nrd_fast"), x_expect[colspec,colspec])
  })


  test_that(paste0("hardcoded data, subset of columns ", colspec[1],",",colspec[2]," ; pairdiff_mode(na=present) vs expected values"), {
    expect_equal(pairdiff_mode(x, cols = colspec, min_value_count = 3, na_mode = "present", n_bin = 100, adjust = 1, bandwidth_method = "nrd"), x_expect[colspec,colspec])
  })

  test_that(paste0("hardcoded data, subset of columns ", colspec[1],",",colspec[2]," ; pairdiff_mode(na=unchecked) vs expected values"), {
    expect_equal(pairdiff_mode(x, cols = colspec, min_value_count = 3, na_mode = "unchecked", n_bin = 100, adjust = 1, bandwidth_method = "nrd"), x_expect[colspec,colspec])
  })

}







### min_value_count=5 should remove pairwise scalings for sample 1 in this example
x = cbind(c(1:4, NA), 1:5 + 0.1, 1:5 + 0.1, 1:5 - 3, 1:5 - 3.3)
x_expect_na = matrix(FALSE, nrow = 5, ncol = 5)
x_expect_na[1,2:5] = TRUE
x_expect_na[2:5,1] = TRUE

test_that("hardcoded data, pairdiff_median() - min_value_count check", {
  expect_equal(is.na(pairdiff_median(x, min_value_count = 5, na_mode = "present")), x_expect_na)
})

test_that("hardcoded data, pairdiff_median() - min_value_count check", {
  expect_equal(is.na(pairdiff_mean(x, min_value_count = 5, na_mode = "present")), x_expect_na)
})

test_that("hardcoded data, pairdiff_trimmedmean() - min_value_count check", {
  expect_equal(is.na(pairdiff_trimmedmean(x, min_value_count = 5, na_mode = "present")), x_expect_na)
})

test_that("hardcoded data, pairdiff_madmean() - min_value_count check", {
  expect_equal(is.na(pairdiff_madmean(x, min_value_count = 5, threshold_std = 3, na_mode = "present")), x_expect_na)
})

test_that("hardcoded data, pairdiff_mode(fastest) - min_value_count check 1", {
  expect_equal(is.na(pairdiff_mode(x, min_value_count = 5, na_mode = "present", n_bins = 100, bandwidth_method = "nrd_fastest")), x_expect_na)
})

test_that("hardcoded data, pairdiff_mode(fast) - min_value_count check 1", {
  expect_equal(is.na(pairdiff_mode(x, min_value_count = 5, na_mode = "present", n_bins = 100, bandwidth_method = "nrd_fast")), x_expect_na)
})

test_that("hardcoded data, pairdiff_mode() - min_value_count check 1", {
  expect_equal(is.na(pairdiff_mode(x, min_value_count = 5, na_mode = "present", n_bins = 100, bandwidth_method = "nrd")), x_expect_na)
})








# with minval filtering, we can artifically create 2 blocks that have no overlap (considering min_value_count)
# -->> test that we get finite scaling factors and there is no runaway
x = cbind(c(1:4, NA), c(1:3 + 0.1, NA, NA), c(NA, NA, 1:3 + 0.1), c(NA, NA, 1:3 + 0.3))
x_expect = matrix(c(0,-0.1,NA,NA, 0.1,0,NA,NA, NA,NA,0,-0.2, NA,NA,0.2,0), ncol=4, byrow = T)

test_that("hardcoded data, partial column overlap, pairdiff_median()", {
  expect_equal(pairdiff_median(x, min_value_count = 3, na_mode = "present"), x_expect)
})

test_that("hardcoded data, partial column overlap, pairdiff_median()", {
  expect_equal(pairdiff_mean(x, min_value_count = 3, na_mode = "present"), x_expect)
})

test_that("hardcoded data, partial column overlap, pairdiff_trimmedmean()", {
  expect_equal(pairdiff_trimmedmean(x, min_value_count = 3, na_mode = "present"), x_expect)
})

test_that("hardcoded data, partial column overlap, pairdiff_madmean()", {
  expect_equal(pairdiff_madmean(x, min_value_count = 3, threshold_std = 3, na_mode = "present"), x_expect)
})

test_that("hardcoded data, partial column overlap, pairdiff_mode(fastest)", {
  expect_equal(pairdiff_mode(x, min_value_count = 3, na_mode = "present", n_bin = 100, adjust = 1, bandwidth_method = "nrd_fastest"), x_expect)
})

test_that("hardcoded data, partial column overlap, pairdiff_mode(fast)", {
  expect_equal(pairdiff_mode(x, min_value_count = 3, na_mode = "present", n_bin = 100, adjust = 1, bandwidth_method = "nrd_fast"), x_expect)
})

test_that("hardcoded data, partial column overlap, pairdiff_mode()", {
  expect_equal(pairdiff_mode(x, min_value_count = 3, na_mode = "present", n_bin = 100, adjust = 1, bandwidth_method = "nrd"), x_expect)
})



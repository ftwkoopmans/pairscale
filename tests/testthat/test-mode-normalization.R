
# data: apply constant scaling factor to each column
# post normalization, all values on each rows should be the same, i.e. the diff computed here should be zero
for(iter in 1:10) {

  tmp = rnorm(100)
  x = cbind(tmp, tmp + 3, tmp - 1)
  l = list()

  y = x + 1  # make a deep copy / copy-on-modify
  s = pairscale_mean(y)
  l$pairscale_mean = y

  y = x + 1
  s = pairscale_median(y)
  l$pairscale_median = y

  y = x + 1
  s = pairscale_madmean(y)
  l$pairscale_madmean = y

  y = x + 1
  s = pairscale_trimmedmean(y)
  l$pairscale_trimmedmean = y

  y = x + 1
  s = pairscale_mode(y, bandwidth_method = "nrd_fastest")
  l$pairscale_mode_fastest = y

  y = x + 1
  s = pairscale_mode(y, bandwidth_method = "nrd_fast")
  l$pairscale_mode_fast = y

  y = x + 1
  s = pairscale_mode(y, bandwidth_method = "nrd")
  l$pairscale_mode = y
  # boxplot(x + 1, main = "before")
  # boxplot(y, main = "after - mode should be equal, median may not be !")


  for(function_name in names(l)) {
    y = l[[function_name]]

    test_that(paste0("3 column fixed scaling, ", function_name, ", columns 1,2"), {
      expect_equal(y[,1], y[,2], tolerance = 0.001)
    })
    test_that(paste0("3 column fixed scaling, ", function_name, ", columns 1,3"), {
      expect_equal(y[,1], y[,3], tolerance = 0.001)
    })
  }

  rm(tmp, x, l, s)
}


r_density_mode = function(x) {
  y = na.omit(x)
  d = stats::density(y[y > quantile(y, probs=0.05) & y < quantile(y, probs=0.95)], bw = "nrd", adjust = 0.9)
  return(d$x[which.max(d$y)])
}

for(iter in 1:10) {

  l = list()
  tmp = rnorm(1000)
  tmp1 = c(tmp + runif(length(tmp), min = -0.1, max = 0.2), rnorm(120, mean = 2, sd = 1))
  tmp2 = c(tmp + runif(length(tmp), min = -0.1, max = 0.2), rnorm(120, mean = 2, sd = 1))
  x = cbind(tmp1, tmp1 + 3, tmp1 - 1, tmp2 + 3, tmp2 - 2, tmp2 + 1)
  x_clusters = c(1,1,1, 9,9,9)
  # plot(density(tmp1)); lines(density(tmp2), col=2)
  # plot(density(x[,1] - x[,5]))


  y = x + 1  # make a deep copy / copy-on-modify
  s = pairscale_mean(y)
  l$pairscale_mean = y

  y = x + 1
  s = pairscale_median(y)
  l$pairscale_median = y

  y = x + 1
  s = pairscale_madmean(y)
  l$pairscale_madmean = y

  y = x + 1
  s = pairscale_trimmedmean(y)
  l$pairscale_trimmedmean = y

  y = x + 1
  s = pairscale_mode(y, bandwidth_method = "nrd_fastest")
  l$pairscale_mode_fastest = y

  y = x + 1
  s = pairscale_mode(y, bandwidth_method = "nrd_fast")
  l$pairscale_mode_fast = y

  y = x + 1
  s = pairscale_mode(y, bandwidth_method = "nrd")
  l$pairscale_mode = y
  # boxplot(x + 1, main = "before")
  # boxplot(y, main = "after - mode should be equal, median may not be !")


  # mean scaling yields errors on some random datasets, it's very sensitive to outliers
  # we do run it while unit testing to ensure it doesn't generate errors
  y = x + 1
  s = pairscale_mean(y, clusters = x_clusters)
  # l$pairscale_mean_withclusters = y

  y = x + 1
  s = pairscale_median(y, clusters = x_clusters)
  l$pairscale_median_withclusters = y

  y = x + 1
  s = pairscale_madmean(y, clusters = x_clusters)
  l$pairscale_madmean_withclusters = y

  y = x + 1
  s = pairscale_trimmedmean(y, clusters = x_clusters)
  l$pairscale_trimmedmean_withclusters = y

  y = x + 1
  s = pairscale_mode(y, clusters = x_clusters, bandwidth_method = "nrd_fastest")
  l$pairscale_mode_fastest_withclusters = y

  y = x + 1
  s = pairscale_mode(y, clusters = x_clusters, bandwidth_method = "nrd_fast")
  l$pairscale_mode_fast_withclusters = y

  y = x + 1
  s = pairscale_mode(y, clusters = x_clusters, bandwidth_method = "nrd")
  l$pairscale_mode_withclusters = y


  for(function_name in names(l)) {
    y = l[[function_name]]
    test_that(paste0("minor noise, ", function_name, ", columns 1,2"), {
      expect_equal(y[,1], y[,2], tolerance = 0.001)
    })
    test_that(paste0("minor noise, ", function_name, ", columns 1,3"), {
      expect_equal(y[,1], y[,3], tolerance = 0.001)
    })
    test_that(paste0("minor noise, ", function_name, ", columns 4,5"), {
      expect_equal(y[,4], y[,5], tolerance = 0.001)
    })
    test_that(paste0("minor noise, ", function_name, ", zero mode in 1,4 post normalization"), {
      expect_equal(r_density_mode(y[,1] - y[,4]), 0, tolerance = 0.1)
    })

  }



  rm(tmp, tmp1, tmp2, x, s)
}


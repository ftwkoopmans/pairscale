// Copyright (C) 2026 Frank Koopmans
//
// This file is part of pairscale.
//
// pairscale is free software: you can redistribute it and/or modify
// it under the terms of the GNU Affero General Public License as published
// by the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.
//
// pairscale is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
// GNU Affero General Public License for more details.
//
// You should have received a copy of the GNU Affero General Public License
// along with this program. If not, see <https://www.gnu.org/licenses/>.

#include <RcppArmadillo.h>
#include "normalize.h"
#include "pairwise_distance.h"
#include "central_tendency.h"
#include <RcppArmadillo.h>

// [[Rcpp::depends(RcppArmadillo)]]


// upstream we already validated that ncol_input > 1
// if no cols are provided, return a vector 0:(ncol_input-1)
// else, validate input (all cols >=0 and <= ncol_input-1)
arma::uvec validate_column_spec_from_r__return_cpp_indices(int ncol_input, Rcpp::Nullable<Rcpp::IntegerVector> cols) {
  if(cols.isUsable()) {
    // typecast. There must be a cleaner way to go from nullable int vector to arma::uvec ...
    Rcpp::IntegerVector r_cols_vector(cols);
    int cols_len = r_cols_vector.length();
    arma::uvec r_cols(cols_len, arma::fill::none); // holds 1-indexed indices, R style
    for(int i = 0; i < cols_len; i++) {
      r_cols(i) = r_cols_vector(i);
    }

    //// enforce that cols are valid integer column indices in X without duplicates
    // calling this from R, NA/NaN/Inf/-Inf may be converted to a tiny integer number (-2147483648) into the arma::ivec data type !!
    // thus, has_nan and has_inf do not work on ivec !!
    if(r_cols.has_nan() || r_cols.has_inf()) {
      throw std::invalid_argument("'cols' must not contain NA or infinite values");
    }
    arma::uvec col_ui = arma::find_unique(r_cols);
    if(col_ui.n_elem != r_cols.n_elem) {
      throw std::invalid_argument("'cols' may not contain duplicates");
    }
    if((int)(r_cols.min()) < 1 || (int)(r_cols.max()) > ncol_input) {
      throw std::invalid_argument("'cols' must contain integer values, NA/NaN/Infinite values not allowed. R-style indexes (1-based) are expected, so the smallest accepted value is 1 and the largest is ncol(x)");
    }
    if(r_cols.n_elem < 2) {
      throw std::invalid_argument("'cols' must contain at least 2 values");
    }

    // result: zero-based indices
    return r_cols - 1;
  }

  // no cols provided; return a sequence along columns
  // result: zero-based indices
  arma::uvec res(ncol_input, arma::fill::none);
  for(int i = 0; i < ncol_input; i++) {
    res(i) = i;
  }
  return res;
}



void validate_matrix_nonempty(const arma::mat& x) {
  if(x.n_rows < 2 || x.n_cols < 2) {
    throw std::invalid_argument("data matrix must have at least 2 columns and rows");
  }
}



bool validate_namode_vec(const arma::vec& x, std::string na_mode) {
  if(na_mode == "unchecked") return false;
  if(na_mode == "present") return true;
  if(na_mode == "check") return x.has_nonfinite();
  throw std::invalid_argument("'na_mode' must be either of; unchecked , present , check");
}



bool validate_namode_mat(const arma::mat& x, std::string na_mode) {
  if(na_mode == "unchecked") return false;
  if(na_mode == "present") return true;
  if(na_mode == "check") return x.has_nonfinite();
  throw std::invalid_argument("'na_mode' must be either of; unchecked , present , check");
}



void validate_mode_bwmethod(std::string s) {
  if(s == "nrd") return;
  if(s == "nrd_fast") return;
  if(s == "nrd_fastest") return;
  throw std::invalid_argument("'bandwidth_method' must be either of; nrd , nrd_fast , nrd_fastest");
}


void validate_trim(double trim) {
  if(!std::isfinite(trim) || trim < 0 || trim >= 0.5) {
    throw std::invalid_argument("'trim' must be a finite number >= 0 and < 0.5. At/near 0.5 you will remove all data (~50% from low- and high-end) ! Typical values are 0.1~0.3");
  }
}


bool validate_minvaluecount(int max, int val) {
  if(!std::isfinite(val) || val <= 0) {
    throw std::invalid_argument("'min_value_count' must be a finite integer >= 1");
  }
  return val > max;
}


void validate_density_nbin(int n) {
  if(!std::isfinite(n) || n < 25 || n > 10000) {
    throw std::invalid_argument("'n_bin' parameter for density computation must be a finite integer >= 25 and <= 10000");
  }
}


void validate_density_adjust(double val) {
  if(!std::isfinite(val) || val < 0.05 || val > 5) {
    throw std::invalid_argument("'adjust' parameter for density computation must be a finite number >= 0.05 and <= 5");
  }
}


void validate_density_kernelwidth(double val) {
  if(!std::isfinite(val) || val < 1 || val > 6) {
    throw std::invalid_argument("'kernel_width_in_sd' parameter for density computation must be a finite number >= 1 and <= 6");
  }
}


void validate_density_fracmaxdens(double val) {
  if(!std::isfinite(val) || val < 0.3 || val > 1) {
    throw std::invalid_argument("'mode_frac_maxdens' parameter for density computation must be a finite number >= 0.3 and <= 1");
  }
}


void validate_madmeanthreshold(double val) {
  if(!std::isfinite(val) || val < 0.1 || val > 50) {
    throw std::invalid_argument("'threshold_std' parameter for MAD-trimmed mean must be a finite number >= 0.1 and <= 50");
  }
}


void validate_niterirls(double val) {
  if(!std::isfinite(val) || val < 0) {
    throw std::invalid_argument("'niter_irls' parameter must be an integer >= 0");
  }
}


arma::uvec validate_clusters(int ncol_input, Rcpp::Nullable<Rcpp::IntegerVector> clusters) {
  arma::uvec clusters_cpp(ncol_input, arma::fill::none);

  // typecast. There must be a cleaner way to go from nullable int vector to arma::uvec ...
  if(clusters.isUsable()) {
    Rcpp::IntegerVector clusters_r(clusters);
    int clusters_len = clusters_r.length();
    if(clusters_len != ncol_input) {
      throw std::invalid_argument("'clusters' must be NULL, indicating there are no clusters / sample groups, or a vector with positive integers of same length as the number of columns in x");
    }
    for(int i = 0; i < ncol_input; i++) {
      clusters_cpp(i) = clusters_r(i);
    }

    // importantly, setting NA to an Armadillo uvec will autocast to value 0
    // so we don't accept zero integers as clusters ID, this requirement will guard against such edge cases (e.g. c(NA, 0, 0, 1, 1))
    arma::uvec uclust = arma::unique(clusters_cpp);
    if(uclust.has_nan() || uclust.has_inf() || uclust.min() <= 0) {
      throw std::invalid_argument("'clusters' must contain only positive integer values (not NA, not negative, not zero)");
    }

  } else {
    // no input, default; fill with a constant
    for(int i = 0; i < ncol_input; i++) {
      clusters_cpp(i) = 1;
    }
  }

  return clusters_cpp;
}







//' @title Compute mean value of a vector, with optional filtering for N datapoints
//' @description Compute measure of central tendency using efficient C++ code
//' @name vector_mean
//' @inheritParams vector_mode
//' @inheritParams pairscale_mode
//' @return a single numeric value representing the mean
//' @export
// [[Rcpp::export]]
double vector_mean(const arma::vec& x, int min_value_count = 1, std::string na_mode = "check") {
  if(validate_minvaluecount(x.n_elem, min_value_count)) {
    return arma::datum::nan;
  }
  bool check_na = validate_namode_vec(x, na_mode);
  arma::vec buffer = x;
  return _mean_vector(buffer, min_value_count, check_na);
}



//' @title Compute median value of a vector, with optional filtering for N datapoints
//' @description Compute measure of central tendency using efficient C++ code
//' @name vector_median
//' @inheritParams vector_mode
//' @inheritParams pairscale_mode
//' @return a single numeric value representing the median
//' @export
// [[Rcpp::export]]
double vector_median(const arma::vec& x, int min_value_count = 1, std::string na_mode = "check") {
  if(validate_minvaluecount(x.n_elem, min_value_count)) {
    return arma::datum::nan;
  }
  bool check_na = validate_namode_vec(x, na_mode);
  arma::vec buffer = x;
  return _median_vector(buffer, min_value_count, check_na);
}



//' @title Compute trimmed mean value of a vector, with optional filtering for N datapoints
//' @description Compute measure of central tendency using efficient C++ code
//' @name vector_trimmedmean
//' @inheritParams vector_mode
//' @inheritParams pairscale_mode
//' @inheritParams pairscale_trimmedmean
//' @return a single numeric value representing the trimmed mean
//' @export
// [[Rcpp::export]]
double vector_trimmedmean(const arma::vec& x, int min_value_count = 1, double trim = 0.2, std::string na_mode = "check") {
  if(validate_minvaluecount(x.n_elem, min_value_count)) {
    return arma::datum::nan;
  }
  bool check_na = validate_namode_vec(x, na_mode);
  validate_trim(trim);
  arma::vec buffer = x;
  return _trimmedmean_vector(buffer, min_value_count, check_na, trim);
}



//' @title Compute MAD-trimmed mean value of a vector, with optional filtering for N datapoints
//' @description Compute measure of central tendency using efficient C++ code
//' @name vector_madmean
//' @inheritParams vector_mode
//' @inheritParams pairscale_mode
//' @inheritParams pairscale_madmean
//' @return a single numeric value representing the MAD-trimmed mean
//' @export
// [[Rcpp::export]]
double vector_madmean(const arma::vec& x, int min_value_count = 1, double threshold_std = 3, std::string na_mode = "check") {
  if(validate_minvaluecount(x.n_elem, min_value_count)) {
    return arma::datum::nan;
  }
  validate_madmeanthreshold(threshold_std);
  bool check_na = validate_namode_vec(x, na_mode);
  arma::vec buffer = x;
  arma::vec buffer2(x.n_elem, arma::fill::none);
  return _madmean_vector(buffer, buffer2, min_value_count, check_na, threshold_std);
}



//' @title Compute mode of a vector, with optional filtering for N datapoints
//' @description Compute measure of central tendency using efficient C++ code
//' @name vector_mode
//' @inheritParams pairscale_mode
//' @param x numeric input vector, may contain non-finite values (removed if `na_mode` is left to default)
//' @return a single numeric value representing the mode
//' @export
// [[Rcpp::export]]
double vector_mode(const arma::vec& x, int min_value_count = 3, int n_bins = 512, double adjust = 1, double kernel_width_in_sd = 3, std::string bandwidth_method = "nrd", double mode_frac_maxdens = 1, std::string na_mode = "check") {
  if(validate_minvaluecount(x.n_elem, min_value_count)) {
    return arma::datum::nan;
  }
  validate_density_nbin(n_bins);
  validate_density_adjust(adjust);
  validate_density_kernelwidth(kernel_width_in_sd);
  validate_density_fracmaxdens(mode_frac_maxdens);
  validate_mode_bwmethod(bandwidth_method);
  bool check_na = validate_namode_vec(x, na_mode);
  arma::vec buffer = x;
  return _mode_vector(buffer, min_value_count, check_na, adjust, kernel_width_in_sd, n_bins, bandwidth_method, mode_frac_maxdens);
}





//' @title Distance matrix between columns of a matrix using their mean differences
//' @description Pairwise differences between all columns in a matrix
//' @name pairdiff_mean
//' @inheritParams pairscale_mode
//' @inheritParams pairdiff_mode
//' @return a N x N numeric matrix (where N is number of column in input `x`) representing the mean difference between each column
//' @export
// [[Rcpp::export]]
arma::mat pairdiff_mean(const arma::mat& x, Rcpp::Nullable<Rcpp::IntegerVector> cols = R_NilValue, int min_value_count = 3, std::string na_mode = "check") {
  validate_minvaluecount(x.n_rows, min_value_count);
  bool check_na = validate_namode_mat(x, na_mode);
  validate_matrix_nonempty(x);
  return _pairwise_distance_mean(x, validate_column_spec_from_r__return_cpp_indices(x.n_cols, cols), min_value_count, check_na);
}



//' @title Distance matrix between columns of a matrix using their trimmed median differences
//' @description Pairwise differences between all columns in a matrix
//' @name pairdiff_median
//' @inheritParams pairscale_mode
//' @inheritParams pairdiff_mode
//' @return a N x N numeric matrix (where N is number of column in input `x`) representing the median difference between each column
//' @export
// [[Rcpp::export]]
arma::mat pairdiff_median(const arma::mat& x, Rcpp::Nullable<Rcpp::IntegerVector> cols = R_NilValue, int min_value_count = 3, std::string na_mode = "check") {
  validate_minvaluecount(x.n_rows, min_value_count);
  bool check_na = validate_namode_mat(x, na_mode);
  validate_matrix_nonempty(x);
  return _pairwise_distance_median(x, validate_column_spec_from_r__return_cpp_indices(x.n_cols, cols), min_value_count, check_na);
}



//' @title Distance matrix between columns of a matrix using their trimmed mean differences
//' @description Pairwise differences between all columns in a matrix
//' @name pairdiff_trimmedmean
//' @inheritParams pairscale_mode
//' @inheritParams pairdiff_mode
//' @inheritParams pairscale_trimmedmean
//' @return a N x N numeric matrix (where N is number of column in input `x`) representing the trimmed mean difference between each column
//' @export
// [[Rcpp::export]]
arma::mat pairdiff_trimmedmean(const arma::mat& x, Rcpp::Nullable<Rcpp::IntegerVector> cols = R_NilValue, int min_value_count = 3, double trim = 0.2, std::string na_mode = "check") {
  validate_minvaluecount(x.n_rows, min_value_count);
  bool check_na = validate_namode_mat(x, na_mode);
  validate_trim(trim);
  validate_matrix_nonempty(x);
  return _pairwise_distance_trimmedmean(x, validate_column_spec_from_r__return_cpp_indices(x.n_cols, cols), min_value_count, trim, check_na);
}



//' @title Distance matrix between columns of a matrix using their MAD-trimmed mean differences
//' @description Pairwise differences between all columns in a matrix
//' @name pairdiff_madmean
//' @inheritParams pairscale_mode
//' @inheritParams pairdiff_mode
//' @inheritParams pairscale_madmean
//' @return a N x N numeric matrix (where N is number of column in input `x`) representing the MAD-trimmed mean difference between each column
//' @export
// [[Rcpp::export]]
arma::mat pairdiff_madmean(const arma::mat& x, Rcpp::Nullable<Rcpp::IntegerVector> cols = R_NilValue, int min_value_count = 3, double threshold_std = 3, std::string na_mode = "check") {
  validate_minvaluecount(x.n_rows, min_value_count);
  validate_madmeanthreshold(threshold_std);
  bool check_na = validate_namode_mat(x, na_mode);
  validate_matrix_nonempty(x);
  return _pairwise_distance_madmean(x, validate_column_spec_from_r__return_cpp_indices(x.n_cols, cols), min_value_count, threshold_std, check_na);
}



//' @title Distance matrix between columns of a matrix using their mode differences
//' @description Pairwise differences between all columns in a matrix
//' @name pairdiff_mode
//' @inheritParams pairscale_mode
//' @param cols optionally, provide an integer vector with column indices (in `x`) that should be used (these should be 1-based indices as per usual in R). Or set to `NULL` or `integer()` to use all columns
//' @return a N x N numeric matrix (where N is number of column in input `x`) representing the mode difference between each column
//' @export
// [[Rcpp::export]]
arma::mat pairdiff_mode(const arma::mat& x, Rcpp::Nullable<Rcpp::IntegerVector> cols = R_NilValue, int min_value_count = 3, int n_bins = 512, double adjust = 1, double kernel_width_in_sd = 3, std::string bandwidth_method = "nrd", double mode_frac_maxdens = 1, std::string na_mode = "check") {
  validate_minvaluecount(x.n_rows, min_value_count);
  validate_density_nbin(n_bins);
  validate_density_adjust(adjust);
  validate_density_kernelwidth(kernel_width_in_sd);
  validate_density_fracmaxdens(mode_frac_maxdens);
  validate_mode_bwmethod(bandwidth_method);
  bool check_na = validate_namode_mat(x, na_mode);
  validate_matrix_nonempty(x);
  return _pairwise_distance_mode(x, validate_column_spec_from_r__return_cpp_indices(x.n_cols, cols), min_value_count, n_bins, adjust, kernel_width_in_sd, bandwidth_method, mode_frac_maxdens, check_na);
}






//' @title Normalize matrix columns using their mode differences
//' @description Pairwise normalization of columns in a matrix, using the mode to define pairwise distances between columns
//' @name pairscale_mode
//' @param x a numeric matrix that needs to be normalized. This variable is changed by reference! i.e. after this function, the original variable is updated
//' @param clusters an integer vector, of same length as `ncol(x)`, that describes cluster identifiers. These values must describe a continuous set of integers between 1 and the total number of unique clusters. For example, a 6 column matrix with 3 columns from some experimental condition followed by 3 columns from another condition would be `c(1,1,1, 2,2,2)`. To disable mode-between rescaling, i.e. treat the entire matrix as one cluster, provide a vector with ones.
//' @param min_value_count the minimum number of values overlapping between a pair of columns (if there are fewer overlapping values, respective rescaling is not computed)
//' @param n_bins grid size for density computation (the resolution / number of data points to use for binning column diffs)
//' @param adjust bandwidth adjust factor for density computation
//' @param kernel_width_in_sd maximum distance in standard deviations at which we'll include data points for the Gaussian kernal. Typically 3 or 4
//' @param na_mode string value that indicated how should we should deal with NA values (default). "check" = test if NA values are present and if so, remove these. "present" = you already know NA values are present so we can skip the check for NA values for efficiency. "unchecked" = you guarantee that there are no NA values in input, thus we'll use the fastest code paths that skip over any NA checks downstream; ONLY select this if you are sure there are no NA/NaN/Inf/-Inf values !
//' @param bandwidth_method method in which this function computes bandwidth and optionally trims the data prior to binning. "nrd" is the robust, safe default. "nrd_fast" is faster and yields similar results for most distributions. Use "nrd_fastest" only when all pairwise distances are known to be near gaussian (i.e. no strong outliers and sd() is a reliable metric). "nrd_subset" is an experimental option that may be removed, it is fast but heavily favors symmetric distributions and is thus biased !
//'
//'   Valid options:
//'   * `"nrd_fastest"` ; use all data points, for bandwidth don't use IQR. Sensitive to outliers !
//'   * `"nrd_fast"` ; trim to quantiles 0.05 and 0.95, then use sd over these data points to compute bandwidth. Not as robust as the original nrd method, which also considers IQR, but much faster.
//'   * `"nrd"` ; adapted version of nrd; instead of IQR and sd we use IQR and sd over data trimmed at quantiles 0.05 and 0.95. Since we already sorted the data, computing sd on trimmed subset is 'free'
//' @param mode_frac_maxdens set to 1 to return the x-coordinate where the density is highest (mode). Setting this to a value < 1 will make this function compute not the mode, but the mean (x) value of the density where the density is some fraction higher than the maximum density. Typical value; 1. Optionally, set to 0.9 or 0.8 for possibly more robust center-finding, depending in your data distribution. Must not be smaller than 0.1 but recommended to never be lower than 0.7
//' @param niter_irls number of iterations to increase robustness (set to zero to disable)
//' @return a numeric vector that represents the normalization factors that were applied to each column in x. Note that x is updated by reference.
//' @export
// [[Rcpp::export]]
Rcpp::NumericVector pairscale_mode(arma::mat& x, Rcpp::Nullable<Rcpp::IntegerVector> clusters = R_NilValue, int min_value_count = 3, int n_bins = 512, double adjust = 1, double kernel_width_in_sd = 3, std::string bandwidth_method = "nrd", double mode_frac_maxdens = 1, int niter_irls = 50, std::string na_mode = "check") {
  const int Ncols = x.n_cols;
  validate_minvaluecount(x.n_rows, min_value_count);
  validate_density_nbin(n_bins);
  validate_density_adjust(adjust);
  validate_density_kernelwidth(kernel_width_in_sd);
  validate_density_fracmaxdens(mode_frac_maxdens);
  validate_mode_bwmethod(bandwidth_method);
  arma::uvec clusters_cpp = validate_clusters(Ncols, clusters);
  validate_niterirls(niter_irls);
  bool check_na = validate_namode_mat(x, na_mode);

  double tmean_trim = 0.2; // default, not used @ mode
  double threshold_std = 3;
  Rcpp::NumericVector result = Rcpp::wrap(_pairscale_normalization(x, clusters_cpp, "mode", min_value_count, n_bins, adjust, kernel_width_in_sd, bandwidth_method, mode_frac_maxdens, tmean_trim, threshold_std, niter_irls, check_na));
  result.attr("dim") = R_NilValue; // strip the dim attribute so the result won't be a 'column vector' aka matrix with 1 column
  return result;
}



//' @title Normalize matrix columns using their median differences
//' @description Pairwise normalization of columns in a matrix, using the median to define pairwise distances between columns
//' @name pairscale_median
//' @inheritParams pairscale_mode
//' @return a numeric vector that represents the normalization factors that were applied to each column in x. Note that x is updated by reference.
//' @export
// [[Rcpp::export]]
Rcpp::NumericVector pairscale_median(arma::mat& x, Rcpp::Nullable<Rcpp::IntegerVector> clusters = R_NilValue, int min_value_count = 3, int niter_irls = 50, std::string na_mode = "check") {
  const int Ncols = x.n_cols;
  validate_minvaluecount(x.n_rows, min_value_count);
  bool check_na = validate_namode_mat(x, na_mode);
  arma::uvec clusters_cpp = validate_clusters(Ncols, clusters);
  validate_niterirls(niter_irls);

  // default, not used @ median
  double tmean_trim = 0.2, density_adjust = 1, density_kernel_width_in_sd = 3, density_mode_frac_maxdens = 1, threshold_std = 3;
  int density_npoints = 200;
  std::string density_bandwidth_method = "nrd";
  Rcpp::NumericVector result = Rcpp::wrap(_pairscale_normalization(x, clusters_cpp, "median", min_value_count, density_npoints, density_adjust, density_kernel_width_in_sd, density_bandwidth_method, density_mode_frac_maxdens, tmean_trim, threshold_std, niter_irls, check_na));
  result.attr("dim") = R_NilValue;
  return result;
}



//' @title Normalize matrix columns using their mean differences
//' @description Pairwise normalization of columns in a matrix, using the mean to define pairwise distances between columns
//' @name pairscale_mean
//' @inheritParams pairscale_mode
//' @return a numeric vector that represents the normalization factors that were applied to each column in x. Note that x is updated by reference.
//' @export
// [[Rcpp::export]]
Rcpp::NumericVector pairscale_mean(arma::mat& x, Rcpp::Nullable<Rcpp::IntegerVector> clusters = R_NilValue, int min_value_count = 3, int niter_irls = 50, std::string na_mode = "check") {
  const int Ncols = x.n_cols;
  validate_minvaluecount(x.n_rows, min_value_count);
  bool check_na = validate_namode_mat(x, na_mode);
  arma::uvec clusters_cpp = validate_clusters(Ncols, clusters);
  validate_niterirls(niter_irls);

  // default, not used @ mean
  double tmean_trim = 0.2, density_adjust = 1, density_kernel_width_in_sd = 3, density_mode_frac_maxdens = 1, threshold_std = 3;
  int density_npoints = 200;
  std::string density_bandwidth_method = "nrd";
  Rcpp::NumericVector result = Rcpp::wrap(_pairscale_normalization(x, clusters_cpp, "mean", min_value_count, density_npoints, density_adjust, density_kernel_width_in_sd, density_bandwidth_method, density_mode_frac_maxdens, tmean_trim, threshold_std, niter_irls, check_na));
  result.attr("dim") = R_NilValue;
  return result;
}



//' @title Normalize matrix columns using their trimmed-mean differences
//' @description Pairwise normalization of columns in a matrix, using the trimmed-mean to define pairwise distances between columns
//' @name pairscale_trimmedmean
//' @inheritParams pairscale_mode
//' @param trim amount of trim to apply to both the lower- and upper-parts of a vector before computing the mean. 0 indicates no trim, 0.5 indicates 100% trim (i.e. 50% of data on both sides) so that value is out of bounds. Typically set to 0.1-0.3
//' @return a numeric vector that represents the normalization factors that were applied to each column in x. Note that x is updated by reference.
//' @export
// [[Rcpp::export]]
Rcpp::NumericVector pairscale_trimmedmean(arma::mat& x, Rcpp::Nullable<Rcpp::IntegerVector> clusters = R_NilValue, int min_value_count = 3, double trim = 0.2, int niter_irls = 50, std::string na_mode = "check") {
  const int Ncols = x.n_cols;
  validate_minvaluecount(x.n_rows, min_value_count);
  bool check_na = validate_namode_mat(x, na_mode);
  validate_trim(trim);
  arma::uvec clusters_cpp = validate_clusters(Ncols, clusters);
  validate_niterirls(niter_irls);

  // default, not used @ trimmed mean
  double density_adjust = 1, density_kernel_width_in_sd = 3, density_mode_frac_maxdens = 1, threshold_std = 3;
  int density_npoints = 200;
  std::string density_bandwidth_method = "nrd";
  Rcpp::NumericVector result = Rcpp::wrap(_pairscale_normalization(x, clusters_cpp, "trimmedmean", min_value_count, density_npoints, density_adjust, density_kernel_width_in_sd, density_bandwidth_method, density_mode_frac_maxdens, trim, threshold_std, niter_irls, check_na));
  result.attr("dim") = R_NilValue;
  return result;
}



//' @title Normalize matrix columns using their MAD-trimmed mean differences
//' @description Pairwise normalization of columns in a matrix, using the MAD-trimmed mean to define pairwise distances between columns
//' @name pairscale_madmean
//' @inheritParams pairscale_mode
//' @param threshold_std ratio of MAD a value has to be away from the median to be considered an outlier (and thus removed/ignored). Note that the MAD thresholds are inclusive, i.e. values at +/- threshold_std*MAD from median are included
//' @return a numeric vector that represents the normalization factors that were applied to each column in x. Note that x is updated by reference.
//' @export
// [[Rcpp::export]]
Rcpp::NumericVector pairscale_madmean(arma::mat& x, Rcpp::Nullable<Rcpp::IntegerVector> clusters = R_NilValue, int min_value_count = 3, double threshold_std = 3, int niter_irls = 50, std::string na_mode = "check") {
  const int Ncols = x.n_cols;
  validate_minvaluecount(x.n_rows, min_value_count);
  validate_madmeanthreshold(threshold_std);
  bool check_na = validate_namode_mat(x, na_mode);
  arma::uvec clusters_cpp = validate_clusters(Ncols, clusters);
  validate_niterirls(niter_irls);

  // default, not used @ madmean
  double tmean_trim = 0.2, density_adjust = 1, density_kernel_width_in_sd = 3, density_mode_frac_maxdens = 1;
  int density_npoints = 200;
  std::string density_bandwidth_method = "nrd";
  Rcpp::NumericVector result = Rcpp::wrap(_pairscale_normalization(x, clusters_cpp, "madmean", min_value_count, density_npoints, density_adjust, density_kernel_width_in_sd, density_bandwidth_method, density_mode_frac_maxdens, tmean_trim, threshold_std, niter_irls, check_na));
  result.attr("dim") = R_NilValue;
  return result;
}





//' @title graph Laplacian approach to finding normalization factors
//' @description find normalization factors for a given distance matrix computed with e.g. `pairdiff_median()`. For increased robustness, this function offers iterative reweighted improvement of the initial estimate.
//' @name solve_graph_laplacian
//' @param M skew-symmetric input matrix, generated with e.g. `pairdiff_median()`
//' @param niter_irls refine the initial estimate using N additional iterative reweighted least squares loops for robust graph laplacian
//' @examples
//' # toy example
//' x = cbind(
//'   c(1,2,3,4),
//'   c(2,3,4,9),
//'   c(1,2,4,5),
//'   c(1,0,1,0)
//' )
//' # compute pairwide median difference between all columns
//' M = pairscale::pairdiff_median(x)
//' # solve matrix M to find scaling factors, without and with reweighting
//' s1 = pairscale::solve_graph_laplacian(M, niter_irls = 0)
//' s2 = pairscale::solve_graph_laplacian(M, niter_irls = 10)
//' # rescaled matrices; only the robust variant correctly aligns columns 1 and 2
//' t(t(x) - s1)
//' t(t(x) - s2)
//' @return a numeric vector of length `ncol(M)` that contains scaling factors for `M`
//' @export
// [[Rcpp::export]]
Rcpp::NumericVector solve_graph_laplacian(arma::mat M, int niter_irls = 1) {
  Rcpp::NumericVector result = Rcpp::wrap(_solve_graph_laplacian(M, niter_irls));
  result.attr("dim") = R_NilValue;
  return result;
}

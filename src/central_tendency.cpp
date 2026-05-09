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

#include "central_tendency.h"



// this function updates x by reference (possible sorting of nonfinite values)
double _mean_vector(arma::vec& x, int min_value_count, bool check_na) {
  int n_valid = x.n_elem;
  if(check_na) {
    n_valid = move_nonfinite_to_end(x);
  }

  if (n_valid < min_value_count) {
    return arma::datum::nan;
  }
  if (n_valid == 1) {
    return x[0];
  }

  // v.head(n) returns a lightweight, zero copy/allocation/copy/overhead view (subview_col) into the original vector's memory
  return arma::mean(x.head(n_valid));
}



// this function updates x by reference (possible sorting of nonfinite values)
double _median_vector(arma::vec& x, int min_value_count, bool check_na) {
  int n_valid = x.n_elem;
  if(check_na) {
    n_valid = move_nonfinite_to_end(x);
  }

  if (n_valid < min_value_count) {
    return arma::datum::nan;
  }
  if (n_valid == 1) {
    return x[0];
  }
  if (n_valid == 2) {
    return (x[0] + x[1]) / 2.0;
  }

  // median over the leading N finite values, faster than arma::median() but partial sorts in place !
  return median_in_partialvector(x, n_valid);
}



// this function updates x by reference (possible sorting of nonfinite values)
double _trimmedmean_vector(arma::vec& x, int min_value_count, bool check_na, double trim) {
  if(trim < 0 || trim >= 0.5) {
    throw std::invalid_argument("'trim' must be a value between 0 and 0.5 (i.e. at 0.5 one would trim 50% from both sides, ergo there are no data to be used)");
  }

  int n_valid = x.n_elem;
  if(check_na) {
    n_valid = move_nonfinite_to_end(x);
  }

  if (n_valid < min_value_count) {
    return arma::datum::nan;
  }
  if (n_valid == 1) {
    return x[0];
  }
  if (n_valid == 2) {
    return (x[0] + x[1]) / 2.0;
  }

  // number of values to cut on each side. If 0, return mean. If too many, return NA
  int cut = static_cast<int>(std::floor(trim * n_valid));
  if(cut == 0) {
    return arma::mean(x.head(n_valid));
  }
  if(cut + cut >= n_valid) {
    return arma::datum::nan;
  }

  // zero-copy implementation using pointers to armadillo buffer
  double* begin = x.memptr();
  double* end   = begin + n_valid;  // importantly, only consider first N elements
  double* lo    = begin + cut;
  double* hi    = begin + n_valid - cut;

  // Partial sort: [begin, lo) ≤ *lo ≤ [lo+1, end)
  std::nth_element(begin, lo, end);
  // Partial sort again on subset past 'lo'; [lo+1, hi) ≤ *hi ≤ [hi+1, end)
  std::nth_element(lo + 1, hi, end);

  // param3; strict = true implies no size check
  // param4; copy_aux_mem = false implies no allocation
  return arma::mean(arma::vec(lo, hi - lo, false, true));
}



arma::vec _trimmedmean_matrixrows(const arma::mat& x, int min_value_count, bool check_na, double trim) {
  if(trim < 0 || trim >= 0.5) {
    throw std::invalid_argument("'trim' must be a value between 0 and 0.5 (i.e. at 0.5 one would trim 50% from both sides, ergo there are no data to be used)");
  }
  const int Ncol = x.n_cols;
  const int Nrow = x.n_rows;
  arma::vec buffer = arma::vec(Ncol, arma::fill::none);
  arma::vec result = arma::vec(Nrow, arma::fill::none);

  // optimize by defering to regular mean computation function if requested trim count is zero (e.g. 0.2 trim for 3 columns)
  bool do_tmeans = Ncol * trim >= 1;
  if(do_tmeans) {
    for(int i = 0; i < Nrow; i++) {
      // copy arma::rowvec data into an arma::vec. We need to copy because downstream code will modify the buffer vector by reference
      buffer = arma::conv_to<arma::vec>::from(x.row(i).t());
      result(i) = _trimmedmean_vector(buffer, min_value_count, check_na, trim);
    }
  } else {
    for(int i = 0; i < Nrow; i++) {
      buffer = arma::conv_to<arma::vec>::from(x.row(i).t());
      result(i) = _mean_vector(buffer, min_value_count, check_na);
    }
  }

  return result;
}



// this function updates x by reference
// the MAD thresholds are inclusive, i.e. values at threshold_std*MAD from median are included
double _madmean_vector(arma::vec& x, arma::vec& buffer, int min_value_count, bool check_na, double threshold_std) {
  int n_valid = x.n_elem;
  if(check_na) {
    n_valid = move_nonfinite_to_end(x);
  }

  if (n_valid < min_value_count) {
    return arma::datum::nan;
  }
  if (n_valid == 1) {
    return x[0];
  }
  if (n_valid == 2) {
    return (x[0] + x[1]) / 2.0;
  }

  // median over the leading N finite values, faster than arma::median() but partial sorts in place !
  const double median = median_in_partialvector(x, n_valid);
  if (n_valid == 3) {
    return median;
  }

  // buffer = absolute deviation from median
  for (int i = 0; i < n_valid; ++i) {
    buffer[i] = std::abs(x[i] - median);
  }

  // MAD = find the median of absolute deviations
  const double mad = median_in_partialvector(buffer, n_valid); // partial sorts in place !

  // rescale MAD to approx equivalent on normally distributed data
  const double scaled_mad = mad * 1.4826;
  // Rcpp::Rcout << "n_valid: " << n_valid << " median: " << median << " mad: " << mad << " scaled_mad: " << scaled_mad << " threshold_std: " << threshold_std << "\n";

  // determine filter thresholds at median +/ MAD*n_std, then compute mean of remaining values
  const double lower_threshold = median - (threshold_std * scaled_mad);
  const double upper_threshold = median + (threshold_std * scaled_mad);
  // edge case: min and max are the same or have tiny difference, so return their mean
  if(upper_threshold - lower_threshold < 0.00001) {
    return median;
  }

  long double sum = 0.0;
  int n_sum = 0;
  for (int i = 0; i < n_valid; ++i) {
    double val = x[i];
    if(val >= lower_threshold && val <= upper_threshold) {
      sum += val;
      n_sum++;
    }
  }
  return static_cast<double>(sum / n_sum);
}



double mode_from_binned_data(const arma::vec& binned, double minval, double step, int n_datapoint_for_kde, double var_estimate, double kernel_width_in_sd, double adjust, double mode_frac_maxdens) {

  // estimate bandwidth parameter
  // analogous to bw.nrd()
  double bandwidth = 1.06 * var_estimate * std::pow(n_datapoint_for_kde, -0.2) * adjust;

  arma::vec kernel = gauss_kernel(bandwidth, n_datapoint_for_kde, kernel_width_in_sd, step);

  // Rcpp::Rcout << "n_valid: " << n_valid << " minval: " << minval << " maxval: " << maxval << " var_estimate: " << var_estimate << " index_start: " << index_start << " index_endincl: " << index_endincl << "\n";

  if(mode_frac_maxdens >= 1) {
    double mode_idx = fast_kde_mode(binned, kernel);
    return minval + mode_idx * step + step/2.0; // offset by half a bin = return bin center
  } else {
    arma::vec binned_y(binned.n_elem, arma::fill::none);
    double mode_idx = fast_kde_modeareamean(binned, binned_y, kernel, mode_frac_maxdens);
    return minval + mode_idx * step + step/2.0; // offset by half a bin = return bin center
  }
}



// this function updates x by reference
// qtl typically is 0.2
bool outlier_scan_progressive(arma::vec& x, int n_valid, double qtl, int& index_lower, int& index_upper) {
  index_lower = 0;
  index_upper = 0;

  int cut_low = static_cast<int>(std::floor(qtl * n_valid));
  int cut_up = static_cast<int>(std::floor((1-qtl) * n_valid));
  if (cut_up == 0) {
    return false;
  }

  double* begin = x.memptr();
  // Partial sort: [begin, lo) ≤ *lo ≤ [lo+1, end)
  std::nth_element(begin, begin + cut_low, begin + n_valid);
  // Partial sort again on subset past 'lo'; [lo+1, hi) ≤ *hi ≤ [hi+1, end)
  std::nth_element(begin + cut_low + 1, begin + cut_up, begin + n_valid);

  index_lower = cut_low;
  index_upper = cut_up;
  double qtl_diff = x[cut_up] - x[cut_low];
  if (qtl_diff < 0.00001) {
    return false;
  }

  // thresholding based on distance to nearest threshold; cannot be more than X times the distance between outer quantiles
  double localdiff_threshold = qtl_diff * 0.5;

  // above code performed partial sorting of x
  // now move outwards from lower quantile and stop when diff to reference value is too large
  for(int i = cut_low - 1; i >= 0; --i) {
    // x[i] is left of x[cut_low], since we partially sorted x[i] is lower or equal
    // once the difference between x[i] and x[cut_low] is too large, stop/break
    if(localdiff_threshold < x[cut_low] - x[i]) {
      break;
    }
    index_lower = i;
  }

  // analogous, but from upper value moving to end of the vector
  for(int i = cut_up + 1; i < n_valid; ++i) {
    if(localdiff_threshold < x[i] - x[cut_up]) {
      break;
    }
    index_upper = i;
  }

  // Rcpp::Rcout << "outlier_scan_progressive() cut_low: " << cut_low << " cut_up: " << cut_up << " x[cut_low]: " << x[cut_low] << " x[cut_up]: " << x[cut_up] << " index_lower: " << index_lower << " index_upper: " << index_upper << "\n";

  return true;
}


// Sorting is a computationally expensive operation, so the slowest approach is where we compute IQR since we have to
// sort most of the data to find the 0.25 and 0.75 quantiles. Higher n_bins also has large impact on computation time
// this function updates x by reference
double _mode_vector(arma::vec& x, int min_value_count, bool check_na, double adjust, double kernel_width_in_sd, int n_bins, std::string bandwidth_method, double mode_frac_maxdens) {
  double* x_ptr_begin = x.memptr();
  int n_valid = x.n_elem;
  int index_lower, index_upper;
  double minval = 0.0, maxval = 0.0;

  if(check_na) {
    n_valid = move_nonfinite_to_end(x);
  }

  if (n_valid < min_value_count) {
    return arma::datum::nan;
  }
  if (n_valid == 1) {
    return x[0];
  }
  if (n_valid == 2) {
    return (x[0] + x[1]) / 2.0;
  }
  if (n_valid <= 5) {
    // median over the leading N finite values, faster than arma::median() but partial sorts in place !
    return median_in_partialvector(x, n_valid);
  }



  //// outlier detection using approximate MAD; robust but slow for large datasets because we have to find median twice (sorting is expensive)
  //// concept/idea, code untested
  //
  // int mid = n_valid / 2; // integer/2; floor() is implicit
  // double* x_begin = x.memptr();
  // std::nth_element(x_begin, x_begin + mid, x_begin + n_valid);
  // double median_est = x[mid];
  // //// we don't care so much about accuracy here, so skip below edge case for even N and use nth element median instead
  // //if (n_valid % 2 == 0) {
  // //  double max_left = *std::max_element(begin, begin + mid);
  // //  median = (median + max_left) / 2.0;
  // //}
  // //// analogous for MAD; also use first index and skip even N edge case
  // // TODO this could go into a buffer we recycle between function calls
  // arma::vec dist = x.head(n_valid) - median_est;
  // double* dist_begin = x.memptr();
  // std::nth_element(dist_begin, dist_begin + mid, dist_begin + n_valid);
  // double mad_est = dist[mid] * 1.4826;
  //
  // const double lower_threshold = median_est - mad_est * 4;
  // const double upper_threshold = median_est + mad_est * 4;



  //// efficient outlier detection using quantiles
  // trim outer 5% when N > 1000, outer 10% when N > 100,
  // test for outlier magnitude when N <= 100 to retain maximum number of datapoints
  // latter case is needed because of a loss of precision when N=10 is high when naively trimming
  // and for large datasets we don't want to trim 10% for performance reasons (sorting is expensive)
  // this is fast because we minimize the amount of sorting we have to do.

  if(bandwidth_method == "nrd_fastest") {
    // fastest; use all data points, for bandwidth don't use IQR. Sensitive to outliers !
    // min/max are only computed upstream if there were NA's
    arma_vec_minmax(x, n_valid, minval, maxval);
    index_lower = 0;
    index_upper = n_valid - 1; // last index that holds valid data

  } else {

    //// bandwidth_method: nrd_fast or nrd

    if(n_valid < 50) {
      bool success = outlier_scan_progressive(x, n_valid, 0.2, index_lower, index_upper); // updated x by reference
      if(!success) {
        return (x[index_lower] + x[index_upper]) / 2.0;
      }
    } else {

      double qtl = 0.025;
      if (n_valid < 500) qtl = 0.05;
      else if (n_valid < 600) qtl = 0.045;
      else if (n_valid < 700) qtl = 0.04;
      else if (n_valid < 800) qtl = 0.035;
      else if (n_valid < 900) qtl = 0.03;

      index_lower = static_cast<int>(std::floor(qtl * n_valid));
      index_upper = static_cast<int>(std::floor((1-qtl) * n_valid));
      // Partial sort: [begin, lo) ≤ *lo ≤ [lo+1, end)
      std::nth_element(x_ptr_begin, x_ptr_begin + index_lower, x_ptr_begin + n_valid);
      // Partial sort again on subset past 'lo'; [lo+1, hi) ≤ *hi ≤ [hi+1, end)
      std::nth_element(x_ptr_begin + index_lower + 1, x_ptr_begin + index_upper, x_ptr_begin + n_valid);
    }

    minval = x[index_lower];
    maxval = x[index_upper];
  }

  // Rcpp::Rcout << "bandwidth_method: " << bandwidth_method << " index_lower: " << index_lower << " index_upper: " << index_upper << " minval: " << minval << " maxval: " << maxval << "\n";


  // edge case: min and max are the same or have tiny difference, so return their mean
  double range = maxval - minval;
  if (range < 0.00001) {
    return (minval + maxval) / 2.0;
  }

  // for tiny vectors, add some padding as a function of total value range.
  // not an issue for distributions that are approx normal and have 25+ data points
  //
  // this is a poor man's fix for lack of out-of-bounds zero value weighting
  // a better solution. is probably to weigh zeros for out of bounds values while applying the KDE
  if(n_valid < 50) {
    minval = 0.7 * minval;
    maxval = 1.4 * maxval;
    range = maxval - minval;
  }

  // standard deviation over selected data range
  // .subvec() should be zero cost/overhead subsetting
  double sd = arma::stddev(x.subvec(index_lower, index_upper));

  // redefine sd as min(sd, iqr). Adds robustness, but finding quantiles is expensive/slow
  if(bandwidth_method == "nrd") {
    // analogous to above nth_element() code
    int idx25 = static_cast<int>(std::floor(0.25 * n_valid));
    int idx75 = static_cast<int>(std::floor(0.75 * n_valid));
    std::nth_element(x_ptr_begin + index_lower + 1, x_ptr_begin + idx25, x_ptr_begin + n_valid);
    std::nth_element(x_ptr_begin + idx25 + 1, x_ptr_begin + idx75, x_ptr_begin + n_valid);
    double q25 = x[idx25];
    double q75 = x[idx75];
    double iqr_div = (q75 - q25) / 1.34;
    sd = std::min(sd, iqr_div);
  }


  // finally, linear binning and KDE
  double step;
  arma::vec binned = linear_binning(x, index_lower, index_upper, minval, range, n_bins, step);
  // Rcpp::Rcout << binned << "\n";
  // Rcpp::Rcout << "range: " << range << " n_bins: " << n_bins << " step: " << step << " index_upper - index_lower + 1 :" << (index_upper - index_lower + 1) << "\n";
  return mode_from_binned_data(binned, minval, step, index_upper - index_lower + 1, sd, kernel_width_in_sd, adjust, mode_frac_maxdens);
}



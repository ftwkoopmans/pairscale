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

#include "pairwise_distance.h"
#include "central_tendency.h"

#ifdef _OPENMP
#include <omp.h>
#endif






// cols must be zero-indexed, C++ style, and validated upstream
arma::mat _pairwise_distance_mean(const arma::mat& x, const arma::uvec& cols_cpp, int min_value_count, bool check_na) {
  const int N_input = x.n_cols;
  const int M = x.n_rows;
  const bool is_valid = N_input > 1 && M >= min_value_count;
  if(!is_valid) {
    return arma::mat(std::min(N_input, (int)cols_cpp.n_elem), std::min(N_input, (int)cols_cpp.n_elem), arma::fill::nan);
  }
  const int N = cols_cpp.n_elem;

  arma::mat dist = arma::mat(N, N, arma::fill::none);
  dist(0,0) = 0.0; // set zeros to the diagonal, coord (0,0) is never reached in below loop

  // create parallel region (if available)
#ifdef _OPENMP
  // 1. Create the parallel region. Threads are spawned here.
#pragma omp parallel if(N > 10)
#endif
{
  // variables declared right here are thread-local;
  // - if openmp is available, each thread instantiates these variables (once)
  // - otherwise without openmp, these variable are created exactly 1 time
  arma::vec ij_diff(M, arma::fill::none);

  // openmp specification for scheduling threads (if available)
#ifdef _OPENMP
#pragma omp for schedule(dynamic, 4)
#endif

  for(int i = 1; i < N; i++) { // start at 1 / skip first column
    dist(i,i) = 0.0; // set zeros to the diagonal
    for(int j = 0; j < i; j++) { // iterate only unique pairs
      ij_diff = x.col(cols_cpp(i)) - x.col(cols_cpp(j));
      double ij_dist = _mean_vector(ij_diff, min_value_count, check_na);
      dist(i,j) = ij_dist;
      dist(j,i) = -1 * ij_dist;
    }
  }

}

return dist;
}





// cols must be zero-indexed, C++ style, and validated upstream
arma::mat _pairwise_distance_median(const arma::mat& x, const arma::uvec& cols_cpp, int min_value_count, bool check_na) {
  const int N_input = x.n_cols;
  const int M = x.n_rows;
  const bool is_valid = N_input > 1 && M >= min_value_count;
  if(!is_valid) {
    return arma::mat(std::min(N_input, (int)cols_cpp.n_elem), std::min(N_input, (int)cols_cpp.n_elem), arma::fill::nan);
  }
  const int N = cols_cpp.n_elem;

  arma::mat dist = arma::mat(N, N, arma::fill::none);
  dist(0,0) = 0.0; // set zeros to the diagonal, coord (0,0) is never reached in below loop

  // openmp code comments; see first function that uses openmp in this file
#ifdef _OPENMP
#pragma omp parallel if(N > 10)
#endif
{
  arma::vec ij_diff(M, arma::fill::none);

#ifdef _OPENMP
#pragma omp for schedule(dynamic, 4)
#endif

  for(int i = 1; i < N; i++) { // start at 1 / skip first column
    dist(i,i) = 0.0; // zero the diagonal
    for(int j = 0; j < i; j++) { // iterate only unique pairs
      ij_diff = x.col(cols_cpp(i)) - x.col(cols_cpp(j));
      double ij_dist = _median_vector(ij_diff, min_value_count, check_na);
      dist(i,j) = ij_dist;
      dist(j,i) = -1 * ij_dist;
    }
  }

}

return dist;
}





//current implementation always checks for NA, i.e. we didn't implement a faster path (yet)
// cols must be zero-indexed, C++ style, and validated upstream
arma::mat _pairwise_distance_trimmedmean(const arma::mat& x, const arma::uvec& cols_cpp, int min_value_count, double trim, bool check_na) {
  const int N_input = x.n_cols;
  const int M = x.n_rows;
  const bool is_valid = N_input > 1 && M >= min_value_count;
  if(!is_valid) {
    return arma::mat(std::min(N_input, (int)cols_cpp.n_elem), std::min(N_input, (int)cols_cpp.n_elem), arma::fill::nan);
  }
  const int N = cols_cpp.n_elem;

  arma::mat dist = arma::mat(N, N, arma::fill::none); // fill none is a minor speedup
  dist(0,0) = 0.0; // set zeros to the diagonal, coord (0,0) is never reached in below loop

  // openmp code comments; see first function that uses openmp in this file
#ifdef _OPENMP
#pragma omp parallel if(N > 10)
#endif
{
  arma::vec ij_diff(M, arma::fill::none);

#ifdef _OPENMP
#pragma omp for schedule(dynamic, 4)
#endif

  for(int i = 1; i < N; i++) { // start at 1 / skip first column
    dist(i,i) = 0.0; // set zeros to the diagonal
    for(int j = 0; j < i; j++) { // iterate only unique pairs
      ij_diff = x.col(cols_cpp(i)) - x.col(cols_cpp(j));
      double ij_dist = _trimmedmean_vector(ij_diff, min_value_count, check_na, trim);
      dist(i,j) = ij_dist;
      dist(j,i) = -1 * ij_dist;
    }
  }

}

return dist;
}





//current implementation always checks for NA, i.e. we didn't implement a faster path (yet)
// cols must be zero-indexed, C++ style, and validated upstream
arma::mat _pairwise_distance_madmean(const arma::mat& x, const arma::uvec& cols_cpp, int min_value_count, double threshold_std, bool check_na) {
  const int N_input = x.n_cols;
  const int M = x.n_rows;
  const bool is_valid = N_input > 1 && M >= min_value_count;
  if(!is_valid) {
    return arma::mat(std::min(N_input, (int)cols_cpp.n_elem), std::min(N_input, (int)cols_cpp.n_elem), arma::fill::nan);
  }
  const int N = cols_cpp.n_elem;

  arma::mat dist = arma::mat(N, N, arma::fill::none); // fill none is a minor speedup
  dist(0,0) = 0.0; // set zeros to the diagonal, coord (0,0) is never reached in below loop

  // openmp code comments; see first function that uses openmp in this file
#ifdef _OPENMP
#pragma omp parallel if(N > 10)
#endif
{
  arma::vec ij_diff(M, arma::fill::none);
  arma::vec buffer(M, arma::fill::none);

#ifdef _OPENMP
#pragma omp for schedule(dynamic, 4)
#endif

  for(int i = 1; i < N; i++) { // start at 1 / skip first column
    dist(i,i) = 0.0; // set zeros to the diagonal
    for(int j = 0; j < i; j++) { // iterate only unique pairs
      ij_diff = x.col(cols_cpp(i)) - x.col(cols_cpp(j));
      double ij_dist = _madmean_vector(ij_diff, buffer, min_value_count, check_na, threshold_std);
      dist(i,j) = ij_dist;
      dist(j,i) = -1 * ij_dist;
    }
  }

}

return dist;
}





// cols must be zero-indexed, C++ style, and validated upstream
arma::mat _pairwise_distance_mode(const arma::mat& x, const arma::uvec& cols_cpp, int min_value_count, int n_bins, double adjust, double kernel_width_in_sd, std::string bandwidth_method, double mode_frac_maxdens, bool check_na) {
  const int N_input = x.n_cols;
  const int M = x.n_rows;
  const bool is_valid = N_input > 1 && M >= min_value_count;
  if(!is_valid) {
    return arma::mat(std::min(N_input, (int)cols_cpp.n_elem), std::min(N_input, (int)cols_cpp.n_elem), arma::fill::nan);
  }
  const int N = cols_cpp.n_elem;

  arma::mat dist = arma::mat(N, N, arma::fill::none);
  dist(0,0) = 0.0; // set zeros to the diagonal, coord (0,0) is never reached in below loop

  // openmp code comments; see first function that uses openmp in this file
#ifdef _OPENMP
#pragma omp parallel if(N > 10)
#endif
{
  arma::vec ij_diff(M, arma::fill::none);

#ifdef _OPENMP
#pragma omp for schedule(dynamic, 4)
#endif


  for(int i = 1; i < N; i++) { // start at 1 / skip first column
    dist(i,i) = 0.0; // set zeros to the diagonal
    for(int j = 0; j < i; j++) { // iterate only unique pairs
      ij_diff = x.col(cols_cpp(i)) - x.col(cols_cpp(j));
      double ij_dist = _mode_vector(ij_diff, min_value_count, check_na, adjust, kernel_width_in_sd, n_bins, bandwidth_method, mode_frac_maxdens);
      dist(i,j) = ij_dist;
      dist(j,i) = -1 * ij_dist;
    }
  }

}

return dist;
}


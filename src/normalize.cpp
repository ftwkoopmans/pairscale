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

#include "normalize.h"
#include "pairwise_distance.h"
#include "central_tendency.h"



arma::vec laplacian_solver(arma::mat L, arma::vec v) {
  int n = L.n_rows;

  // 1. Drop the last row and column of L, and last element of v
  arma::mat L_sub = L.submat(0, 0, n - 2, n - 2);
  arma::vec v_sub = v.subvec(0, n - 2);

  // 2. Solve using standard, ultra-fast linear solver
  // arma::solve uses LAPACK algorithms internally (much faster than pinv)
  arma::vec s_sub = arma::solve(L_sub, v_sub);

  // 3. Re-append the dropped element (as 0)
  arma::vec s(n, arma::fill::zeros);
  s.subvec(0, n - 2) = s_sub;

  // 4. Mean-center to match exactly what MASS::ginv does
  s = s - arma::mean(s);

  return s;
}



//' @title graph Laplacian approach to finding normalization factors
//' @description find normalization factors for a given distance matrix computed with e.g. `pairdiff_median()`. For increased robustness, this function offers iterative reweighted improvement of the initial estimate.
//' @name solve_graph_laplacian
//' @param M skew-symmetric input matrix, generated with e.g. `pairdiff_median()`
//' @param niter_irls refine the initial estimate using N additional iterative reweighted least squares loops for robust graph laplacian
//' @examples
//' \dontrun{
//' # toy example
//' x = cbind(
//'   c(1,2,3,4),
//'   c(2,3,4,9),
//'   c(1,2,4,5),
//'   c(1,0,1,0)
//' )
//' # compute pairwide median difference between all columns
//' M = pairdiff_median(x)
//' # solve matrix M to find scaling factors, without and with reweighting
//' s1 = pairscale::solve_graph_laplacian(M, niter_irls = 0)
//' s2 = pairscale::solve_graph_laplacian(M, niter_irls = 10)
//' # rescaled matrices; only the robust variant correctly aligns columns 1 and 2
//' t(t(x) - s1[,1])
//' t(t(x) - s2[,1])
//' }
//' @export
// [[Rcpp::export]]
arma::vec solve_graph_laplacian(arma::mat M, int niter_irls = 1) {
  const int n = M.n_rows;

  // Adjacency matrix A: 1 if M(i,j) is finite (not NA/NaN), 0 otherwise
  // Matrix M_zero is a copy of M, with NAs replaced by zero for safe matrix math
  arma::mat A(n, n, arma::fill::none);
  arma::mat M_zero(n, n, arma::fill::none);

  // iterate all values in M; prepare matrices A and M_zero for efficient matrix math
  for (int i = 0; i < n; i++) {
    for (int j = 0; j < n; j++) {
      if (i != j && std::isfinite(M(i, j))) {
        A(i, j)      = 1.0;
        M_zero(i, j) = M(i, j);
      } else {
        // zero the diagonal and NA values
        A(i, j)      = 0.0;
        M_zero(i, j) = 0.0;
      }
    }
  }

  // Degree matrix D
  // R code: D = diag(rowSums(A))
  arma::vec deg = arma::sum(A, 1); // column vector of row sums
  arma::mat D   = arma::diagmat(deg);
  // Graph Laplacian L = D - A
  arma::mat L = D - A;
  // Target vector v = rowSums(M_zero)
  arma::vec v = arma::sum(M_zero, 1);
  // Solve
  arma::vec s_l2 = laplacian_solver(L, v);


  // if robust reweighting was not requested, immediately exit
  if (niter_irls < 1) {
    return s_l2;
  }


  //// iterative reweighted least squares loop for robust graph laplacian
  arma::vec s_irls = s_l2;      // warm-start from the L2 solution
  const double epsilon  = 1e-4; // Small constant to prevent division by zero

  arma::mat S_diff(n, n, arma::fill::none);
  arma::mat residuals(n, n, arma::fill::none);
  arma::mat W(n, n, arma::fill::none);
  arma::mat L_w(n, n, arma::fill::none);
  arma::vec deg_w(n, arma::fill::none);
  arma::vec v_w(n, arma::fill::none);

  for (int iter = 0; iter < niter_irls; iter++) {

    // Step 1: Calculate expected measurements based on current estimates
    // R code: S_diff <- outer(s_irls, s_irls, "-")
    // Allocation-free outer difference, written directly into pre-allocated S_diff
    for (int i = 0; i < n; i++) {
      for (int j = 0; j < n; j++) {
        S_diff(i, j) = s_irls(i) - s_irls(j);
      }
    }

    // Step 2: absolute residuals
    residuals = arma::abs(M_zero - S_diff);

    // Step 3: weights masked by A
    // values are zeroed where A == 0 ; this is important when dealing with NA values !
    // i.e. M_zero will have 0 for NA, updating with S_diff witll make these nonzero
    // thus impacting our metrics. Multiplying by A prevents this
    W = (1.0 / (residuals + epsilon)) % A;   // % = element-wise multiply

    // Step 4: weighted graph Laplacian
    // R code: L_w <- diag(rowSums(W)) - W
    deg_w = arma::sum(W, 1);
    L_w   = -W;                          // copy W negated
    L_w.diag() += deg_w;                 // add degree on diagonal in-place

    // Step 5: weighted target vector
    v_w = arma::sum(W % M_zero, 1);

    // Step 6: solve
    s_irls = laplacian_solver(L_w, v_w);
  }

  return s_irls;
}



// assumes upstream code performed input validation !
// clusters; positive, finite integers
// check_na = check_na && (x.has_nan() || x.has_inf()); // this check should be performed upstream
arma::vec _pairscale_normalization(arma::mat& x, const arma::uvec& clusters, std::string centeral_tendency_measure, int min_value_count, int density_npoints, double density_adjust, double density_kernel_width_in_sd, std::string bandwidth_method, double mode_frac_maxdens, double tmean_trim, double threshold_std, int niter_irls, bool check_na) {

  arma::uvec uclust = arma::unique(clusters);

  // init vars
  int nclust = uclust.n_elem;
  int i_ncol = 0;
  arma::mat profiles;
  arma::uvec cols; // a vector of zero-indexed column indices @ matrix x
  arma::mat distmat;
  arma::vec factors;
  arma::vec all_factors = arma::vec(x.n_cols, arma::fill::zeros);

  // if there are multiple clusters, create a matrix with mean values per row per cluster for between-cluster scaling
  if(nclust > 1) {
    profiles = arma::mat(x.n_rows, nclust, arma::fill::none);
  }

  // iterate clusters
  for(int i = 0; i < nclust; i++) {
    // use Armadillo find() to return indices within clusters that match current cluster ID
    cols = arma::find(clusters == uclust(i));
    i_ncol = cols.n_elem;

    // Rcpp::Rcout << "cols: " << cols << "\n";

    // current cluster ID has only 1 column: store profile as the data as-is and done
    if(i_ncol == 1) {
      profiles.col(i) = x.col(cols(0));
    } else {

      //// R equivalent of this section;
      // pm = pairwise_distance_median(x, min_value_count = 3, check_na = TRUE)
      // s = graph_laplacian(pm, niter_reweighted = 5)
      // t(t(x) - s)

      if(centeral_tendency_measure == "mean") {
        distmat = _pairwise_distance_mean(x, cols, min_value_count, check_na);
      } else {
        if(centeral_tendency_measure == "median") {
          distmat = _pairwise_distance_median(x, cols, min_value_count, check_na);
        } else {
          if(centeral_tendency_measure == "madmean") {
            distmat = _pairwise_distance_madmean(x, cols, min_value_count, threshold_std, check_na);
          } else {
            if(centeral_tendency_measure == "trimmedmean") {
              distmat = _pairwise_distance_trimmedmean(x, cols, min_value_count, tmean_trim, check_na);
            } else {
              if(centeral_tendency_measure == "mode") {
                distmat = _pairwise_distance_mode(x, cols, min_value_count, density_npoints, density_adjust, density_kernel_width_in_sd, bandwidth_method, mode_frac_maxdens, check_na);
                // Rcpp::Rcout << "min_value_count: " << min_value_count << " density_npoints: " << density_npoints << " density_adjust: " << density_adjust << " density_kernel_width_in_sd: " << density_kernel_width_in_sd << " bandwidth_method: " << bandwidth_method << "\n";
                // Rcpp::Rcout << distmat << "\n";
              } else {
                throw std::invalid_argument("'centeral_tendency_measure' unknown option/metric is provided");
              }
            }
          }
        }
      }

      //distmat = _pairwise_distance_robustmode(x, cols, min_value_count, density_npoints, density_adjust, density_kernel_width_in_sd, check_na);
      // compute scaling factors for each column
      factors = solve_graph_laplacian(distmat, niter_irls);
      // store factors in overall result vector
      all_factors(cols) = factors;
      // Rcpp::Rcout << "factors: " << factors << "\n";


      // rescale subset of cols in x
      for(int j = 0; j < i_ncol; j++) {
        x.col(cols(j)) -= factors(j);
      }

      if(nclust > 1) {
        const int min_value_count = std::min(3, (int)(cols.n_elem / 2));
        // Rcpp::Rcout << "min_value_count: " << min_value_count << " cols.n_elem: " << cols.n_elem << "\n";
        // compute profile, but only if there are at least 2 clusters (i.e. we need these for between-cluster scaling)
        profiles.col(i) = _trimmedmean_matrixrows(x.cols(cols), min_value_count, check_na, 0.2);
      }
    }
  }

  // normalize between clusters
  if(nclust > 1) {
    arma::uvec cols_profiles(nclust, arma::fill::none);
    for(int i = 0; i < nclust; i++) {
      cols_profiles(i) = i;
    }

    if(centeral_tendency_measure == "mean") {
      distmat = _pairwise_distance_mean(profiles, cols_profiles, min_value_count, check_na);
    } else {
      if(centeral_tendency_measure == "median") {
        distmat = _pairwise_distance_median(profiles, cols_profiles, min_value_count, check_na);
      } else {
        if(centeral_tendency_measure == "madmean") {
          double threshold_std = 1;
          distmat = _pairwise_distance_madmean(profiles, cols_profiles, min_value_count, threshold_std, check_na);
        } else {
          if(centeral_tendency_measure == "trimmedmean") {
            distmat = _pairwise_distance_trimmedmean(profiles, cols_profiles, min_value_count, tmean_trim, check_na);
          } else {
            if(centeral_tendency_measure == "mode") {
              distmat = _pairwise_distance_mode(profiles, cols_profiles, min_value_count, density_npoints, density_adjust, density_kernel_width_in_sd, bandwidth_method, mode_frac_maxdens, check_na);
            } else {
              throw std::invalid_argument("'centeral_tendency_measure' unknown option/metric is provided");
            }
          }
        }
      }
    }

    factors = solve_graph_laplacian(distmat, niter_irls);

    // rescale matrix; iterate clusters and subtract factor to all columns
    for(int i = 0; i < nclust; i++) {
      cols = arma::find(clusters == uclust(i));
      x.cols(cols) -= factors(i);

      //update overall result vector; increase as opposed to subtract @ data matrix update
      all_factors(cols) += factors(i);
    }
  }

  return all_factors;
}



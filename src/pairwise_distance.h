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

#ifndef PAIRWISEDISTANCE_H
#define PAIRWISEDISTANCE_H

#include <RcppArmadillo.h>

arma::mat _pairwise_distance_mean(const arma::mat& x, const arma::uvec& cols, int min_value_count, bool check_na);
arma::mat _pairwise_distance_median(const arma::mat& x, const arma::uvec& cols, int min_value_count, bool check_na);
arma::mat _pairwise_distance_trimmedmean(const arma::mat& x, const arma::uvec& cols, int min_value_count, double trim, bool check_na);
arma::mat _pairwise_distance_madmean(const arma::mat& x, const arma::uvec& cols, int min_value_count, double threshold_std, bool check_na);
arma::mat _pairwise_distance_mode(const arma::mat& x, const arma::uvec& cols_cpp, int min_value_count, int n_bins, double adjust, double kernel_width_in_sd, std::string bandwidth_method, double mode_frac_maxdens, bool check_na);
//arma::mat _pairwise_distance_thresholdedmode(const arma::mat& x, const arma::uvec& cols, int min_value_count, int density_npoints, double density_adjust, double density_kernel_width_in_sd, double threshold_lower, double threshold_upper);

#endif

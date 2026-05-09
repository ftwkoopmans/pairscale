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

#ifndef CENTRALTENDENCY_H
#define CENTRALTENDENCY_H

#include <RcppArmadillo.h>
#include "util.h"
#include "kde.h"

double _mean_vector(arma::vec& x, int min_value_count, bool check_na);
double _median_vector(arma::vec& x, int min_value_count, bool check_na);
double _trimmedmean_vector(arma::vec& x, int min_value_count, bool check_na, double trim);
arma::vec _trimmedmean_matrixrows(const arma::mat& x, int min_value_count, bool check_na, double trim);
double _madmean_vector(arma::vec& x, arma::vec& buffer, int min_value_count, bool check_na, double threshold_std);
double _mode_vector(arma::vec& x, int min_value_count, bool check_na, double adjust, double kernel_width_in_sd, int n_bins, std::string bandwidth_method, double mode_frac_maxdens);

#endif

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

#ifndef NORMALIZE_H
#define NORMALIZE_H

#include <RcppArmadillo.h>

arma::vec _pairscale_normalization(arma::mat& x, const arma::uvec& clusters, std::string centeral_tendency_measure, int min_value_count, int density_npoints, double density_adjust, double density_kernel_width_in_sd, std::string bandwidth_method, double mode_frac_maxdens, double tmean_trim, double threshold_std, int niter_irls, bool check_na);

#endif

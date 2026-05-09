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

#ifndef KDE_H
#define KDE_H

#include <RcppArmadillo.h>

double fast_kde_mode(const arma::vec& x, const arma::vec& kernel);
double fast_kde_modeareamean(const arma::vec& x, arma::vec& y, const arma::vec& kernel, double frac_max_density);
arma::vec linear_binning(const arma::vec& x, int index_start, int index_endincl, double minval, double range, int n_bins, double& step);
arma::vec gauss_kernel(double bandwidth, int n, double max_dist, double bin_width);

#endif

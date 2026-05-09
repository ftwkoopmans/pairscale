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

#ifndef UTIL_H
#define UTIL_H

#include <RcppArmadillo.h>

int move_nonfinite_to_end(arma::vec& v);
void find_quantile_in_partialvector(arma::vec& x, int values_in_x, double qtl, int& index, double& value);
void find_quantile_skipn_in_partialvector(arma::vec& x, int values_in_x, double qtl, int skip_n_elements, int& index, double& value);
double median_in_partialvector(arma::vec& x, int values_in_x);
void arma_vec_minmax(const arma::vec& v, int n_valid, double& min_val, double& max_val);

#endif


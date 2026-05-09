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

#include "util.h"



// move finite values to the front, efficiently in O(N) time
// returns the number of finite values in v
int move_nonfinite_to_end(arma::vec& v) {
  auto it = std::partition(v.begin(), v.end(), [](double val) {
    return std::isfinite(val);
  });
  return (int)std::distance(v.begin(), it);
}



// compute median from an arma::vec
// this functions sorts input vector in-place
// O(N) calculation of nth percentile (avoid full sorting)
// assumes; values_in_x <= length(x)
// values_in_x indicates number of valid elements in x, starting from index 0 (i.e. subset of x may be in use)
double median_in_partialvector(arma::vec& x, int values_in_x) {
  const size_t mid_index = values_in_x / 2;
  std::nth_element(x.begin(), x.begin() + mid_index, x.begin() + values_in_x);
  double median = x[mid_index];

  // if even number of elements, median is average of mid_index and (mid_index - 1)
  if (values_in_x % 2 == 0) {
    double max_left = *std::max_element(x.begin(), x.begin() + mid_index);
    median = (median + max_left) / 2.0;
  }

  return median;
}



// efficient implementation to find min and max values in a vector in 1 pass
// there is no such function in Armadillo at the moment and calling both v.min() and v.max() is wasteful
void arma_vec_minmax(const arma::vec& v, int n_valid, double& min_val, double& max_val) {
  const arma::uword n = n_valid;
  const double* data = v.memptr();  // raw pointer; no bounds-check overhead

  min_val = std::numeric_limits<double>::max();
  max_val = std::numeric_limits<double>::lowest();

  // unroll by 4 to hint SIMD vectorization for the compiler
  arma::uword i = 0;
  for (; i + 4 <= n; i += 4) {
    const double a = data[i],   b = data[i+1];
    const double c = data[i+2], d = data[i+3];

    if (a < min_val) min_val = a;
    if (a > max_val) max_val = a;
    if (b < min_val) min_val = b;
    if (b > max_val) max_val = b;
    if (c < min_val) min_val = c;
    if (c > max_val) max_val = c;
    if (d < min_val) min_val = d;
    if (d > max_val) max_val = d;
  }

  // remaining elements
  for (; i < n; ++i) {
    const double x = data[i];
    if (x < min_val) min_val = x;
    if (x > max_val) max_val = x;
  }
}




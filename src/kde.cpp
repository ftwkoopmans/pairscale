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

#include "kde.h"

// returns index of (first occurrence of) KDE max value
double fast_kde_mode(const arma::vec& x, const arma::vec& kernel) { // , arma::vec& y, bool return_y
  double maxdens = 0;
  int maxdens_index = 0;
  int n = x.n_elem;
  int m = kernel.n_elem;
  if (n == 0 || m == 0)  {
    return arma::datum::nan;
  }

  int h = m / 2; // kernel offset = center

  // raw pointers for maximum memory access speed
  const double* px = x.memptr();
  const double* pk = kernel.memptr();

  // default/typical input data: x is longer than the kernel/mask
  if (n >= m) {
    int center_end = n + h - m + 1;

    for (int i = 0; i < n; ++i) {
      double sum = 0.0;

      if (i < h) {
        // Left boundary; mask overlaps with left out-of-bounds
        int j_start = h - i;
        for (int j = j_start; j < m; ++j) {
          sum += px[i + j - h] * pk[j];
        }
      }
      else if (i < center_end) {
        // Center region; mask is fully inside bounds
        // offset pointer for zero-overhead access (fast, but dangerous)
        const double* px_offset = px + i - h;

        // instruct the compiler to use SIMD if possible
        // only do this for the center region (sufficient datapoints to make this worthwhile / overhead tradeoff)
#pragma omp simd
        for (int j = 0; j < m; ++j) {
          sum += px_offset[j] * pk[j];
        }
      }
      else {
        // Right boundary; mask overlaps with right out-of-bounds
        int j_end = n + h - i;
        for (int j = 0; j < j_end; ++j) {
          sum += px[i + j - h] * pk[j];
        }
      }

      //if(return_y) {
      //  y[i] = sum;
      //}
      if(sum > maxdens) {
        maxdens = sum;
        maxdens_index = i;
      }
    }
  } else {
    // fallback for edge case where input is shorter than the kernel
    for (int i = 0; i < n; ++i) {
      double sum = 0.0;
      int j_start = std::max(0, h - i);
      int j_end   = std::min(m, n + h - i);
      for (int j = j_start; j < j_end; ++j) {
        sum += px[i + j - h] * pk[j];
      }
      //if(return_y) {
      //  y[i] = sum;
      //}
      if(sum > maxdens) {
        maxdens = sum;
        maxdens_index = i;
      }
    }
  }

  return maxdens_index;
}



// finds mean around the mode using the area (of x) that lies within x% density of mode
// most code analogous to fast_kde_mode()
double fast_kde_modeareamean(const arma::vec& x, arma::vec& y, const arma::vec& kernel, double frac_max_density) {
  double maxdens = 0;
  int maxdens_index = 0;
  int n = x.n_elem;
  int m = kernel.n_elem;
  if (n == 0 || m == 0)  {
    return arma::datum::nan;
  }

  int h = m / 2; // kernel offset = center

  // raw pointers for maximum memory access speed
  const double* px = x.memptr();
  const double* pk = kernel.memptr();
  //const double* py = y.memptr();

  // default/typical input data: x is longer than the kernel/mask
  if (n >= m) {
    int center_end = n + h - m + 1;

    for (int i = 0; i < n; ++i) {
      double sum = 0.0;

      if (i < h) {
        // Left boundary; mask overlaps with left out-of-bounds
        int j_start = h - i;
        for (int j = j_start; j < m; ++j) {
          sum += px[i + j - h] * pk[j];
        }
      }
      else if (i < center_end) {
        // Center region; mask is fully inside bounds
        // offset pointer for zero-overhead access (fast, but dangerous)
        const double* px_offset = px + i - h;

        // instruct the compiler to use SIMD if possible
        // only do this for the center region (sufficient datapoints to make this worthwhile / overhead tradeoff)
        for (int j = 0; j < m; ++j) {
          sum += px_offset[j] * pk[j];
        }
      }
      else {
        // Right boundary; mask overlaps with right out-of-bounds
        int j_end = n + h - i;
        for (int j = 0; j < j_end; ++j) {
          sum += px[i + j - h] * pk[j];
        }
      }

      y[i] = sum;
      if(sum > maxdens) {
        maxdens = sum;
        maxdens_index = i;
      }
    }
  } else {
    // fallback for edge case where input is shorter than the kernel
    for (int i = 0; i < n; ++i) {
      double sum = 0.0;
      int j_start = std::max(0, h - i);
      int j_end   = std::min(m, n + h - i);
      for (int j = j_start; j < j_end; ++j) {
        sum += px[i + j - h] * pk[j];
      }
      y[i] = sum;
      if(sum > maxdens) {
        maxdens = sum;
        maxdens_index = i;
      }
    }
  }

  // weighted mean within threshold
  if(frac_max_density < 1) {
    double dens_threshold = maxdens * frac_max_density;
    double area_sum = maxdens;
    double weighted_index = maxdens * maxdens_index;
    for (int i = maxdens_index + 1; i < n && y[i] > dens_threshold; ++i) {
      area_sum += y[i];
      weighted_index += i * y[i];
    }
    for (int i = maxdens_index - 1; i >= 0 && y[i] > dens_threshold; --i) {
      area_sum += y[i];
      weighted_index += i * y[i];
    }
    return weighted_index / area_sum;
  }

  return maxdens_index + 0.0;
}



// max_dist maximum number of standard deviations
arma::vec gauss_kernel(double bandwidth, int n, double max_dist, double bin_width) {
  const double inv_bandwidth = 1.0 / bandwidth;

  // at (3 or) 4 standard deviation there is barely any contribution from data points in a gaussian kernel
  // window size: stop at (3 or) 4 times bandwidth. Divide by bin_width to get distance in terms of indices over X and Y
  // min(n-1, ...) is a bugfix for edge-cases
  const int kernel_bin_width = std::min(n - 1, (int)std::ceil((max_dist * bandwidth) / bin_width));

  // kernel vector; at each element it holds the respective weight values at and around a data point
  arma::vec kernel = arma::vec(1 + kernel_bin_width * 2, arma::fill::none);

  kernel(kernel_bin_width) = 0; // center index is zero distance
  double val;
  for(int k = 1; k <= kernel_bin_width; k++) {
    // k-th bin from center; what is the distance in "standard deviations"?  (to be used as x in dnorm())
    val = k * bin_width * inv_bandwidth;
    kernel(kernel_bin_width + k) = val; // k-th element from the center, to the right
    kernel(kernel_bin_width - k) = val; // analogous for moving left from center
  }

  return arma::normpdf(kernel); // don't provide mu and sigma, default is equivalent to R's dnorm()
}



arma::vec linear_binning(const arma::vec& x, int index_start, int index_endincl, double minval, double range, int n_bins, double& step) {
  arma::vec binned(n_bins, arma::fill::zeros);
  step = range / (n_bins - 1);
  const double inv_step = 1.0 / step;
  const double* ptr_data = x.memptr();

  // linear binning
  for (int i = index_start; i <= index_endincl; ++i) {
    double val = ptr_data[i];

    double grid_idx = (val - minval) * inv_step;
    int idx = static_cast<int>(grid_idx);

    if (idx < 0) {
      binned(0) += 1.0;
    } else if (idx >= n_bins - 1) {
      binned(n_bins - 1) += 1.0;
    } else {
      double frac = grid_idx - idx;
      binned(idx)     += (1.0 - frac);
      binned(idx + 1) += frac;
    }
  }

  return binned;
}


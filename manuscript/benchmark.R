# Copyright (C) 2026 Frank Koopmans
#
# This file is part of pairscale.
#
# pairscale is free software: you can redistribute it and/or modify
# it under the terms of the GNU Affero General Public License as published
# by the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# pairscale is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
# GNU Affero General Public License for more details.
#
# You should have received a copy of the GNU Affero General Public License
# along with this program. If not, see <https://www.gnu.org/licenses/>.


output_log_filename = "~/temp/log_time.txt"
rdata_result_filename = "~/temp/benchmark_dataset"

# parse commandline arguments
args = commandArgs(trailingOnly = TRUE)
stopifnot(length(args) == 3)
param_algo = args[1]
param_nrow = args[2]
param_ncol = args[3]

# load RData file that contains the data matrix with exact dimensions and log2 tranformed yes/no for these parameters
#
# since we don't execute any other R code or perform data manipulations prior to calling the normalization functions,
# we here strictly measure computation time and RAM of the normalization functions below
param_islog2 = param_algo != "vsn"
f_cache = sprintf("%s_nrow=%s_ncol=%s_log2=%s.RData", rdata_result_filename, param_nrow, param_ncol, tolower(as.character(param_islog2)))
load(f_cache)



time_start = Sys.time()


if(param_algo == "pairscale_mean") {
  s = pairscale::pairscale_mean(data_matrix)
}

if(param_algo == "pairscale_median") {
  s = pairscale::pairscale_median(data_matrix)
}

if(param_algo == "pairscale_trimmedmean") {
  s = pairscale::pairscale_trimmedmean(data_matrix)
}

if(param_algo == "pairscale_mode_fastest") {
  s = pairscale::pairscale_mode(data_matrix, bandwidth_method = "nrd_fastest")
}

if(param_algo == "pairscale_mode_fast") {
  s = pairscale::pairscale_mode(data_matrix, bandwidth_method = "nrd_fast")
}

if(param_algo == "pairscale_mode_default") {
  s = pairscale::pairscale_mode(data_matrix, bandwidth_method = "nrd")
}

if(param_algo == "msdap_mwmb") {
  data_matrix = msdap::normalize_vwmb(data_matrix, metric_within = "mode", metric_between = "mode")
}

if(param_algo == "limma_loess") {
  data_matrix = limma::normalizeCyclicLoess(data_matrix, iterations = 10, method = "fast")
}

if(param_algo == "vsn") {
  data_matrix = vsn::justvsn(data_matrix) # note; justvsn takes non-log input
}


# report timing to console and logfile
time_stop = Sys.time()
time_delta = as.numeric(difftime(time_stop, time_start, units = "secs"))
cat(sprintf("time_delta; %.1f sec, %.1f min, %.1f hour\n", time_delta, time_delta/60, time_delta/3600))
cat(sprintf("%s\t%s\t%s\t%s\n", param_algo, param_nrow, param_ncol, time_delta), file = output_log_filename, append = TRUE)



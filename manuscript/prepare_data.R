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



### this script:
#
# 1) generate synthetic datasets for evaluation of computational performance
# - load the LFQbench dataset (PMID:27701404) bundled with MS-DAP as a readily-available source of proteomics data
#   (we're already using the msdap R package to benchmark the MWMB algorithm, so no need to download the data table separately)
# - bootstrap these data into large synthetic dataset by adding some random noise, while at the core preserving characteristics of a typical proteomics dataset
# - the generated datasets are stored as .RData files for later use in benchmarking RAM and CPU footprints
#
# 2) apply all normalization algorithms to LFQbench and compute log2fc error afterwards
#
### next steps:
# - after generating synthetic datasets using this script, run benchmark.sh (which uses the benchmark.R script)
# - finally, use figures.R to generate data visualizations using output files from this script and benchmark.sh



library(tidyverse)

# fix seed so we generate the same synthetic dataset every time we run this script
set.seed(123)




normalization_algorithm_names = tribble(
  ~from, ~to,
  "vsn", "VSN",
  "limma_loess", "Limma Loess",
  "msdap_mwmb", "MWMB (original)",
  "pairscale_mean", "pairscale mean",
  "pairscale_median", "pairscale median",
  "pairscale_trimmedmean", "pairscale tmean",
  "pairscale_mode_fastest", "pairscale mode fastest",
  "pairscale_mode_fast", "pairscale mode fast",
  "pairscale_mode_default", "pairscale mode"
)



######################################## load LFQbench dataset




rdata_result_filename = "~/temp/benchmark_dataset"
lfqbench_filename = system.file("extdata", "Skyline_HYE124_TTOF5600_64var_it2.tsv.gz", package = "msdap")

# read TSV file from disk
lfqbench_asis = read_tsv(lfqbench_filename)

# only retain precursors observed in at least 4/6 samples
lfqbench = lfqbench_asis |>
  mutate(precursor = paste(ModifiedSequence, PrecursorCharge)) |>
  filter(is.finite(TotalArea)) |>
  add_count(precursor) |>
  filter(n >= 4)

n_distinct(lfqbench$FileName)
n_distinct(lfqbench$precursor)


# convert from long-format tibble to a precursor*sample matrix with log2 abundance values
tmp = lfqbench |> pivot_wider(id_cols = "precursor", names_from = "FileName", values_from = "TotalArea")
lfqbench_matrix = as.matrix(tmp |> select(-precursor))
rownames(lfqbench_matrix) = tmp$precursor
lfqbench_matrix = log2(lfqbench_matrix)
lfqbench_matrix__condition = LETTERS[1 + grepl("007|009|011", colnames(lfqbench_matrix))]
lfqbench_matrix__proteingroup = lfqbench$ProteinName[match(tmp$precursor, lfqbench$precursor)]
lfqbench_matrix__classification = msdap::regex_classification(lfqbench_matrix__proteingroup, regex=c(background="_HUMA", foreground="_YEAS", discard="_ECOL"))
print(table(lfqbench_matrix__classification))
# proportion of NA values per FileName; colSums(!is.na(lfqbench_matrix), na.rm=T) / nrow(lfqbench_matrix)



######################################## generate synthetic datasets



# expand filtered input matrix from ~20k rows to 100k
mat_100k = matrix(NA_real_, nrow = 100000, ncol = ncol(lfqbench_matrix))
for(iter_row in 0:4) {
  # random set of 20k rows, with columns in random order
  iter_data = lfqbench_matrix[sample(seq_len(nrow(lfqbench_matrix)), 20000), sample(seq_len(ncol(lfqbench_matrix)), ncol(lfqbench_matrix))]
  # map to indices in mat_100k
  row_start = 1 + iter_row * 20000
  row_end = row_start + 20000 - 1
  for(j in 1:ncol(iter_data)) {
    mat_100k[row_start:row_end, j] = iter_data[ , j]
  }
}

cat("lfqbench_matrix fraction of NA values:", sum(is.na(lfqbench_matrix)) / length(lfqbench_matrix), "\n")
cat("mat_100k fraction of NA values:", sum(is.na(mat_100k)) / length(mat_100k), "\n")


# some cleanup to keep RAM footprint minimal
rm(lfqbench_asis, lfqbench, tmp)
gc()
gc()

# expand to large dataset; 100k rows, 4k columns
NCOL = 4000
mat = matrix(NA_real_, nrow = 100000, ncol = NCOL)
mat_col_index = 0
frac_intensity_variation = 0.1
# ballpark number for noise-level per row
sd_target = apply(mat_100k, 1, sd, na.rm = TRUE)
sd_target = sd_target + quantile(sd_target, probs = 0.1, na.rm = TRUE)


while(mat_col_index < NCOL) {
  for(j in 1:ncol(mat_100k)) {
    mat_col_index = mat_col_index + 1
    if(mat_col_index <= NCOL) {
      # measurement noise across rows + single scaling factor to simulate sample loading differences
      noise_col = runif(n = 1, min = -1, max = 1)
      noise_row = rnorm(nrow(mat_100k), mean = 0, sd = sd_target)
      mat[,mat_col_index] = mat_100k[,j] + noise_row + noise_col
    }
  }
}

cat("mat fraction of NA values:", sum(is.na(mat)) / length(mat), "\n")


for(iter_ncol in c(100, 500, 1000, 2000, 4000)) {
  for(iter_nrow in c(10000, 100000)) {
    data_matrix = mat[1:iter_nrow, 1:iter_ncol]
    save(data_matrix, file = sprintf("%s_nrow=%d_ncol=%d_log2=true.RData", rdata_result_filename, iter_nrow, iter_ncol), compression_level = 3)
    data_matrix = 2^data_matrix
    save(data_matrix, file = sprintf("%s_nrow=%d_ncol=%d_log2=false.RData", rdata_result_filename, iter_nrow, iter_ncol), compression_level = 3)
  }
}


### basic QC plot to validate above data generation; data distributions per sample/filename/column

# initial matrix that represents the input data extrapolated to 100k rows
for(j in 1:ncol(mat_100k)) {
  if(j == 1) plot(density(mat_100k[,j], na.rm = TRUE), type = "l")
  else lines(density(mat_100k[,j], na.rm = TRUE), col = j+1)
}

# plot 25 random samples from the final dataset, which includes row and column noise added to each sample/filename/column
plot(density(mat[,sample(seq_len(NCOL), size = 1)], na.rm = TRUE), type = "l", xlim = quantile(mat, probs = c(0.01, 0.99), na.rm = TRUE))
for(j in sample(seq_len(NCOL), size = 25)) {
    lines(density(mat[,j], na.rm = TRUE), col = j)
}








######################################## apply all normalization algorithms to LFQbench and compute log2fc error afterwards

norm_function_wrapper = function(data_matrix, param_algo) {
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
    data_matrix = vsn::justvsn(2^data_matrix) # note; justvsn takes non-log input
  }

  return(data_matrix)
}


compute_error = function(param_algo, mat_input, col_condition, row_classifications) {
  # add a constant to force local deep copy. Should not be needed in R, but enforce anyway
  mat_norm = norm_function_wrapper(mat_input + 0.0, param_algo)
  cond_a = mat_norm[,col_condition == "A"]
  cond_b = mat_norm[,col_condition == "B"]
  log2fc = apply(cond_b, 1, mean, na.rm = TRUE) - apply(cond_a, 1, mean, na.rm = TRUE)
  rows = is.finite(log2fc) & rowSums(is.finite(cond_a)) >= 2 & rowSums(is.finite(cond_b)) >= 2
  return(tibble(classification = row_classifications[rows], log2fc = log2fc[rows], algo = param_algo))
}


# add random scaling factors to define the 'unnormalized' data
# if we don't and the input data is already normalized, an algorithm that doesn't do anything performs well
set.seed(123)
data_pre_norm = lfqbench_matrix
for(j in 1:ncol(data_pre_norm)) {
  data_pre_norm[,j] = data_pre_norm[,j] + runif(1, min = -2, max = 2)
}

# for each normalization algorithm; apply normalization, then compute log2fc errors for background/foreground protein species
log2fc_error = bind_rows(lapply(
  normalization_algorithm_names$from,
  FUN = compute_error,
  mat_input = data_pre_norm,
  col_condition = lfqbench_matrix__condition,
  row_classifications = lfqbench_matrix__classification
)) |>
  # remove ecoli and ambiguous proteingroups, as well as the pairscale_mean approach (not enough space in the figure, but also not interesting; not accurate for omics data)
  # since it's very naive and not applicable to such data (huge errors, as expected)
  filter(classification != "discard", algo != "pairscale_mean") |>
  mutate(
    expect = as.numeric(classification != "background"),
    error = log2fc - expect,
    algo_label = replace_values(algo, from = normalization_algorithm_names$from, to = normalization_algorithm_names$to),
    algo_label = factor(algo_label, levels = unique(c(normalization_algorithm_names$to, algo_label)) )
  )


save(log2fc_error, normalization_algorithm_names, file = "~/temp/log2fc_error.RData")



# QC: plot distributions
p = ggplot(log2fc_error, aes(error, colour = algo_label)) +
  geom_density() +
  facet_wrap(.~classification) +
  coord_cartesian(c(-1,1))
print(p)

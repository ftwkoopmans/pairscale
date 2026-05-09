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


library(tidyverse)
library(patchwork)


# throughout this script, we'll assume all files are in the ~/temp dir
load("~/temp/log2fc_error.RData")




theme_custom_textsize = function() {
  theme(
    axis.text = element_text(size = 6),
    strip.text = element_text(size = 7),
    axis.title = element_text(size = 7),
    legend.text = element_text(size = 7),
    legend.title = element_text(size = 7),
    plot.title = element_text(size = 7),
    plot.subtitle = element_text(size = 6),
    plot.background = element_blank()
  )
}



plotdata = log2fc_error |> filter(algo == "pairscale_mode_default")
p_log2fc_distributions = ggplot(plotdata, aes(log2fc, colour = classification)) +
  geom_density(key_glyph = "path") +
  geom_vline(xintercept = c(0, 1), colour = "grey") +
  geom_vline(xintercept = c(median(plotdata$log2fc)), colour = "black", linetype = "dashed") +
  coord_cartesian(xlim = c(-1, 3)) +
  labs(x = "Log2fc", y = "Density") +
  theme_classic(base_size = 9) +
  theme_custom_textsize() +
  theme(
    legend.title = element_blank(),
    legend.position = "inside",
    legend.position.inside = c(1,1),
    legend.justification.inside = c(1,0.9),
    legend.margin = ggplot2::margin(),
    legend.box.margin = ggplot2::margin(),
    legend.key.spacing.y = ggplot2::unit(-2, "pt")
  )
# print(p_log2fc_distributions)



# QC: plot distributions
p = ggplot(log2fc_error, aes(error, colour = algo_label)) +
  geom_density(key_glyph = "path") +
  facet_wrap(.~classification) +
  coord_cartesian(c(-1,1))
print(p)




p_norm_error = log2fc_error |>
  group_by(algo_label, classification) |>
  summarise(mean = mean(error), sd = sd(error), .groups = "drop") |>
  ggplot(aes(x = mean, y = algo_label, colour = classification)) +
  geom_vline(xintercept = 0, colour = "darkgrey") +
  geom_linerange(aes(xmin = mean - sd, xmax = mean + sd), position = position_dodge2(0.5), show.legend = FALSE) +
  geom_point(size = 0.75, position = position_dodge2(0.5)) +
  labs(y = NULL, x = "Log2fc error") +
  theme_bw(base_size = 9) +
  theme_custom_textsize() +
  theme(
    panel.grid.major.y = element_blank(),
    legend.title = element_blank(),
    legend.position = "none"
  )
# print(p_norm_error)




########################################################################################################################


tib_time = read.delim("~/temp/log_time.txt", sep = "\t", header = FALSE, col.names = c("algo", "nrow", "ncol", "value")) |>
  group_by(algo, ncol, nrow) |>
  summarise(value = median(value), .groups = "drop") |>
  mutate(
    algo_label = replace_values(algo, from = normalization_algorithm_names$from, to = normalization_algorithm_names$to),
    algo_label = factor(algo_label, levels = unique(c(normalization_algorithm_names$to, algo_label)) )
  )

tib_ram = read.delim("~/temp/log_ram.txt", sep = "\t", header = FALSE, col.names = c("algo", "nrow", "ncol", "value")) |>
  group_by(algo, ncol, nrow) |>
  summarise(value = median(value), .groups = "drop") |>
  mutate(
    algo_label = replace_values(algo, from = normalization_algorithm_names$from, to = normalization_algorithm_names$to),
    algo_label = factor(algo_label, levels = unique(c(normalization_algorithm_names$to, algo_label)) )
  )

print(tib_time, n = Inf)
print(tib_ram, n = Inf)





plotdata = tib_time |>
  filter( ! algo %in% c("none", "pairscale_mean")) |>
  mutate(row_label = factor(sprintf("%dk rows", nrow/1000L), levels = c("10k rows", "100k rows")) )


p_time_norm = ggplot(plotdata, aes(x = log10(ncol), y = value / 60, colour = algo_label, label = ifelse(value / 3600 >= 1, sprintf("%.1f hour", value / 3600), sprintf("%.1f", value / 60)) )) +
  geom_line() +
  geom_point(size = 0.5, show.legend = FALSE) +
  geom_text(data = plotdata |> arrange(desc(ncol)) |> distinct(nrow, algo, .keep_all = TRUE), size = 6, size.unit = "pt", hjust = -0.3, show.legend = FALSE) +
  scale_x_continuous(breaks = log10(sort(unique(tib_time$ncol))), labels = sapply(sort(unique(tib_time$ncol)), function(x) ifelse(x < 1000, x, paste0(x/1000, "k"))), expand = expansion(mult = c(0.05, 0.35)) ) +
  scale_y_continuous(transform = "log10", labels = scales::comma, expand = expansion(mult = c(0.01, 0.05))) +
  scale_color_discrete("Algorithm") +
  labs(x = "Dataset size (N samples)", y = "Computation time (minutes)") +
  facet_wrap(.~row_label) +
  theme_bw(base_size = 9) +
  theme_custom_textsize() +
  theme(
    legend.position = "bottom",
    panel.grid.minor.x = element_blank()
  )
# print(p_time_norm)




plotdata = tib_ram |>
  filter( ! algo %in% c("none", "pairscale_mean")) |>
  left_join(
    tib_ram |> filter(algo == "none") |> select(nrow, ncol, refvalue = value),
    by = c("nrow", "ncol")
  ) |>
  mutate(
    row_label = factor(sprintf("%dk rows", nrow/1000L), levels = c("10k rows", "100k rows"))
    ### optionally, show RAM without overhead; subtract the amount of RAM taken up by the input data matrix
    ## value = value - refvalue
  )

p_ram_norm = ggplot(plotdata, aes(x = log10(ncol), y = (value / 10^9), colour = algo_label, label = sprintf("%.1f", value / 10^9))) +
  geom_line() +
  geom_point(size = 0.5, show.legend = FALSE) +
  geom_text(data = plotdata |> arrange(desc(ncol)) |> distinct(nrow, algo, .keep_all = TRUE), size = 6, size.unit = "pt", hjust = -0.3, show.legend = FALSE) +
  scale_x_continuous(breaks = log10(sort(unique(tib_time$ncol))), labels = sapply(sort(unique(tib_time$ncol)), function(x) ifelse(x < 1000, x, paste0(x/1000, "k"))), expand = expansion(mult = c(0.05, 0.25)) ) +
  scale_y_continuous(transform = "log10", labels = scales::comma, expand = expansion(mult = c(0.01, 0.05))) +
  scale_color_discrete("Algorithm") +
  labs(x = "Dataset size (N samples)", y = "RAM (gigabytes)") +
  facet_wrap(.~row_label) +
  theme_bw(base_size = 9) +
  theme_custom_textsize() +
  theme(
    legend.position = "bottom",
    panel.grid.minor.x = element_blank()
  )
# print(p_ram_norm)






layout <- "
ACD
BCD
"
p = patchwork::wrap_plots(
  patchwork::free(p_log2fc_distributions) + theme(plot.margin = margin()),
  p_norm_error + theme(plot.margin = margin()),
  p_time_norm, p_ram_norm + theme(legend.position = "none"),
  design = layout,
  widths = c(1, 1.5, 1.5),
  heights = c(1, 1.5)) +
  patchwork::plot_annotation(tag_levels = "A") &
  theme(plot.tag = element_text(size = 9, face = "bold"))

print(p)
ggsave(filename = "~/temp/figure.pdf", plot = p, width = 7.5, height = 3.5)


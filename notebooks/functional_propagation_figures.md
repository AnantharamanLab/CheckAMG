Functional label propagation
================
James C. Kosmopoulos
2026-09-22

``` r
knitr::opts_chunk$set(echo = TRUE, warning = FALSE, message = FALSE, fig.width = 10, fig.height = 6, dpi = 150)

library("arrow");packageVersion("arrow")
```

    ## 
    ## Attaching package: 'arrow'

    ## The following object is masked from 'package:utils':
    ## 
    ##     timestamp

    ## [1] '13.0.0'

``` r
library("tidyverse");packageVersion("tidyverse")
```

    ## ── Attaching core tidyverse packages ──────────────────────────────────────────────────────────────────────────────────────────────────────── tidyverse 2.0.0 ──
    ## ✔ dplyr     1.1.4     ✔ readr     2.1.5
    ## ✔ forcats   1.0.0     ✔ stringr   1.5.1
    ## ✔ ggplot2   3.5.2     ✔ tibble    3.2.1
    ## ✔ lubridate 1.9.4     ✔ tidyr     1.3.1
    ## ✔ purrr     1.0.4

    ## ── Conflicts ────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────── tidyverse_conflicts() ──
    ## ✖ lubridate::duration() masks arrow::duration()
    ## ✖ dplyr::filter()       masks stats::filter()
    ## ✖ dplyr::lag()          masks stats::lag()
    ## ℹ Use the conflicted package (<http://conflicted.r-lib.org/>) to force all conflicts to become errors

    ## [1] '2.0.0'

``` r
library("cowplot");packageVersion("cowplot")
```

    ## 
    ## Attaching package: 'cowplot'
    ## 
    ## The following object is masked from 'package:lubridate':
    ## 
    ##     stamp

    ## [1] '1.1.3'

``` r
library("scales");packageVersion("scales")
```

    ## 
    ## Attaching package: 'scales'
    ## 
    ## The following object is masked from 'package:purrr':
    ## 
    ##     discard
    ## 
    ## The following object is masked from 'package:readr':
    ## 
    ##     col_factor

    ## [1] '1.4.0'

``` r
library("patchwork");packageVersion("patchwork")
```

    ## 
    ## Attaching package: 'patchwork'
    ## 
    ## The following object is masked from 'package:cowplot':
    ## 
    ##     align_plots

    ## [1] '1.3.2'

``` r
library("jsonlite");packageVersion("jsonlite")
```

    ## 
    ## Attaching package: 'jsonlite'
    ## 
    ## The following object is masked from 'package:purrr':
    ## 
    ##     flatten

    ## [1] '2.0.0'

``` r
PROP_DIR  <- file.path("./tables/propagation")
NOVEL_DIR <- file.path("./tables/novel_avgs")
PLOT_DIR  <- file.path("./plots/functional_propagation")
dir.create(PLOT_DIR, showWarnings = FALSE, recursive = TRUE)

SCHEMA_VERSION <- "propagation-1.2"
MANIFEST <- jsonlite::fromJSON(file.path(PROP_DIR, "manifest.json"))
stopifnot(identical(MANIFEST$schema_version, SCHEMA_VERSION))
EVAL <- jsonlite::fromJSON(file.path(PROP_DIR, "evaluation_summary.json"))
TARGET <- MANIFEST$primary_target

prop  <- function(name) read_parquet(file.path(PROP_DIR, paste0(name, ".parquet")))
novel <- function(name) read_parquet(file.path(NOVEL_DIR, paste0(name, ".parquet")))
save_plot_both <- function(p, stem, w, h, dpi = 600, bg = "white") {
  ggsave(file.path(PLOT_DIR, paste0(stem, ".png")), p, width = w, height = h, units = "in", dpi = dpi, bg = bg)
  invisible(p)
}
```

# Palette and theme

``` r
# level colors run dark to light with depth of information, so specific reads as the strongest claim
level_lvl <- c("specific", "L1", "category", "unassigned")
level_lab <- c("specific" = "Specific function", "L1" = "L1 function category",
               "category" = "Broad category", "unassigned" = "No label at target")
level_pal <- c("specific" = "#1a9850", "L1" = "#74add1", "category" = "#fdae61", "unassigned" = "grey70")

category_pal <- c("metabolic" = "#2c7fb8", "physiological" = "#7fbc41", "regulatory" = "#d95f0e")
rank_pal <- c("calibrated" = "#1a9850", "d1_only" = "grey55")
rank_lab <- c("calibrated" = "Calibrated probability", "d1_only" = "Embedding distance only")

# tier scales: cluster tier purple, embedding tier the same green the calibrated curves use
reach_lvl <- c("both", "cluster_only", "embedding_only", "neither")
reach_lab <- c("both" = "Both tiers", "cluster_only" = "Cluster tier only",
               "embedding_only" = "Embedding tier only", "neither" = "Neither tier")
method_lvl <- c("cluster_only", "embedding_only", "tiered")
method_lab <- c("cluster_only" = "Cluster tier alone", "embedding_only" = "Calibrated embedding alone",
                "tiered" = "Deployed tiered rule")
method_pal <- c("cluster_only" = "#762a83", "embedding_only" = "#1a9850", "tiered" = "#0f3b5f")
source_lvl <- c("cluster", "embedding", "none")
source_lab <- c("cluster" = "Cluster tier", "embedding" = "Calibrated embedding tier",
                "none" = "No label at target")
source_pal <- c("cluster" = "#762a83", "embedding" = "#1a9850", "none" = "grey70")

base_theme <- theme_cowplot(font_size = 12) +
  theme(axis.text = element_text(size = 11), axis.title = element_text(size = 13),
        legend.title = element_text(size = 13, face = "bold"), legend.text = element_text(size = 11),
        strip.text = element_text(size = 12, face = "bold"),
        legend.key.size = unit(0.42, "cm"), panel.grid.major = element_blank())
theme_set(base_theme)
```

# Method comparison at the deployed target

# Propagation coverage and precision at different sequence identity thresholds, for each propagation tier

``` r
tstrat <- prop("fig_tier_by_identity") %>%
  mutate(stratum = factor(stratum, levels = c("<50%", "50-90%", ">=90%", "all"),
                          labels = c("< 50%", "50-90%", ">= 90%", "All")),
         method = factor(method, levels = method_lvl))
xlabs <- tstrat %>% distinct(stratum, n_report) %>%
  mutate(lab = paste0(stratum, "\nn = ", comma(n_report))) %>% arrange(stratum)
dg <- position_dodge(width = 0.78)

p_prop_cov_precision <- ggplot(tstrat, aes(stratum, coverage, fill = method)) +
  geom_col(position = dg, width = 0.72, alpha = 0.45, color = NA) +
  geom_line(aes(y = precision, color = method, group = method), position = dg, linewidth = 0.7) +
  geom_point(aes(y = precision, color = method), position = dg, size = 2.4) +
  scale_fill_manual(values = method_pal, labels = method_lab, name = NULL) +
  scale_color_manual(values = method_pal, labels = method_lab, name = NULL) +
  scale_x_discrete(labels = xlabs$lab) +
  scale_y_continuous(labels = percent_format(accuracy = 1), limits = c(0, 1),
                     expand = expansion(c(0, 0.02))) +
  labs(x = "Seq. identity to nearest training protein",
       y = "Coverage (bars) and\nprecision (points)") +
  theme(legend.position = "bottom", legend.direction = "vertical", legend.location = "plot",
        plot.subtitle = element_text(size = 9.5), axis.text.x = element_text(size = 9.5),
        axis.title.y = element_text(size = 11.5),
        plot.margin = margin(t = 12, r = 4, b = 2, l = 4, unit = "pt"))
print(p_prop_cov_precision)
```

# Which propagation tier could reach each held-out protein

``` r
ovl <- prop("fig_tier_overlap")
stopifnot(nrow(ovl) == 12, setequal(ovl$reachable_by, reach_lvl))
N_REPORT <- round(unique(ovl$n / (ovl$pct_of_slice / 100))[1])
pc1 <- function(x) sprintf("%.1f", 100 * x)

p_propagation_reach <- ovl %>%
  mutate(level = factor(level, levels = c("specific", "L1", "category")),
         reachable_by = factor(reachable_by, levels = rev(reach_lvl)),
         lab = paste0(comma(n), " (", pc1(pct_of_slice / 100), "%)\n",
                      case_when(
                        reachable_by == "both" ~
                          paste0("C ", pc1(tier1_precision), " | E ", pc1(tier2_precision)),
                        reachable_by == "cluster_only" ~
                          paste0("C ", pc1(tier1_precision), " | 1NN ", pc1(tier2_precision)),
                        reachable_by == "embedding_only" ~ paste0("E ", pc1(tier2_precision)),
                        TRUE ~ paste0("1NN ", pc1(tier2_precision))))) %>%
  ggplot(aes(level, reachable_by, fill = pct_of_slice)) +
  geom_tile(color = "white", linewidth = 0.5) +
  geom_text(aes(label = lab), size = 3.1, lineheight = 0.95) +
  scale_fill_gradient(low = "#f7f7f7", high = "#9ecae1", name = "Share of held-out proteins",
                      labels = function(x) paste0(x, "%"), guide = guide_colorbar(barwidth = 8)) +
  scale_x_discrete(labels = c("specific" = "Specific", "L1" = "L1 category",
                              "category" = "Broad category")) +
  scale_y_discrete(labels = reach_lab) +
  labs(x = paste0("Level at which a tier could assign, ", comma(N_REPORT), " held-out proteins"),
       y = NULL,
       subtitle = paste0("C = cluster tier precision; E = admitted embedding precision;\n1NN = precision the uncalibrated neighbor would have delivered")) +
  theme(
    legend.position = "bottom",
    plot.subtitle = element_text(size = 9.5),
    axis.text.y = element_text(size = 10.5)
  )
print(p_propagation_reach)
```

# Validation controls of label propagation

``` r
ctrl <- prop("fig_validation_controls")
ctrl$panel <- gsub("Contribution of finetuning", "Embedding space", ctrl$panel)
ctrl <- ctrl %>%
  mutate(panel = factor(panel, levels = c("Unseen held-out proteins", "Unseen genomic context",
                                          "Specificity", "Embedding space"),
                        labels = c("Unseen held-out proteins\n(category precision)",
                                   "Unseen genomic context\n(category precision)",
                                   "Specificity\n(confident-call rate)",
                                   "Embedding space\n(category precision)")),
         condition = fct_rev(fct_inorder(condition)))

p_ctrl <- ggplot(ctrl, aes(value, condition, fill = panel)) +
  geom_col(width = 0.68, color = "black", linewidth = 0.25) +
  geom_text(aes(label = ifelse(value < 0.01, sprintf("%.4f%%", 100 * value),
                               sprintf("%.1f%%", 100 * value))),
            hjust = -0.08, size = 3) +
  facet_grid(panel ~ ., scales = "free_y", space = "free_y", switch = "y") +
  scale_x_continuous(labels = percent, expand = expansion(mult = c(0, 0.12)), breaks = c(0, 0.25, 0.5, 0.75, 1)) +
  scale_fill_brewer(palette = "Set2", guide = "none") +
  labs(x = "Proportion of proteins", y = NULL) +
  theme(strip.placement = "outside", strip.text.y.left = element_text(angle = 0, size = 9, face = "bold"),
        axis.text.y = element_text(size = 9), panel.spacing.y = unit(0.25, "lines"))

print(p_ctrl)
```

![](functional_propagation_figures_files/figure-gfm/plot-propagation-validation-controls-1.png)<!-- -->

# What propagation delivers on the applied auxiliary viral genes across prediction targets

``` r
hi_sweep <- prop("fig_target_sweep_applied")

p_prop_comp_by_target <- hi_sweep %>%
  mutate(level = factor(level, levels = level_lvl)) %>%
  ggplot(aes(target, pct, fill = level)) +
  geom_col(width = 0.045) +
  scale_fill_manual(values = level_pal, labels = level_lab, name = NULL) +
  scale_y_continuous(labels = function(x) paste0(x, "%"), expand = expansion(c(0, 0.02))) +
  scale_x_continuous(breaks = unique(hi_sweep$target), labels = percent_format(accuracy = 1)) +
  labs(x = "Precision target", y = "Sequence-similarity invisible AVGs") +
  theme(
    axis.text.x = element_text(size = 10, angle = 45, hjust = 1),
    legend.position = "bottom",
    legend.justification = "center",
    legend.location = "plot",
    legend.direction = "vertical"
    ) +
  guides(fill = guide_legend(nrow = 3))

print(p_prop_comp_by_target)
```

# Assignment depth of each propagation tier for sequence similarity invisible AVG slices

``` r
SUB_LVL <- c("All invisible AVGs", "Already identified", "Unresolved")
nm <- novel("novel_master_per_protein") %>% select(known_function, label_source, final_level)
app <- bind_rows(
    nm %>% mutate(subset = SUB_LVL[1]),
    nm %>% filter(known_function)  %>% mutate(subset = SUB_LVL[2]),
    nm %>% filter(!known_function) %>% mutate(subset = SUB_LVL[3])) %>%
  count(subset, final_level, label_source, name = "n") %>%
  group_by(subset) %>% mutate(n_subset = sum(n), frac = n / n_subset) %>% ungroup() %>%
  mutate(final_level = factor(final_level, levels = level_lvl),
         label_source = factor(label_source, levels = source_lvl))

stopifnot(identical(
  app %>% filter(subset == SUB_LVL[1]) %>%
    transmute(level = as.character(final_level), label_source = as.character(label_source), n) %>%
    arrange(level, label_source) %>% pull(n),
  TIER_APP %>% transmute(level, label_source, n = as.integer(n)) %>%
    arrange(level, label_source) %>% pull(n)))

sub_labs <- app %>% distinct(subset, n_subset) %>%
  mutate(lab = paste0(subset, "\nn = ", comma(n_subset)))
sub_labs <- setNames(sub_labs$lab, sub_labs$subset)
app <- app %>% mutate(subset = factor(subset, levels = SUB_LVL))

p_tier_assignment_depth <- ggplot(app, aes(frac, fct_rev(final_level), fill = label_source)) +
  geom_col(width = 0.72, color = "black", linewidth = 0.2) +
  geom_text(aes(label = ifelse(frac >= 0.02, comma(n), ""),
                # vjust = case_when(label_source == "cluster" ~ -0.28,
                #                   label_source == "embedding" ~ 1.28,
                #                   TRUE ~ 0.5)
                ),
            position = position_stack(vjust = 0.5), size = 2.6) +
  facet_wrap(~ subset, nrow = 1, labeller = as_labeller(sub_labs)) +
  scale_fill_manual(values = source_pal, labels = source_lab, name = NULL) +
  scale_y_discrete(labels = c("specific" = "Specific", "L1" = "L1", "category" = "Broad",
                              "unassigned" = "None")) +
  scale_x_continuous(labels = percent_format(accuracy = 1), expand = expansion(c(0.03, 0.04)),
                     breaks = c(0, 0.2, 0.4)) +
  coord_cartesian(clip = "off") +
  coord_flip() +
  labs(y = "Assignment depth at the deployed target", x = "Share of the subset") +
  theme(legend.position = "bottom", legend.location = "plot", legend.direction = "vertical",
        strip.text = element_text(size = 8.8), strip.clip = "off",
        plot.subtitle = element_text(size = 9),
        panel.spacing.x = unit(0.7, "lines"), axis.text.x = element_text(size = 9))

print(p_tier_assignment_depth)
```

# Precision of the uncalibrated neighbor against sequence identity

``` r
joint <- prop("fig_identity_by_distance")
strat <- prop("fig_ladder_by_identity")

p_precision_vs_seq_id <- strat %>%
  pivot_longer(c(category, L1, specific), names_to = "level", values_to = "precision") %>%
  mutate(level = factor(level, levels = c("specific", "L1", "category"))) %>%
  ggplot(aes(band, precision, color = level, group = level)) +
  geom_line(linewidth = 0.9) + geom_point(size = 2) +
  scale_color_manual(values = level_pal, labels = level_lab, name = NULL) +
  scale_y_continuous(labels = percent_format(accuracy = 1), limits = c(0, 1)) +
  labs(x = "Seq. identity to nearest training protein", y = "Raw embedding 1-NN precision") +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    legend.position = "bottom",
    # legend.location = "plot",
    legend.direction = "vertical",
    legend.text = element_text(size = 9.5),
    legend.key.size = unit(0.32, "cm"),
    plot.margin = margin(t = 2, r = 4, b = 8, l = 4, unit = "pt")
  )

print(p_precision_vs_seq_id)
```

# Realized precision vs. precision target at each propagation level

``` r
sweep <- prop("fig_target_sweep_heldout")

p_realized_v_target <- sweep %>%
  filter(!is.na(realized_precision), level != "unassigned") %>%
  mutate(level = factor(level, levels = level_lvl)) %>%
  ggplot(aes(target, realized_precision, color = level)) +
  geom_abline(slope = 1, intercept = 0, linetype = 2, color = "grey60") +
  geom_line(linewidth = 0.8) + geom_point(size = 2) +
  scale_color_manual(values = level_pal, labels = level_lab, name = NULL) +
  scale_x_continuous(labels = percent_format(accuracy = 1)) +
  scale_y_continuous(labels = percent_format(accuracy = 1)) +
  labs(x = "Precision target", y = "Realized precision") +
  theme(legend.position = "bottom", legend.direction = "vertical", legend.location = "plot", plot.margin = margin(l = 0, r =))

print(p_realized_v_target)
```

# Precision vs. propagation coverage at each level

``` r
rc <- prop("fig_risk_coverage")

p_precison_vs_coverage <- rc %>%
  mutate(level = factor(level, levels = c("specific", "L1", "category")),
         ranking = factor(ranking, levels = names(rank_pal))) %>%
  ggplot(aes(coverage, precision, color = ranking)) +
  geom_hline(yintercept = TARGET, linetype = 2, color = "grey60") +
  geom_line(linewidth = 0.9) +
  facet_wrap(~ level, nrow = 1, labeller = as_labeller(
    c("specific" = "Specific function", "L1" = "L1 function category", "category" = "Broad category"))) +
  scale_color_manual(values = rank_pal, labels = rank_lab, name = NULL) +
  scale_x_continuous(labels = percent_format(accuracy = 1)) +
  scale_y_continuous(labels = percent_format(accuracy = 1), limits = c(0, 1)) +
  labs(x = "Coverage", y = "Precision") +
  theme(legend.position = "bottom")

print(p_precison_vs_coverage)
```

# Accuraacy of the calibrated probabilities at each level

``` r
rel <- prop("fig_reliability")
ece <- tibble(level = names(EVAL$calibration_ece), ece = unlist(EVAL$calibration_ece)) %>%
  mutate(lab = paste0("ECE = ", format(round(ece, 4), nsmall = 4)))

p_accuracy_vs_prob <- rel %>%
  mutate(level = factor(level, levels = c("specific", "L1", "category"))) %>%
  ggplot(aes(confidence, accuracy)) +
  geom_abline(slope = 1, intercept = 0, linetype = 2, color = "grey60") +
  geom_point(aes(size = n), color = "#1a9850", alpha = 0.8) +
  geom_line(color = "#1a9850", linewidth = 0.6) +
  geom_text(data = ece %>% mutate(level = factor(level, levels = c("specific", "L1", "category"))),
            aes(x = 0.03, y = 0.96, label = lab), hjust = 0, size = 4, inherit.aes = FALSE) +
  facet_wrap(~ level, nrow = 1, labeller = as_labeller(
    c("specific" = "Specific function", "L1" = "L1 function category", "category" = "Broad category"))) +
  scale_size_continuous(range = c(0.8, 4.5), guide = "none") +
  scale_x_continuous(labels = percent_format(accuracy = 1), limits = c(0, 1)) +
  scale_y_continuous(labels = percent_format(accuracy = 1), limits = c(0, 1)) +
  labs(x = "Calibrated probability", y = "Observed accuracy")

print(p_accuracy_vs_prob)
```

# Composite figure

``` r
p_row1 <- ((p_prop_cov_precision | p_propagation_reach | p_ctrl) + plot_layout(widths = c(4.25, 6, 5.75)))
p_row2 <- ((p_prop_comp_by_target | p_tier_assignment_depth | p_precision_vs_seq_id) + plot_layout(widths = c(4.25, 7.25, 4)))
p_row3 <- ((p_realized_v_target | p_precison_vs_coverage | p_accuracy_vs_prob) + plot_layout(widths = c(2, 7, 7))) &
  theme(axis.text.x = element_text(size = 10))

p_comp <- (p_row1 / p_row2 / p_row3) + plot_annotation(tag_levels = "A")
p_comp <- p_comp +
  plot_layout(heights = c(1, 1, 1)) &
  theme(plot.tag = element_text(size = 16, face = "bold", hjust = 0, vjust = 1),
        plot.tag.position = "topleft", plot.tag.location = "margin")

save_plot_both(p_comp, "functional_propagation_supplement", w = 16, h = 16)
print(p_comp)
```

![](functional_propagation_figures_files/figure-gfm/composite-1.png)<!-- -->

# Supplemental table

``` r
supp <- prop("fig_method_summary")
knitr::kable(supp, digits = 4)
```

| level | quantity | value | note |
|:---|:---|---:|:---|
| category | raw 1-NN micro precision | 0.8952 | deployed rule, distinct held-out proteins |
| category | raw 1-NN macro precision | 0.8491 | unweighted over 3 labels |
| L1 | raw 1-NN micro precision | 0.5866 | deployed rule, distinct held-out proteins |
| L1 | raw 1-NN macro precision | 0.5595 | unweighted over 28 labels |
| specific | raw 1-NN micro precision | 0.3813 | deployed rule, distinct held-out proteins |
| specific | raw 1-NN macro precision | 0.2255 | unweighted over 6,181 labels |
| L1 | coverage at 0.90 precision | 0.3364 | calibrated ranking |
| L1 | AURC | 0.8065 | calibrated ranking |
| category | coverage at 0.90 precision | 0.9182 | calibrated ranking |
| category | AURC | 0.9501 | calibrated ranking |
| specific | coverage at 0.90 precision | 0.2339 | calibrated ranking |
| specific | AURC | 0.6801 | calibrated ranking |
| L1 | coverage at 0.90 precision | 0.0276 | distance-only baseline |
| category | coverage at 0.90 precision | 0.8162 | distance-only baseline |
| specific | coverage at 0.90 precision | 0.0001 | distance-only baseline |
| specific | expected calibration error | 0.0077 | reporting slice, 20 bins |
| L1 | expected calibration error | 0.0055 | reporting slice, 20 bins |
| category | expected calibration error | 0.0057 | reporting slice, 20 bins |
| L1 | precision, rare label tercile | 0.2975 | n=1,136 |
| category | precision, rare label tercile | 0.8574 | n=1,136 |
| specific | precision, rare label tercile | 0.0537 | n=1,136 |

# Source Data and Supplemental Tables

``` r
SRC_DIR <- file.path(PROP_DIR, "source_data")
dir.create(SRC_DIR, showWarnings = FALSE, recursive = TRUE)

readr::write_csv(sweep,    file.path(SRC_DIR, "edfig9A_target_sweep.csv"))
readr::write_csv(rc,       file.path(SRC_DIR, "edfig9B_risk_coverage.csv"))
readr::write_csv(rel,      file.path(SRC_DIR, "edfig9C_reliability.csv"))
readr::write_csv(strat,    file.path(SRC_DIR, "edfig9F_precision_by_identity_band.csv"))
readr::write_csv(ovl,      file.path(SRC_DIR, "edfig9D_tier_reach_overlap.csv"))
readr::write_csv(tstrat,   file.path(SRC_DIR, "edfig9E_tier_by_identity.csv"))
readr::write_csv(app,      file.path(SRC_DIR, "edfig9G_applied_label_source.csv"))
readr::write_csv(hi_sweep, file.path(SRC_DIR, "edfig9H_applied_target_sweep.csv"))
readr::write_csv(ctrl,     file.path(SRC_DIR, "edfig9I_validation_controls.csv"))

readr::write_csv(joint,    file.path(SRC_DIR, "edfig9_extra_identity_by_distance.csv"))
cat(sprintf("source data written to %s\n", SRC_DIR))
```

    ## source data written to ./tables/propagation/source_data

``` r
print(tibble::tibble(
  file = list.files(SRC_DIR, pattern = "^edfig9"),
  rows = vapply(list.files(SRC_DIR, pattern = "^edfig9", full.names = TRUE),
                function(p) nrow(readr::read_csv(p, show_col_types = FALSE)), integer(1))))
```

    ## # A tibble: 10 × 2
    ##    file                                    rows
    ##    <chr>                                  <int>
    ##  1 edfig9_extra_identity_by_distance.csv     23
    ##  2 edfig9A_target_sweep.csv                  36
    ##  3 edfig9B_risk_coverage.csv               2424
    ##  4 edfig9C_reliability.csv                   51
    ##  5 edfig9D_tier_reach_overlap.csv            12
    ##  6 edfig9E_tier_by_identity.csv              12
    ##  7 edfig9F_precision_by_identity_band.csv     8
    ##  8 edfig9G_applied_label_source.csv          20
    ##  9 edfig9H_applied_target_sweep.csv          36
    ## 10 edfig9I_validation_controls.csv            9

``` r
supp_dir <- PROP_DIR

lad  <- prop("fig_ladder_micro_macro"); ladf <- prop("fig_ladder_by_frequency")
ladi <- prop("fig_ladder_by_identity"); rcs <- prop("fig_risk_coverage_summary")
ops  <- prop("fig_operating_points")

summary_table <- bind_rows(
  lad %>% pivot_longer(c(micro, macro), names_to = "metric", values_to = "value") %>%
    transmute(section = "performance", level, stratum_type = "overall", stratum = "all",
              metric = paste0(metric, "_precision"), value, n),
  lad %>% transmute(section = "performance", level, stratum_type = "overall", stratum = "all",
                    metric = "n_distinct_labels", value = as.numeric(n_labels), n),
  ladf %>% transmute(section = "performance", level, stratum_type = "reference_frequency_tercile",
                     stratum = freq_tercile, metric = "precision", value = precision, n),
  ladi %>% pivot_longer(c(category, L1, specific), names_to = "level", values_to = "value") %>%
    transmute(section = "performance", level, stratum_type = "identity_to_training",
              stratum = band, metric = "precision", value, n),
  tibble(level = names(EVAL$calibration_ece), value = unlist(EVAL$calibration_ece)) %>%
    transmute(section = "performance", level, stratum_type = "overall", stratum = "all",
              metric = "expected_calibration_error", value, n = NA_integer_),
  ops %>% filter(abs(target - MANIFEST$primary_target) < 1e-9) %>%
    pivot_longer(c(min_seq_id, p_threshold, coverage, precision),
                 names_to = "metric", values_to = "value") %>%
    transmute(section = "performance", level, stratum_type = "operating_point",
              stratum = sprintf("target_%.2f_%s", target, tier), metric, value, n = n_report),
  rcs %>% pivot_longer(c(aurc, `cov@0.90`, `cov@0.80`, `cov@0.70`),
                       names_to = "metric", values_to = "value") %>%
    transmute(section = "performance", level, stratum_type = "risk_coverage",
              stratum = ranking, metric, value, n = NA_integer_),
  prop("fig_validation_controls") %>%
    transmute(section = "validation_control", level = NA_character_, stratum_type = panel,
              stratum = condition, metric, value, n = NA_integer_),
  prop("fig_target_sweep_heldout") %>%
    pivot_longer(c(rate_heldout, realized_precision), names_to = "metric", values_to = "value") %>%
    transmute(section = "target_sweep", level, stratum_type = "held-out test proteins",
              stratum = sprintf("target_%.2f", target), metric, value, n),
  prop("fig_target_sweep_applied") %>%
    transmute(section = "target_sweep", level, stratum_type = "sequence-similarity invisible AVGs",
              stratum = sprintf("target_%.2f", target), metric = "rate", value = pct / 100, n)
) %>% relocate(section, level, stratum_type, stratum, metric, value, n)

readr::write_csv(summary_table, file.path(PROP_DIR, "propagation_evaluation.csv"))
```

# Session info

``` r
sessionInfo()
```

    ## R version 4.3.3 (2024-02-29)
    ## Platform: x86_64-conda-linux-gnu (64-bit)
    ## Running under: Ubuntu 20.04.6 LTS
    ## 
    ## Matrix products: default
    ## BLAS/LAPACK: /storage2/scratch/kosmopoulos/miniconda3/envs/peatlands_env/lib/libopenblasp-r0.3.21.so;  LAPACK version 3.9.0
    ## 
    ## locale:
    ##  [1] LC_CTYPE=en_US.UTF-8       LC_NUMERIC=C              
    ##  [3] LC_TIME=en_US.UTF-8        LC_COLLATE=en_US.UTF-8    
    ##  [5] LC_MONETARY=en_US.UTF-8    LC_MESSAGES=en_US.UTF-8   
    ##  [7] LC_PAPER=en_US.UTF-8       LC_NAME=C                 
    ##  [9] LC_ADDRESS=C               LC_TELEPHONE=C            
    ## [11] LC_MEASUREMENT=en_US.UTF-8 LC_IDENTIFICATION=C       
    ## 
    ## time zone: America/Chicago
    ## tzcode source: system (glibc)
    ## 
    ## attached base packages:
    ## [1] stats     graphics  grDevices utils     datasets  methods   base     
    ## 
    ## other attached packages:
    ##  [1] jsonlite_2.0.0  patchwork_1.3.2 scales_1.4.0    cowplot_1.1.3  
    ##  [5] lubridate_1.9.4 forcats_1.0.0   stringr_1.5.1   dplyr_1.1.4    
    ##  [9] purrr_1.0.4     readr_2.1.5     tidyr_1.3.1     tibble_3.2.1   
    ## [13] ggplot2_3.5.2   tidyverse_2.0.0 arrow_13.0.0   
    ## 
    ## loaded via a namespace (and not attached):
    ##  [1] utf8_1.2.4         generics_0.1.3     stringi_1.8.7      hms_1.1.3         
    ##  [5] digest_0.6.37      magrittr_2.0.3     evaluate_1.0.3     grid_4.3.3        
    ##  [9] timechange_0.3.0   RColorBrewer_1.1-3 fastmap_1.2.0      textshaping_0.3.7 
    ## [13] cli_3.6.4          rlang_1.2.0        crayon_1.5.3       bit64_4.6.0-1     
    ## [17] withr_3.0.2        yaml_2.3.10        parallel_4.3.3     tools_4.3.3       
    ## [21] tzdb_0.5.0         assertthat_0.2.1   vctrs_0.6.5        R6_2.6.1          
    ## [25] lifecycle_1.0.4    bit_4.6.0          vroom_1.6.5        ragg_1.3.3        
    ## [29] pkgconfig_2.0.3    pillar_1.10.2      gtable_0.3.6       glue_1.8.0        
    ## [33] systemfonts_1.2.1  xfun_0.52          tidyselect_1.2.1   knitr_1.50        
    ## [37] farver_2.1.2       htmltools_0.5.8.1  rmarkdown_2.29     labeling_0.4.3    
    ## [41] compiler_4.3.3

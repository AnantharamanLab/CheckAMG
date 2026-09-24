LGBM figures
================
James C. Kosmopoulos
2026-09-10

``` r
knitr::opts_chunk$set(echo = TRUE, warning = FALSE, message = FALSE, fig.width = 10, fig.height = 6, dpi = 150)

library("arrow");packageVersion("arrow")
```

    ## Warning: package 'arrow' was built under R version 4.4.3

    ## 
    ## Attaching package: 'arrow'

    ## The following object is masked from 'package:utils':
    ## 
    ##     timestamp

    ## [1] '22.0.0.1'

``` r
library("tidyverse");packageVersion("tidyverse")
```

    ## Warning: package 'tibble' was built under R version 4.4.3

    ## Warning: package 'tidyr' was built under R version 4.4.3

    ## Warning: package 'purrr' was built under R version 4.4.3

    ## ── Attaching core tidyverse packages ──────────────────────── tidyverse 2.0.0 ──
    ## ✔ dplyr     1.1.4     ✔ readr     2.1.5
    ## ✔ forcats   1.0.0     ✔ stringr   1.5.1
    ## ✔ ggplot2   3.5.2     ✔ tibble    3.3.1
    ## ✔ lubridate 1.9.3     ✔ tidyr     1.3.2
    ## ✔ purrr     1.2.1

    ## ── Conflicts ────────────────────────────────────────── tidyverse_conflicts() ──
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
library("patchwork");packageVersion("patchwork")
```

    ## Warning: package 'patchwork' was built under R version 4.4.1

    ## 
    ## Attaching package: 'patchwork'
    ## 
    ## The following object is masked from 'package:cowplot':
    ## 
    ##     align_plots

    ## [1] '1.3.2'

``` r
library("scales");packageVersion("scales")
```

    ## Warning: package 'scales' was built under R version 4.4.1

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
library("RColorBrewer");packageVersion("RColorBrewer")
```

    ## [1] '1.1.3'

``` r
library("svglite");packageVersion("svglite")
```

    ## [1] '2.1.3'

``` r
library("ggh4x");packageVersion("ggh4x")
```

    ## Warning: package 'ggh4x' was built under R version 4.4.1

    ## [1] '0.3.1'

``` r
TABLES_DIR  <- file.path("./tables/lgbm")
PLOT_DIR <- file.path("./plots/lgbm")
dir.create(PLOT_DIR, showWarnings = FALSE, recursive = TRUE)

save_plot_both <- function(p, stem, w, h, dpi = 600) {
  ggsave(file.path(PLOT_DIR, paste0(stem, ".png")), p,
         width = w, height = h, units = "in", dpi = dpi, bg = "white")
  # ggsave(file.path(PLOT_DIR, paste0(stem, ".svg")), p,
  #        width = w, height = h, units = "in", bg = "white", device = svglite::svglite)
  invisible(p)
}

read_if <- function(f) {
  if (file.exists(f)) suppressMessages(readr::read_tsv(f, show_col_types = FALSE, progress = FALSE)) else NULL
}
```

# Palettes and theme

``` r
dataset_pal <- c(
  "Training"              = "grey35",
  "Virus enriched"        = "#5E4FA2",
  "Near all virus"        = "#3288BD",
  "Half viral/host"       = "#66C2A5",
  "Equal viral/nonviral"  = "#ABDDA4",
  "Training distribution" = "#E6F598",
  "Equal viral/MGE/host"  = "#FEE08B",
  "Host enriched"         = "#FDAE61",
  "Integrated proviruses" = "#F46D43",
  "MGE enriched"          = "#D53E4F",
  "Near all host"         = "#9E0142"
)

confidence_shapes <- c(
  "High confidence"   = 21,
  "Medium confidence" = 22,
  "Low confidence"    = 25
)

confidence_order  <- c("High confidence", "Medium confidence", "Low confidence")

strict_pal <- c("True" = "#377EB8", "False" = "#E41A1C", "Either" = "#4DAF4A")

source_order <- c("Virus", "MGE", "Host")

# Training has no curve in the performance panels, so drop = FALSE leaves an empty key;
# this invisible layer supplies the glyph, made opaque again by the guide's override.aes
training_key_layer <- geom_line(
  data = tibble(
    x = c(0, 1), y = c(0, 1),
    dataset = factor("Training", levels = names(dataset_pal))
  ),
  aes(x = x, y = y, color = dataset),
  alpha = 0,
  inherit.aes = FALSE
)

base_theme <- theme_cowplot(font_size = 12) +
  theme(plot.title    = element_text(face = "bold", size = 14),
        axis.text     = element_text(size = 12),
        axis.title    = element_text(size = 14),
        strip.background = element_rect(fill = "grey92", color = NA),
        strip.text    = element_text(face = "bold"),
        legend.position = "right",
        legend.title    = element_text(face = "bold"),
        panel.grid.major.x = element_blank(),
        panel.grid.major.y = element_blank(),
        )
theme_set(base_theme)
```

# LGBM performance

## Precision-recall curves

``` r
pr_points <- read_parquet(file.path(TABLES_DIR, "lgbm_pr_points.parquet")) %>%
  mutate(
    dataset = factor(dataset, levels =  names(dataset_pal)),
    confidence = str_replace(confidence, "Confidence", tolower("Confidence"))
  )

p_pr <- ggplot() +
  geom_hline(yintercept = 0.90, linetype = 3, color = "grey40") +
  geom_hline(yintercept = 0.95, linetype = 2, color = "grey40") +
  geom_line(
    data = pr_points %>% filter(type == "curve"),
    aes(x = recall, y = precision, color = dataset),
    alpha = 0.7
  ) +
  geom_point(
    data = pr_points %>% filter(type == "point"),
    aes(x = recall, y = precision, fill = dataset, shape = confidence),
    color = "black",
    alpha = 0.8,
    size = 3,
    stroke = 0.7
  ) +
  training_key_layer +
  scale_color_manual(values = dataset_pal, name = "Dataset", drop = FALSE) +
  scale_fill_manual(values = dataset_pal, guide = "none") +
  scale_shape_manual(values = confidence_shapes, name = "Confidence cutoff") +
  guides(
    color = guide_legend(order = 1, override.aes = list(alpha = 1, linewidth = 1)), # dataset first
    shape = guide_legend(order = 2, override.aes = list(alpha = 1))  # confidence second
  ) +
  labs(x = "Recall", y = "Precision") +
  base_theme
save_plot_both(p_pr, "lgbm_precision_recall_plot", w = 6, h = 4)
print(p_pr)
```

![](lgbm_figures_files/figure-gfm/pr-curves-1.png)<!-- -->

## ROC curves

``` r
roc_points <- read_parquet(file.path(TABLES_DIR, "lgbm_roc_points.parquet")) %>%
  mutate(
    dataset = factor(dataset, levels =  names(dataset_pal)),
    confidence = str_replace(confidence, "Confidence", tolower("Confidence"))
  )

p_roc <- ggplot() +
  geom_abline(slope = 1, xintercept = 0, yintercept = 0, linetype = 2, color = "grey40") +
  geom_line(
    data = roc_points %>% filter(type == "curve"),
    aes(x = fpr, y = tpr, color = dataset),
    alpha = 0.7
  ) +
  geom_point(
    data = roc_points %>% filter(type == "point"),
    aes(x = fpr, y = tpr, fill = dataset, shape = confidence),
    color = "black",
    alpha = 0.8,
    size = 3,
    stroke = 0.7
  ) +
  training_key_layer +
  scale_color_manual(values = dataset_pal, name = "Dataset", drop = FALSE) +
  scale_fill_manual(values = dataset_pal, guide = "none") +
  scale_shape_manual(values = confidence_shapes, name = "Confidence cutoff") +
  guides(
    color = guide_legend(order = 1, override.aes = list(alpha = 1, linewidth = 1)), # dataset first
    shape = guide_legend(order = 2, override.aes = list(alpha = 1))  # confidence second
  ) +
  labs(x = "False positive rate", y = "True positive rate") +
  base_theme
save_plot_both(p_roc, "lgbm_roc_plot", w = 6, h = 4)
print(p_roc)
```

![](lgbm_figures_files/figure-gfm/roc-curves-1.png)<!-- -->

## True viral fraction over model probabilotiy (“reliability”)

``` r
reliability_points <- read_parquet(file.path(TABLES_DIR, "lgbm_reliability_points.parquet")) %>%
  mutate(
    dataset = factor(dataset, levels =  names(dataset_pal))
  )

p_rel <- ggplot() +
  geom_abline(slope = 1, xintercept = 0, yintercept = 0, linetype = 2, color = "grey40") +
  geom_line(
    data = reliability_points,
    aes(x = predicted_prob, y = observed_fraction, color = dataset),
    alpha = 0.7
  ) +
  geom_point(
    data = reliability_points,
    aes(x = predicted_prob, y = observed_fraction, fill = dataset),
    shape = 21,
    color = "black",
    alpha = 0.8,
    size = 2,
    stroke = 0.7
  ) +
  training_key_layer +
  scale_color_manual(values = dataset_pal, name = "Dataset", drop = FALSE) +
  scale_fill_manual(values = dataset_pal, guide = "none") +
  guides(
    color = guide_legend(override.aes = list(alpha = 1, linewidth = 1)),
  ) +
  labs(x = "LightGBM predicted probability", y = "True viral fraction") +
  base_theme
save_plot_both(p_rel, "lgbm_reliability_plot", w = 6, h = 4)
print(p_rel)
```

![](lgbm_figures_files/figure-gfm/reliability-curves-1.png)<!-- -->

# Train/test dataset composition

``` r
composition <- read_if(file.path(TABLES_DIR, "train_test_dataset_composition.tsv")) %>%
  mutate(
    dataset = factor(dataset, levels = names(dataset_pal)),
    unit = factor(unit, levels = c("Proteins", "Scaffolds")),
    source = factor(source, levels = source_order)
  ) %>%
  pivot_longer(
    cols = c(percent, count),
    names_to = "metric",
    values_to = "value"
  ) %>%
  mutate(
    metric = factor(metric, levels = c("percent", "count"))
  ) %>%
  arrange(unit, metric, dataset, source)

comp_dodge <- position_dodge(width = 0.8)

comp_layers <- function() {
  list(
    geom_line(linewidth = 0.7, alpha = 0.7, position = comp_dodge),
    geom_point(
      aes(fill = dataset), shape = 21, color = "black",
      size = 2.6, stroke = 0.5, alpha = 0.9, position = comp_dodge
    ),
    scale_color_manual(values = dataset_pal, name = "Dataset", drop = FALSE),
    scale_fill_manual(values = dataset_pal, guide = "none"),
    base_theme,
    theme(
      panel.grid.major.y = element_line(color = "grey90"),
      legend.position = "bottom"
    )
  )
}

no_x <- theme(axis.text.x = element_blank(), axis.ticks.x = element_blank(), axis.title.x = element_blank())

count_protein <- composition %>% filter(metric == "count", unit == "Proteins")

low_max  <- count_protein %>% filter(dataset != "Training") %>% pull(value) %>% max(na.rm = TRUE)
high_min <- count_protein %>% filter(dataset == "Training") %>% pull(value) %>% min(na.rm = TRUE)

bottom_lim <- low_max * 1.2
top_lim    <- high_min * 0.9

pct_labels   <- scales::label_percent(scale = 1)
count_labels <- scales::label_number(scale = 1e-6, suffix = "M", accuracy = 0.1)

pct_breaks <- seq(0, 100, 25)

pct_scale <- scale_y_continuous(
  limits = c(0, 105),
  breaks = pct_breaks,
  labels = scales::label_percent(scale = 1),
  expand = expansion(mult = c(0, 0))
)

# Percent row
p_pct_proteins <- ggplot(
  composition %>% filter(metric == "percent", unit == "Proteins"),
  aes(x = source, y = value, color = dataset, group = dataset)
) +
  comp_layers() +
  facet_wrap(~ unit) +
  pct_scale +
  labs(x = NULL, y = "Percent of dataset") +
  no_x

p_pct_scaffolds <- ggplot(
  composition %>% filter(metric == "percent", unit == "Scaffolds"),
  aes(x = source, y = value, color = dataset, group = dataset)
) +
  comp_layers() +
  facet_wrap(~ unit) +
  pct_scale +
  labs(x = NULL, y = NULL) +
  no_x

# Count / Scaffolds
p_count_scaffolds <- ggplot(
  composition %>% filter(metric == "count", unit == "Scaffolds"),
  aes(x = source, y = value, color = dataset, group = dataset)
) +
  comp_layers() +
  scale_y_continuous(labels = count_labels, expand = expansion(mult = c(0.02, 0.08))) +
  labs(x = "Sequence source", y = NULL)

# Count / Proteins, top
p_count_proteins_top <- ggplot(
  count_protein,
  aes(x = source, y = value, color = dataset, group = dataset)
) +
  comp_layers() +
  coord_cartesian(ylim = c(top_lim, NA)) +
  scale_y_continuous(breaks = scales::pretty_breaks(n = 3), labels = count_labels, expand = expansion(mult = c(0.1, 0.15))) +
  labs(x = NULL, y = NULL) +
  no_x +
  theme(plot.margin = margin(t = 5.5, r = 5.5, b = 2, l = 5.5))

# Count / Proteins, bottom
p_count_proteins_bottom <- ggplot(
  count_protein,
  aes(x = source, y = value, color = dataset, group = dataset)
) +
  comp_layers() +
  coord_cartesian(ylim = c(0, bottom_lim)) +
  scale_y_continuous(breaks = scales::pretty_breaks(n = 4), labels = count_labels, expand = expansion(mult = c(0.02, 0.05))) +
  labs(x = "Sequence source", y = "Count (millions)") +
  theme(plot.margin = margin(t = 2, r = 5.5, b = 5.5, l = 5.5))

design <- "
AB
CD
ED
"

p_comp <- p_pct_proteins + p_pct_scaffolds +
  p_count_proteins_top + p_count_scaffolds + p_count_proteins_bottom +
  plot_layout(guides = "collect", design = design, heights = c(2, 0.7, 1.3)) &
  theme(
    legend.position       = "bottom",
    legend.box            = "horizontal",
    legend.box.just       = "center",
    legend.justification  = c(0.5, 0),
    legend.title.position = "top",
    legend.title.align    = 0,
    legend.text           = element_text(size = 12),
    legend.title          = element_text(size = 14, face = "bold")
  )

print(p_comp)
```

![](lgbm_figures_files/figure-gfm/dataset-composition-1.png)<!-- -->

``` r
save_plot_both(p_comp, "lgbm_dataset_composition", w = 10, h = 6.5)
```

## Composite figure (LGBM)

``` r
composite <- (wrap_elements(full = (p_comp & theme(legend.position = "none"))) / (p_roc | p_pr | p_rel)) +
  plot_annotation(tag_levels = 'A') +
  plot_layout(guides = "collect", heights = c(6.5, 3.5)) &
  theme(plot.tag = element_text(size = 16)) &
  scale_fill_manual(values = dataset_pal, guide = "none") &
  guides(
    color = guide_legend(order = 1, nrow = 4, override.aes = list(alpha = 1, linewidth = 1)),
    shape = guide_legend(order = 2, nrow = 4, override.aes = list(alpha = 1))
  ) &
  theme(
    legend.position       = "bottom",
    legend.box            = "horizontal",
    legend.box.just       = "center",
    legend.justification  = c(0.5, 0),
    legend.title.position = "top",
    legend.title.align    = 0,
    legend.text           = element_text(size = 12),
    legend.title          = element_text(size = 14, face = "bold")
  )

save_plot_both(composite, "lgbm_composite", w = 10, h = 10)
print(composite)
```

![](lgbm_figures_files/figure-gfm/composite-1.png)<!-- -->

# Strict/ambiguous viral region predictions vs LGBM

## Compare performace of confidence levels when stratified by strict viral region labels

``` r
benchmark_strict <- read_if(file.path(TABLES_DIR, "lgbm_benchmark_summary_table_strict_strat.tsv")) %>%
  mutate(
    Dataset = factor(Dataset, levels =  names(dataset_pal)),
    Confidence_Level = str_replace(Confidence_Level, "Confidence", tolower("Confidence"))
  ) %>%
  mutate(Confidence_Level = factor(Confidence_Level, levels = confidence_order)) %>%
  mutate(Strict_Viral_Region = factor(Strict_Viral_Region, levels = c("True", "False", "Either")))

num_cols_after_threshold <- benchmark_strict %>%
  select(which(names(.) == "Threshold"):last_col()) %>%
  select(-Threshold) %>%
  select(where(is.numeric)) %>%
  names()

benchmark_long <- benchmark_strict %>%
  pivot_longer(
    cols = all_of(num_cols_after_threshold),
    names_to = "Metric",
    values_to = "Value"
  ) %>%
  distinct()

metrics_to_plot <- c(
  "TPR", "TNR", "Precision", "Recall", "Predictions_Made"
)

new_labels <- c(
  "TPR" = "Sensitivity", "TNR" = "Specificity", "Precision" = "Precision", "Recall" = "Recall", "Predictions_Made" = "N predictions",
  "High confidence" = "High confidence", "Medium confidence" = "Medium confidence", "Low confidence" = "Low confidence"
)

benchmark_plotting <- benchmark_long %>%
  filter(Metric %in% metrics_to_plot) %>%
  mutate(Metric = factor(Metric, levels = metrics_to_plot))

p_strict_compare <- ggplot() +
  geom_col(
    data = benchmark_plotting,
    aes(x=Dataset, y = Value, fill = Strict_Viral_Region),
    stat = "identity",
    position = position_dodge2(padding = 0.2),
    color = "black",
    
  ) +
  facet_grid(
    Metric ~ Confidence_Level,
    labeller = as_labeller(new_labels),
    scales = "free_y"
  ) +
  facetted_pos_scales(
    y = list(
      Metric == "Predictions_Made" ~ scale_y_log10()
    )
  ) +
  scale_fill_manual(values = strict_pal) +
  guides(
    fill = guide_legend(title = "Strict viral region enforced")
  ) +
  base_theme +
  theme(
    axis.text.x           = element_text(angle = 30, hjust = 1),
    panel.grid.major.y    = element_line(color = "grey60"),
    legend.position       = "bottom",
    legend.box            = "horizontal",
    legend.box.just       = "center",
    legend.justification  = c(0.5, 0),
    legend.title.position = "left",
    legend.title.align    = 0,
    legend.text           = element_text(size = 12),
    legend.title          = element_text(size = 14, face = "bold")
  )
  labs(x = "Dataset", y = "Value")
```

    ## $x
    ## [1] "Dataset"
    ## 
    ## $y
    ## [1] "Value"
    ## 
    ## attr(,"class")
    ## [1] "labels"

``` r
save_plot_both(p_strict_compare, "lgbm_strict_viral_region_compare", w = 12, h = 8)
print(p_strict_compare)
```

![](lgbm_figures_files/figure-gfm/strict-regions-compare-1.png)<!-- -->

# Session info

``` r
sessionInfo()
```

    ## R version 4.4.0 (2024-04-24)
    ## Platform: aarch64-apple-darwin20
    ## Running under: macOS 26.6.2
    ## 
    ## Matrix products: default
    ## BLAS:   /Library/Frameworks/R.framework/Versions/4.4-arm64/Resources/lib/libRblas.0.dylib 
    ## LAPACK: /Library/Frameworks/R.framework/Versions/4.4-arm64/Resources/lib/libRlapack.dylib;  LAPACK version 3.12.0
    ## 
    ## locale:
    ## [1] en_US.UTF-8/en_US.UTF-8/en_US.UTF-8/C/en_US.UTF-8/en_US.UTF-8
    ## 
    ## time zone: America/Chicago
    ## tzcode source: internal
    ## 
    ## attached base packages:
    ## [1] stats     graphics  grDevices utils     datasets  methods   base     
    ## 
    ## other attached packages:
    ##  [1] ggh4x_0.3.1        svglite_2.1.3      RColorBrewer_1.1-3 scales_1.4.0      
    ##  [5] patchwork_1.3.2    cowplot_1.1.3      lubridate_1.9.3    forcats_1.0.0     
    ##  [9] stringr_1.5.1      dplyr_1.1.4        purrr_1.2.1        readr_2.1.5       
    ## [13] tidyr_1.3.2        tibble_3.3.1       ggplot2_3.5.2      tidyverse_2.0.0   
    ## [17] arrow_22.0.0.1    
    ## 
    ## loaded via a namespace (and not attached):
    ##  [1] utf8_1.2.4        generics_0.1.3    stringi_1.8.4     hms_1.1.3        
    ##  [5] digest_0.6.37     magrittr_2.0.4    evaluate_1.0.5    grid_4.4.0       
    ##  [9] timechange_0.3.0  fastmap_1.2.0     fansi_1.0.6       textshaping_0.4.0
    ## [13] cli_3.6.5         crayon_1.5.3      rlang_1.1.7       bit64_4.5.2      
    ## [17] withr_3.0.2       yaml_2.3.10       parallel_4.4.0    tools_4.4.0      
    ## [21] tzdb_0.4.0        assertthat_0.2.1  vctrs_0.6.5       R6_2.6.1         
    ## [25] lifecycle_1.0.5   bit_4.5.0         vroom_1.6.5       ragg_1.5.1       
    ## [29] pkgconfig_2.0.3   pillar_1.9.0      gtable_0.3.6      glue_1.8.0       
    ## [33] systemfonts_1.3.1 highr_0.11        xfun_0.48         tidyselect_1.2.1 
    ## [37] rstudioapi_0.16.0 knitr_1.48        farver_2.1.2      htmltools_0.5.8.1
    ## [41] labeling_0.4.3    rmarkdown_2.28    compiler_4.4.0

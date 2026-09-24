Visualize strict/ambiguous viral regions
================
James C. Kosmopoulos
2026-08-04

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

    ## ── Attaching core tidyverse packages ────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────── tidyverse 2.0.0 ──
    ## ✔ dplyr     1.1.4     ✔ readr     2.1.5
    ## ✔ forcats   1.0.0     ✔ stringr   1.5.1
    ## ✔ ggplot2   3.5.2     ✔ tibble    3.2.1
    ## ✔ lubridate 1.9.4     ✔ tidyr     1.3.1
    ## ✔ purrr     1.0.4

    ## ── Conflicts ──────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────── tidyverse_conflicts() ──
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

    ## [1] '2.1.2'

``` r
library("gggenomes");packageVersion("gggenomes")
```

    ## gggenomes v1.1.3
    ## 
    ## If you use 'gggenomes' in published research, please cite:
    ## 
    ## Hackl T, Ankenbrand M, van Adrichem B, Wilkins D, Haslinger K (2024).
    ## "gggenomes: effective and versatile visualizations for comparative
    ## genomics." _arXiv_. doi:10.48550/arXiv.2411.13556
    ## <https://doi.org/10.48550/arXiv.2411.13556>.
    ## 
    ## Attaching package: 'gggenomes'
    ## 
    ## The following object is masked from 'package:dplyr':
    ## 
    ##     pick
    ## 
    ## The following object is masked from 'package:graphics':
    ## 
    ##     layout

    ## [1] '1.1.3'

``` r
library("ggnewscale");packageVersion("ggnewscale")
```

    ## [1] '0.5.2'

``` r
TABLES_DIR  <- file.path("./tables/genome_viz")
PLOT_DIR <- file.path("./plots/genome_viz")
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

safe_pull <- function(df) {
  x <- pull(df, contig)
  if (length(x) == 0) NA_character_ else x[[1]]
}
```

# Load genes/genome data

``` r
seqs <- read_parquet(file.path(TABLES_DIR, "genes_sequence_table.parquet"))
genes <- read_parquet(file.path(TABLES_DIR, "genes_feature_table.parquet"))
```

# Candidate selection

Pick candidates:

1.  A whole viral contig (from geNomad)
2.  A “whole viral” contig from proGenomes with contamination
3.  A contig containing a proviral sequence integrated in a host
    chromosome (without relabeled proteins)
4.  A contig with very mixed composition

and

5.  A contig that shows the progression of each step in the strict viral
    region algorithm

``` r
cand1 <- seqs %>%
  filter(region_contig_type == "standalone_virus") %>%
  filter(dataset == "geNomad") %>%
  filter(pct_ptns_strict_viral_region > 0.75) %>%
  filter(length >= 50000) %>%
  arrange(length) %>%
  pull(contig) %>%
  first()

cand2 <- seqs %>%
  filter(region_contig_type == "standalone_virus") %>%
  filter(dataset == "progenomes") %>%
  filter(pct_ptns_strict_viral_region > 0.5) %>%
  filter(pct_ptns_relabeled > 0.1) %>%
  filter(length >= 50000) %>%
  arrange(length) %>%
  pull(contig) %>%
  first()

cand3 <- seqs %>%
  filter(region_contig_type == "chromosome_mixed") %>%
  filter(n_viral_regions >= 1) %>%
  filter(length > 50000) %>%
  filter(pct_true_viral < 0.75 & pct_true_viral > 0.30) %>%
  filter(pct_source_host_ptns > 0.5) %>%
  arrange(length) %>%
  pull(contig) %>%
  first()

cand4 <- seqs %>%
  filter(region_contig_type == "chromosome_mixed") %>%
  filter(n_viral_regions < 3) %>%
  filter(pct_source_viral_ptns > 0.2) %>%
  filter(pct_source_mge_ptns > 0.2) %>%
  # filter(pct_source_host_ptns > 0.2) %>%
  # filter(pct_true_viral > 0.5) %>%
  filter(pct_ptns_strict_viral_region > 0.25 & pct_ptns_strict_viral_region < 0.8) %>%
  filter(length > 100000) %>%
  arrange(length) %>%
  pull(contig) %>%
  first()

# Pick a contig where every consecutive pair of displayed steps actually changes some gene's membership
displayed_step_cols <- c(
  "step1_seed_strict_wvl",
  "step2_in_refined_core",
  "step3_in_walkback_extended",
  "step4_in_snapped_region",
  "step5_in_merged_region"
)

# Per-contig: did each consecutive pair of displayed steps actually differ anywhere along the contig?
diff_pair <- function(a, b) {
  a <- coalesce(a, FALSE)
  b <- coalesce(b, FALSE)
  any(a != b)
}

step_summary <- genes %>%
  group_by(contig) %>%
  summarize(
    n_genes = n(),
    n_seed  = sum(step1_seed_strict_wvl,      na.rm = TRUE),
    n_core  = sum(step2_in_refined_core,      na.rm = TRUE),
    n_walk  = sum(step3_in_walkback_extended, na.rm = TRUE),
    n_snap  = sum(step4_in_snapped_region,    na.rm = TRUE),
    n_merge = sum(step5_in_merged_region,     na.rm = TRUE),
    diff_seed_core  = diff_pair(step1_seed_strict_wvl,      step2_in_refined_core),
    diff_core_walk  = diff_pair(step2_in_refined_core,      step3_in_walkback_extended),
    diff_walk_snap  = diff_pair(step3_in_walkback_extended, step4_in_snapped_region),
    diff_snap_merge = diff_pair(step4_in_snapped_region,    step5_in_merged_region),
    .groups = "drop"
  )

cand5 <- step_summary %>%
  filter(n_seed  >= 1,                # has a seed region
         n_merge >= 1,                # a final merged region exists
         diff_seed_core,              # 1a -> 2b changes
         diff_core_walk,              # 2b -> 3b changes
         diff_walk_snap,              # 3b -> 4 changes
         diff_snap_merge,             # 4  -> 5 changes
         between(n_genes, 20, 55)) %>%
  inner_join(select(seqs, contig), by = "contig") %>%
  arrange(desc(n_genes), desc(n_merge)) %>%
  pull(contig) %>%
  first()
```

# Genome plot building functions

## Data-prep helpers

``` r
V_SCORE_NAME        <- "V-score"
WIN_VL_SCORE_NAME   <- "Window-average\nVL-score"
V_SCORE_LIMITS      <- c(0, 10)
WIN_VL_SCORE_LIMITS <- c(0, 3)
V_SCORE_OPT         <- "D" # viridis
WIN_VL_SCORE_OPT    <- "D"

DB_LABELS <- c(kegg = "KEGG", pfam = "Pfam", phrog = "PHROG")

REGION_COLORS <- c("TRUE" = "#2166AC", "FALSE" = "#D73027")

# Seqs track: optional window restricts the displayed region (for focus view)
make_seqs_gg <- function(seqs_df, contig_id, window = NULL) {
  s <- seqs_df %>%
    filter(contig == contig_id) %>%
    transmute(seq_id = contig, length)
  if (!is.null(window))
    s <- mutate(s, start = window$start, end = window$end)
  s
}
```

## Gene track

Gene track for geom_gene:

- `name_col` picks which annotation column drives gene labels, falls
  back to Function/top_hit_description when NULL
- `score_override` lets the caller supply a precomputed score vector
  (used by the step plot, where the score is the per-gene max across
  DBs)
- Label selection is priority-ranked then greedily pruned to enforce a
  minimum spacing in bp, so labels are spread across the plot rather
  than clumping in annotation-dense regions:
  - priority 4 — genes with V-score \>= label_score_high (strong viral
    hit)
  - priority 4 — genes with V-score \< label_score_low (strong non-viral
    hit)
  - priority 3 — genes flanking a strict-viral-region boundary
    (transition in `label_region_col`, default step5_in_merged_region)
  - priority 1 — any other gene with an annotation (fills sparse zones)
- `label_target_per_kb` controls how many labels per kb of plotted
  region, the cap is bounded between `label_min_n` and `label_max_n`.
  Min spacing in bp is plotted_length / chosen_n

``` r
make_gene_track <- function(genes_df, score_col = NULL, score_override = NULL,
                             name_col = NULL,
                             label_score_high = 10,
                             label_score_low  = 2,
                             label_region_col = "step5_in_merged_region",
                             label_target_per_kb = 0.75,
                             label_min_n = 6L,
                             label_max_n = 16L) {
  ord <- order(genes_df$contig_pos_start)
  g   <- genes_df[ord, , drop = FALSE]

  name_src <- if (!is.null(name_col)) g[[name_col]]
              else                      coalesce(g$Function, g$top_hit_description)
  name_src <- str_trunc(name_src, 36, ellipsis = "...")

  score <- if (!is.null(score_override)) score_override
           else if (!is.null(score_col)) g[[score_col]]
           else NA_real_

  region_vals <- g[[label_region_col]]
  region_vals[is.na(region_vals)] <- FALSE
  prev_v <- c(NA, region_vals[-length(region_vals)])
  next_v <- c(region_vals[-1L], NA)
  is_boundary <- (region_vals != prev_v & !is.na(prev_v)) |
                 (region_vals != next_v & !is.na(next_v))

  is_high  <- !is.na(score) & score >= label_score_high
  is_low   <- !is.na(score) & score <  label_score_low
  has_name <- !is.na(name_src) & nchar(name_src) > 0L

  # Per-gene label priority (higher = picked first by the greedy pass)
  priority <- integer(nrow(g))
  priority[has_name]                       <- 1L
  priority[has_name & is_boundary]         <- 3L
  priority[has_name & (is_high | is_low)]  <- 4L

  cand_idx <- which(priority > 0L)
  if (length(cand_idx) > 0L) {
    midpts    <- (g$contig_pos_start[cand_idx] + g$contig_pos_end[cand_idx]) / 2
    total_len <- max(g$contig_pos_end) - min(g$contig_pos_start)
    n_target  <- max(label_min_n,
                     min(label_max_n,
                         as.integer(round(total_len / 1000 * label_target_per_kb))))
    min_gap   <- if (n_target > 0L) total_len / n_target else 0

    # Tie-break: prefer scores far from neutral (more informative arrow color)
    extr      <- ifelse(is.na(score[cand_idx]), 0, abs(score[cand_idx] - 5))
    ord_cands <- order(-priority[cand_idx], -extr, midpts)

    selected <- logical(length(cand_idx))
    chosen   <- numeric(0L)
    for (k in ord_cands) {
      m <- midpts[k]
      if (length(chosen) == 0L || min(abs(chosen - m)) >= min_gap) {
        selected[k] <- TRUE
        chosen <- c(chosen, m)
        if (length(chosen) >= n_target) break
      }
    }

    final_label <- logical(nrow(g))
    final_label[cand_idx[selected]] <- TRUE
  } else {
    final_label <- logical(nrow(g))
  }

  df  <- tibble(
    seq_id = g$contig,
    start  = g$contig_pos_start,
    end    = g$contig_pos_end,
    strand = if_else(g$frame >= 0, 1L, -1L),
    name   = if_else(final_label, name_src, NA_character_)
  )
  df$score   <- score
  df$introns <- vector("list", nrow(df))
  df
}
```

## Strict viral region feat track and focus window

- `region_col` selects which step-prefixed column drives the T/F
  coloring
- The step plot passes a different column per step so the band shifts as
  the algorithm walks through its stages
- Consecutive same-value genes are merged so the band is continuous
  across runs
- Focus window: center on the viral regions with `flank` non-viral genes
  on each side
  - Returns NULL when all genes are the same `Strict viral region`

``` r
make_region_feats <- function(genes_df, region_col = "step5_in_merged_region") {
  ord <- order(genes_df$contig_pos_start)
  g   <- genes_df[ord, , drop = FALSE]
  if (nrow(g) == 0L) {
    return(tibble(seq_id = character(), start = integer(),
                  end = integer(), `Strict viral region` = character()))
  }
  vals <- g[[region_col]]
  vals[is.na(vals)] <- FALSE

  rl      <- rle(vals)
  end_i   <- cumsum(rl$lengths)
  start_i <- end_i - rl$lengths + 1L

  tibble(
    seq_id                = g$contig[start_i],
    start                 = g$contig_pos_start[start_i],
    end                   = g$contig_pos_end[end_i],
    `Strict viral region` = as.character(rl$values)
  )
}

get_focus_window <- function(genes_df, contig_id, flank = 8L) {
  g <- genes_df %>%
    filter(contig == contig_id) %>%
    arrange(contig_pos_start)

  viral_idx    <- which(g$`Strict viral region` == TRUE)
  nonviral_idx <- which(g$`Strict viral region` != TRUE)
  if (length(viral_idx) == 0L || length(nonviral_idx) == 0L) return(NULL)

  first_v <- max(1L,      min(viral_idx) - flank)
  last_v  <- min(nrow(g), max(viral_idx) + flank)

  list(
    start = g$contig_pos_start[[first_v]],
    end   = g$contig_pos_end[[last_v]]
  )
}
```

## Shared track-panel builder

Each panel is its own gggenomes plot

``` r
build_track_panel <- function(gene_trk, src_feats, s,
                               scale_name, scale_option, scale_limits,
                               subtitle           = NULL,
                               # `subtitle_position`:
                               # "top"  — render `subtitle` as plot.subtitle (top-left, default)
                               # "left" — render `subtitle` as a horizontal y-axis title on the left
                               subtitle_position  = c("top", "left"),
                               # `left_label`: independent left y-axis metric label shown like DB names
                               # (KEGG/Pfam/PHROG) in the genome map. When set, always shown as the
                               # y-axis title regardless of subtitle_position. subtitle still goes to
                               # its own position (top subtitle or y-label) as normal
                               left_label         = NULL,
                               show_seq           = FALSE,
                               show_xaxis         = FALSE,
                               show_labels        = FALSE,
                               show_region_legend = FALSE,
                               # When FALSE the fill colorbar is suppressed, set FALSE on all but
                               # one designated panel so guides are collected only once
                               show_fill_legend   = TRUE,
                               pad_for_labels     = show_labels,
                               # Extra vertical room above the gene track for staggered labels
                               label_expand_y     = 2.1,
                               # ggrepel repulsion strength; increase for denser label zones
                               repel_force        = 0.4,
                               repel_force_pull   = 0.1,
                               # Font size for gene-name annotation labels (ggrepel)
                               annot_size         = 2.0,
                               # Font size for the step subtitle (plot.subtitle)
                               subtitle_size      = 7.5,
                               # Upper y-limit for ggrepel label placement
                               ylim_upper         = NULL,
                               ylim_lower         = NULL,
                               # How far to nudge labels upward from the gene track anchor
                               nudge_y_val        = 0.25) {
  subtitle_position <- match.arg(subtitle_position)

  # vjust for the left y-axis title: aligns the label with the gene track (y = 1)
  # Base value 0.27 was calibrated for label_expand_y = 2.1 (genome-map default)
  # When the panel is taller (step-plot labels) we scale down proportionally
  y_title_vjust <- if (pad_for_labels) {
    max(0.10, 0.27 * 2.1 / label_expand_y)
  } else {
    0.5
  }

  # Show the y-axis title style when using subtitle_position == "left" (genome map)
  # OR when an explicit left_label is provided (step plot)
  show_y_title <- (!is.null(subtitle) && subtitle_position == "left") || !is.null(left_label)

  gg <- gggenomes(seqs = s, genes = gene_trk, feats = list(src = src_feats))
  p  <- gg

  if (show_seq) p <- p + geom_seq()

  p <- p +
    geom_feat(
      aes(color = `Strict viral region`), data = feats(src),
      alpha = 0.5, linewidth = 10, position = "identity",
      show.legend = show_region_legend
    ) +
    scale_color_manual(
      values = REGION_COLORS, name = "Strict viral region",
      breaks = c("TRUE", "FALSE"),
      labels = c("True", "False"),
      guide  = if (show_region_legend)
                 guide_legend(order = 1, override.aes = list(linewidth = 4, alpha = 0.8))
               else "none"
    ) +
    geom_gene(aes(fill = score), color = "black", linewidth = 0.2) +
    scale_fill_viridis_c(
      name = scale_name, option = scale_option,
      limits = scale_limits, oob = squish, na.value = "grey85",
      guide = if (show_fill_legend)
                guide_colorbar(order = 2, barwidth = 0.5, barheight = 3)
              else "none"
    ) +
    theme_gggenomes_clean() +
    theme(
      legend.position = "right",
      legend.box      = "vertical",
      legend.key.size = unit(0.38, "cm"),
      legend.text     = element_text(size = 7),
      legend.title    = element_text(size = 8, face = "bold"),
      plot.subtitle   = element_text(size = subtitle_size, color = "grey20",
                                     margin = margin(b = 1)),
      axis.title.x    = if (show_xaxis) element_text(size = 8, color = "grey30") else element_blank(),
      axis.text.x     = if (show_xaxis) element_text(size = 7) else element_blank(),
      axis.ticks.x    = if (show_xaxis) element_line() else element_blank(),
      axis.line.x     = if (show_xaxis) element_line() else element_blank(),
      axis.title.y    = if (show_y_title)
                          element_text(size = 9, face = "bold", color = "grey15",
                                       angle = 0, hjust = 1, vjust = y_title_vjust,
                                       margin = margin(r = 4))
                        else element_blank(),
      plot.margin     = margin(0, 2, 0, 2, "mm")
    )

  # subtitle: goes to plot.subtitle (top) or y-axis label (left DB-name style)
  if (!is.null(subtitle)) {
    if (subtitle_position == "left") {
      p <- p + labs(y = subtitle)
    } else {
      p <- p + labs(subtitle = subtitle)
    }
  }
  # left_label overrides/supplements the y-axis label (step-plot metric names)
  if (!is.null(left_label)) {
    p <- p + labs(y = left_label)
  }
  if (show_xaxis) p <- p + labs(x = "Contig position (bp)")

  # Pad y so this panel reserves vertical space for staggered labels
  if (pad_for_labels) p <- p + expand_limits(y = label_expand_y)

  if (show_labels) {
    label_data <- gene_trk[!is.na(gene_trk$name), , drop = FALSE]
    if (nrow(label_data) > 0L) {
      p <- p +
        ggrepel::geom_text_repel(
          aes(x = (start + end) / 2, y = 1, label = name),
          data           = label_data,
          inherit.aes    = FALSE,
          size           = annot_size,
          angle          = 0,
          direction      = "y",
          # Ceiling == label_expand_y so labels fill the full allocated space
          ylim           = c(if (is.null(ylim_lower)) label_expand_y - 0.3 else ylim_lower, if (is.null(ylim_upper)) label_expand_y + 0.8 else ylim_upper),
          nudge_y        = nudge_y_val,
          force          = repel_force,
          force_pull     = repel_force_pull,
          segment.size   = 0.25,
          segment.color  = "grey40",
          segment.alpha  = 0.7,
          min.segment.length = 0,
          box.padding    = unit(0.1, "lines"),
          point.padding  = unit(0.03, "lines"),
          max.overlaps   = Inf
        )
    }
  }
  p
}
```

# Algorithm step progression plot

- 8 steps x 2 metrics (V-score / window-avg VL)
- All genes are shown with their full per-gene max-across-DBs score
- Only the region band shifts per step, by switching `region_col` so its
  TRUE/FALSE boundaries follow the algorithm’s evolving region call
- Layout overview:
  - Left column: Steps 1-2 (each step = V-score row + WV-score row
    stacked)
  - Right column: Steps 3-5 (same stacked pair, no gene-name labels)
  - Within each step the V and WV rows are tightly packed and a thin
    plot_spacer() separates successive steps
  - Metric names (“V-score” / “Window-average-score”) appear as
    left-side y-axis titles on every track
  - Legends are collected once at the top level

``` r
build_step_plot <- function(contig_id, seqs, genes) {
  if (is.na(contig_id)) return(ggplot() + theme_void())

  g_all <- genes %>%
    filter(contig == contig_id) %>%
    arrange(contig_pos_start)

  s <- make_seqs_gg(seqs, contig_id)

  max_v  <- pmax(g_all$`KEGG_V-score`,
                 g_all$`Pfam_V-score`,
                 g_all$`PHROG_V-score`,
                 na.rm = TRUE)
  max_wv <- pmax(g_all$`window_avg_KEGG_VL-score`,
                 g_all$`window_avg_Pfam_VL-score`,
                 g_all$`window_avg_PHROG_VL-score`,
                 na.rm = TRUE)

  # V track: more labels allowed, step 1 carries annotations, extra labels fill the label zone
  # WV track: same settings (labels only shown on step 1 V, not here, but track needs names for geom)
  gene_v  <- make_gene_track(g_all, score_override = max_v,
                              label_target_per_kb = 0.67,
                              label_min_n = 4L, label_max_n = 24L)
  gene_wv <- make_gene_track(g_all, score_override = max_wv,
                              label_target_per_kb = 0.5,
                              label_min_n = 4L, label_max_n = 12L)

  step_defs <- list(
    list(col = "step1_seed_strict_wvl",      label = "Step 1: Seed (window-average VL-score \u2265 3)"),
    list(col = "step2_in_refined_core",      label = "Step 2: Refined core (V-score = 10 gate)"),
    list(col = "step3_in_walkback_extended", label = "Step 3: Walkback extended (relaxed window-average VL-score)"),
    list(col = "step4_in_snapped_region",    label = "Step 4: Snapped region (all three database V-score = 10)"),
    list(col = "step5_in_merged_region",     label = "Step 5: Final merged regions (discard small regions, merge close regions)")
  )

  LEFT_STEPS  <- 1:2
  RIGHT_STEPS <- 3:5

  # Row-height constants
  # STEP_V_H + STEP_WV_H = one "step unit"; the step1 label zone adds one more
  # step unit, so step1 total = 2 step units (matching any 2 right-column steps)
  HEADER_H    <- 0.38   # col-header row (tall enough to separate header from subtitle)
  STEP_V_H    <- 0.90   # V-score track (no labels)
  STEP_WV_H   <- 0.90   # WV-score track
  SPACER_H    <- 0.25   # gap between successive steps
  # Label zone height = one full step unit, so step1_V patchwork height = 1 step
  # + label zone = STEP_V_H + STEP_V_H + STEP_WV_H = 2.28
  STEP1_ANNOT_H <- STEP_V_H + STEP_WV_H # label space = one step unit above track
  STEP1_V_H     <- STEP_V_H + STEP1_ANNOT_H

  # label_expand_y for step1: fills the full patchwork row tightly
  STEP1_EXPAND_Y <- STEP1_V_H

  # Only the very first panel (step 1 V-score, left column) shows legends,
  # all others suppress them so plot_layout(guides="collect") yields one set
  make_step_pair <- function(i,
                              show_labels_v,
                              show_legend,
                              is_col_bottom,
                              label_expand_y_v = 3.0) {
    def      <- step_defs[[i]]
    src_step <- make_region_feats(g_all, region_col = def$col)

    # Step label slightly larger than genome-map default for readability
    sub_size <- 8.5

    v_panel <- build_track_panel(
      gene_v, src_step, s,
      scale_name         = V_SCORE_NAME,
      scale_option       = V_SCORE_OPT,
      scale_limits       = V_SCORE_LIMITS,
      subtitle           = def$label,
      left_label         = "V-score",
      show_seq           = FALSE,
      show_xaxis         = FALSE,
      show_labels        = show_labels_v,
      show_region_legend = show_legend,
      show_fill_legend   = show_legend,
      pad_for_labels     = show_labels_v,
      label_expand_y     = label_expand_y_v,
      repel_force        = if (show_labels_v) 1.5 else 0.4,
      repel_force_pull   = if (show_labels_v) 0.01 else 0.1,
      annot_size         = if (show_labels_v) 2.5 else 2.0,
      subtitle_size      = sub_size,
      # For the annotated step-1 panel: push labels all the way to the y ceiling
      # and nudge strongly upward so they spread across the full label zone
      ylim_upper         = if (show_labels_v) STEP1_EXPAND_Y * 3.0 else NULL,
      ylim_lower         = if (show_labels_v) STEP1_EXPAND_Y * 0.55 else NULL,
      nudge_y_val        = if (show_labels_v) 0.2 else 0.25
    )
    # Negative bottom margin on V panel pulls the WV row up against it
    v_panel <- v_panel + theme(plot.margin = margin(0, 2, -3, 2, "mm"))

    wv_panel <- build_track_panel(
      gene_wv, src_step, s,
      scale_name         = WIN_VL_SCORE_NAME,
      scale_option       = WIN_VL_SCORE_OPT,
      scale_limits       = WIN_VL_SCORE_LIMITS,
      left_label         = "Window-average\nVL-score",
      show_seq           = is_col_bottom,
      show_xaxis         = is_col_bottom,
      show_labels        = FALSE,
      show_region_legend = FALSE,
      show_fill_legend   = FALSE,
      pad_for_labels     = FALSE,
      subtitle_size      = sub_size
    )
    # Negative top margin on WV panel mirrors the pull from the V side
    wv_panel <- wv_panel + theme(plot.margin = margin(-3, 2, 0, 2, "mm"))

    list(v = v_panel, wv = wv_panel)
  }

  # Build left column: header + step pairs with spacers between steps
  left_header  <- ggplot() +
    annotate("text", x = 0.5, y = 0.5, label = "Steps 1\u20132",
             fontface = "bold", size = 3.5) +
    theme_void()

  left_panels  <- list()
  left_heights <- c(HEADER_H)

  for (i in LEFT_STEPS) {
    is_top    <- (i == min(LEFT_STEPS))
    is_bottom <- (i == max(LEFT_STEPS))

    if (!is_top) {
      left_panels  <- c(left_panels, list(plot_spacer()))
      left_heights <- c(left_heights, SPACER_H)
    }

    pair <- make_step_pair(
      i,
      show_labels_v    = is_top,
      show_legend      = is_top,
      is_col_bottom    = is_bottom,
      label_expand_y_v = if (is_top) STEP1_EXPAND_Y else 2.1
    )
    left_panels  <- c(left_panels, list(pair$v, pair$wv))
    left_heights <- c(left_heights,
                      if (is_top) STEP1_V_H else STEP_V_H,
                      STEP_WV_H)
  }

  # Build right column: header + step pairs with spacers between steps
  right_header  <- ggplot() +
    annotate("text", x = 0.5, y = 0.5, label = "Steps 3\u20135",
             fontface = "bold", size = 3.5) +
    theme_void()

  right_panels  <- list()
  right_heights <- c(HEADER_H)

  for (i in RIGHT_STEPS) {
    is_first  <- (i == min(RIGHT_STEPS))
    is_bottom <- (i == max(RIGHT_STEPS))

    if (!is_first) {
      right_panels  <- c(right_panels, list(plot_spacer()))
      right_heights <- c(right_heights, SPACER_H)
    }

    pair <- make_step_pair(
      i,
      show_labels_v = FALSE,
      show_legend   = FALSE,
      is_col_bottom = is_bottom
    )
    right_panels  <- c(right_panels, list(pair$v, pair$wv))
    right_heights <- c(right_heights, STEP_V_H, STEP_WV_H)
  }

  # Assemble columns then combine side-by-side
  left_col <- wrap_plots(
    c(list(left_header), left_panels),
    ncol = 1, heights = left_heights
  ) + plot_layout(guides = "collect")

  right_col <- wrap_plots(
    c(list(right_header), right_panels),
    ncol = 1, heights = right_heights
  ) + plot_layout(guides = "collect")

  (left_col | right_col) +
    plot_layout(guides = "collect") +
    plot_annotation(
      title    = "Algorithm step progression",
      subtitle = contig_id,
      theme = theme(
        plot.title    = element_text(size = 10, face = "bold", hjust = 0.5,
                                     margin = margin(b = 1)),
        plot.subtitle = element_text(size = 7.5, face = "plain", hjust = 0.5,
                                     color = "grey25", margin = margin(b = 4))
      )
    )
}
```

# Main genome map building function

- Three V-score tracks (KEGG, Pfam, PHROG), each with DB-specific gene
  labels
- The window-avg VL panels live only in the algorithm-progression plot

``` r
build_genome_map <- function(contig_id, seqs, genes,
                              title     = contig_id,
                              use_focus = TRUE) {
  if (is.na(contig_id) || length(contig_id) == 0) {
    return(ggplot() +
             annotate("text", x = 0.5, y = 0.5, label = paste("No candidate\n", title), size = 5, color = "grey50") +
             theme_void())
  }

  win   <- if (use_focus) get_focus_window(genes, contig_id) else NULL
  g_all <- genes %>% filter(contig == contig_id)
  if (!is.null(win))
    g_all <- filter(g_all, contig_pos_start >= win$start, contig_pos_end <= win$end)

  s         <- make_seqs_gg(seqs, contig_id, window = win)
  src_feats <- make_region_feats(g_all)

  track_defs <- list(
    list(db = "KEGG",  col = "KEGG_V-score",  name_col = "KEGG_Description"),
    list(db = "Pfam",  col = "Pfam_V-score",  name_col = "Pfam_Description"),
    list(db = "PHROG", col = "PHROG_V-score", name_col = "PHROG_Description")
  )
  n_tracks <- length(track_defs)

  panels <- lapply(seq_along(track_defs), function(i) {
    def      <- track_defs[[i]]
    gene_trk <- make_gene_track(g_all,
                                score_col = def$col,
                                name_col  = def$name_col)
    is_top    <- (i == 1L)
    is_bottom <- (i == n_tracks)
    build_track_panel(
      gene_trk, src_feats, s,
      scale_name = V_SCORE_NAME, scale_option = V_SCORE_OPT,
      scale_limits = V_SCORE_LIMITS,
      subtitle = def$db,
      subtitle_position = "left",
      show_seq = is_bottom, show_xaxis = is_bottom,
      show_labels = TRUE,
      show_region_legend = is_top
    )
  })

  wrap_plots(panels, ncol = 1) +
    plot_layout(guides = "collect") +
    plot_annotation(
      title    = title,
      subtitle = contig_id,
      theme = theme(
        plot.title    = element_text(size = 10,  face = "bold",  hjust = 0.5,
                                     margin = margin(b = 1)),
        plot.subtitle = element_text(size = 7.5, face = "plain", hjust = 0.5,
                                     color = "grey25",
                                     margin = margin(b = 4))
      )
    )
}
```

# Build genome maps

## Plot 1: Whole viral contig (geNomad dataset)

``` r
p1 = build_genome_map(
  cand1, seqs, genes,
  title = "Whole viral contig (geNomad dataset)"
)
save_plot_both(p1, "genome_p1", 12, 3)
print(p1)
```

![](strict_viral_regions_figures_files/figure-gfm/p1-1.png)<!-- -->

## Plot 2: Whole viral contig (proGenomes dataset)

``` r
p2 = build_genome_map(
  cand2, seqs, genes,
  title = "Whole viral contig (proGenomes dataset)"
)
save_plot_both(p2, "genome_p2", 12, 3)
print(p2)
```

![](strict_viral_regions_figures_files/figure-gfm/p2-1.png)<!-- -->

## Plot 3: Integrated provirus

``` r
p3 = build_genome_map(
  cand3, seqs, genes,
  title = "Integrated provirus"
)
save_plot_both(p3, "genome_p3", 12, 3)
print(p3)
```

![](strict_viral_regions_figures_files/figure-gfm/p3-1.png)<!-- -->

## Plot 4: Mixed viral/MGE/host composition

``` r
p4 = build_genome_map(
  cand4, seqs, genes,
  title = "Mixed viral/MGE/host composition"
)
save_plot_both(p4, "genome_p4", 12, 3)
print(p4)
```

![](strict_viral_regions_figures_files/figure-gfm/p4-1.png)<!-- -->

## Plot 5: Algorithm progression

``` r
p5 = build_step_plot(cand5, seqs, genes)
save_plot_both(p5, "genome_strict_progression", 12, 5.5)
print(p5)
```

![](strict_viral_regions_figures_files/figure-gfm/p5-1.png)<!-- -->

# Combined panel

``` r
p5_composite <- build_step_plot(cand5, seqs, genes)

plot_list <- list(p5_composite, p1, p2, p3)
plot_tags <- c("A", "B", "C", "D")
plot_h    <- c(2.0, 1, 1, 1)

# Strip each plot's legend
strip_legend <- function(p) p & theme(legend.position = "none")
panels_naked <- lapply(plot_list, strip_legend)

# Build a horizontal shared legend
extract_legend <- function(p) {
  parts <- cowplot::get_plot_component(p, "guide-box", return_all = TRUE)
  if (!is.list(parts)) parts <- list(parts)
  is_real <- vapply(parts, function(g) inherits(g, "gtable"), logical(1))
  if (any(is_real)) parts[[which(is_real)[1]]] else parts[[1]]
}

make_shared_legend <- function() {
  d_region <- tibble(x = c(1, 2),
                     region = factor(c("TRUE", "FALSE"),
                                     levels = c("TRUE", "FALSE")))
  d_v  <- tibble(x = c(1, 2), v  = V_SCORE_LIMITS)
  d_wv <- tibble(x = c(1, 2), wv = WIN_VL_SCORE_LIMITS)

  p_legend <- ggplot() +
    geom_tile(data = d_region, aes(x = x, y = 1, fill = region),
              width = 1, height = 0.5, alpha = 0.5) +
    scale_fill_manual(values = REGION_COLORS, name = "Strict viral region",
                      breaks = c("TRUE", "FALSE"), labels = c("True", "False"),
                      guide  = guide_legend(order = 1,
                                            override.aes = list(alpha = 0.85))) +
    ggnewscale::new_scale_fill() +
    geom_tile(data = d_v, aes(x = x, y = 2, fill = v),
              width = 1, height = 0.5) +
    scale_fill_viridis_c(name = V_SCORE_NAME, option = V_SCORE_OPT,
                         limits = V_SCORE_LIMITS, oob = squish, na.value = "grey85",
                         guide  = guide_colorbar(order = 2,
                                                  barwidth = 6, barheight = 0.5)) +
    ggnewscale::new_scale_fill() +
    geom_tile(data = d_wv, aes(x = x, y = 3, fill = wv),
              width = 1, height = 0.5) +
    scale_fill_viridis_c(name = WIN_VL_SCORE_NAME, option = WIN_VL_SCORE_OPT,
                         limits = WIN_VL_SCORE_LIMITS, oob = squish, na.value = "grey85",
                         guide  = guide_colorbar(order = 3,
                                                  barwidth = 6, barheight = 0.5)) +
    theme_void() +
    theme(legend.position  = "bottom",
          legend.box       = "horizontal",
          legend.direction = "horizontal",
          legend.box.just  = "center",
          legend.text      = element_text(size = 9),
          legend.title     = element_text(size = 10, face = "bold", hjust = 1),
          legend.spacing.x = unit(1, "cm")
         )
  extract_legend(p_legend)
}

shared_legend <- make_shared_legend()

# Stack the five panels with cowplot tags
stacked <- cowplot::plot_grid(
  plotlist       = panels_naked,
  ncol           = 1,
  rel_heights    = plot_h,
  labels         = plot_tags,
  label_size     = 16,
  label_fontface = "bold",
  align          = "v",
  axis           = "lr"
)

# Add the shared legend underneath
composite <- cowplot::plot_grid(
  stacked, shared_legend,
  ncol        = 1,
  rel_heights = c(1, 0.05)
)

save_plot_both(composite, "genome_maps_combined", 12, 12)
print(composite)
```

![](strict_viral_regions_figures_files/figure-gfm/composite-map-1.png)<!-- -->

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
    ##  [1] ggnewscale_0.5.2   gggenomes_1.1.3    svglite_2.1.2      RColorBrewer_1.1-3
    ##  [5] scales_1.4.0       patchwork_1.3.2    cowplot_1.1.3      lubridate_1.9.4   
    ##  [9] forcats_1.0.0      stringr_1.5.1      dplyr_1.1.4        purrr_1.0.4       
    ## [13] readr_2.1.5        tidyr_1.3.1        tibble_3.2.1       ggplot2_3.5.2     
    ## [17] tidyverse_2.0.0    arrow_13.0.0      
    ## 
    ## loaded via a namespace (and not attached):
    ##  [1] generics_0.1.3    stringi_1.8.7     hms_1.1.3         digest_0.6.37    
    ##  [5] magrittr_2.0.3    evaluate_1.0.3    grid_4.3.3        timechange_0.3.0 
    ##  [9] fastmap_1.2.0     jsonlite_2.0.0    ggrepel_0.9.8     viridisLite_0.4.2
    ## [13] textshaping_0.3.7 cli_3.6.4         rlang_1.2.0       bit64_4.6.0-1    
    ## [17] withr_3.0.2       yaml_2.3.10       tools_4.3.3       tzdb_0.5.0       
    ## [21] colorspace_2.1-1  assertthat_0.2.1  vctrs_0.6.5       R6_2.6.1         
    ## [25] lifecycle_1.0.4   bit_4.6.0         ragg_1.3.3        pkgconfig_2.0.3  
    ## [29] pillar_1.10.2     gtable_0.3.6      Rcpp_1.0.14       glue_1.8.0       
    ## [33] systemfonts_1.2.1 xfun_0.52         tidyselect_1.2.1  knitr_1.50       
    ## [37] farver_2.1.2      htmltools_0.5.8.1 labeling_0.4.3    rmarkdown_2.29   
    ## [41] compiler_4.3.3    S7_0.2.2

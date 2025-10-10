# Enriched Heatmaps with Metaplots (R)

Generate publication ready **side-by-side** heatmap columns for multiple groups of matrices.  
Each column shows a **single metaplot** (column means) above its stacked heatmaps.  
All heatmaps share **one global color scale**. A **single legend column** on the right shows the signal palette and metaplot line keys.

---

## Features
- Any number of groups; each group contains 1+ matrices.
- All matrices must have the **same number of columns**.
- One metaplot per group (lines = per-matrix column means).
- Global percentile-based color scaling with robust legend options:
  - **capped** legend (ticks between min and a chosen high quantile)
  - **quantile** legend (ticks at chosen quantiles)
- Optional clipping of extreme highs for stable palettes.
- Uses `cairo_pdf` for crisp output.

---

## Install

```r
install.packages("circlize")
if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")
BiocManager::install("ComplexHeatmap")
install.packages("magick")  # optional; used in some environments
```
---
# Usage

```r
suppressPackageStartupMessages({
  library(ComplexHeatmap)
  library(circlize)
  library(grid)
  library(magick)
})

g1_mats   <- list(ipa_matrix, polyA_matrix, utr3_matrix)
g1_titles <- c("IPA Sites", "Poly(A) Sites", "3' UTR Sites")

g2_mats   <- list(ipa_matrix, polyA_matrix, utr3_matrix)
g2_titles <- c("IPA Sites 2", "Poly(A) Sites 2", "3' UTR Sites 2")

doEnrichedHeatmaps(
  mats_groups   = list(g1_mats, g2_mats),
  titles_groups = list(g1_titles, g2_titles),
  group_names   = c("Condition A", "Condition B"),
  out_dir       = "out/",
  file_name     = "sample_output.pdf",
  palette       = c("white", "blue"),
  legend_mode   = "capped",          # or "quantile"
  legend_q_hi   = 0.90,              # for capped legend
  legend_quantiles = c(0,.1,.25,.5,.75,.9,.95), # for quantile legend
  clip_high_q   = 0.95,              # clamp extreme highs in palette computation
  width_per_group = 2.3,
  height_per_panel = 2,
  meta_height_mm = 25,
  meta_gap_mm = 1
)

```
---

# Function

```r
doEnrichedHeatmaps(
  mats_groups,                # list(list(M1,...), list(M1b,...), ...)
  titles_groups,              # list(c("t1","t2",...), c("t1b","t2b",...))
  out_dir,
  file_name = "two_groups_side_by_side.pdf",
  group_names = NULL,
  order_rows_by_rowmean = TRUE,
  palette = c("white","blue"),
  width_per_group = 2,
  height_per_panel = 2.0,
  meta_height_mm = 25,        # annotation height (label + metaplot)
  meta_gap_mm = 1,            # gap between metaplot and first heatmap
  show_signal_legend = TRUE,
  legend_mode = c("capped","quantile"),
  legend_q_hi = 0.95,
  legend_quantiles = c(0, .1, .25, .5, .75, .9, .95),
  clip_high_q = NA_real_      # e.g. 0.95 to clip outliers for palette calc
)

```
---

# Notes

**Constraint:** all matrices in all groups must share the same number of columns.

Rows can be auto-sorted by decreasing row mean per matrix.

Color scale is computed globally across all groups.

To avoid faint lines at the bottom of metaplots, the function uses `cairo_pdf` and padded ranges.


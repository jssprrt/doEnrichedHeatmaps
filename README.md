# Side-by-Side Heatmaps with Metaplots (R)

Generate publication-ready **side-by-side** heatmap columns for multiple groups of matrices.  
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

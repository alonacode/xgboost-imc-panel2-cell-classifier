# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## What this is

A supervised **XGBoost cell-type classifier for Panel 2 Imaging Mass Cytometry (IMC)** data. It is
an interactive R/knitr analysis workflow — *not* an application or installable package. It is built
on the Bioconductor IMC stack (`imcRtools`, `cytomapper`, `CATALYST`, `SingleCellExperiment` /
`SpatialExperiment`, `xgboost`, `caret`). The classifier is trained on manually gated cells from a
curated set of images and assigns biological cell types (T/B subsets, macrophage/DC populations,
neutrophil/MDSC, stroma, tumor, etc.) to every segmented cell.

## Running the workflow

There is no build, lint, or test system. The entire pipeline lives in
`XGBoost_classifier/code/01_read_data.Rmd` and is run **chunk-by-chunk in RStudio** (or via
`rmarkdown::render()`). Heavy compute (grid-search CV) is in the `train-model-cv` chunk.

Dependencies are *not* declared in any manifest — they are the `library()` calls at the top of the
Rmd. Install via Bioconductor/CRAN before running.

Upstream cell measurements are produced by the `steinbock` CLI before any R code runs (see the
`shell` chunk):

```
steinbock preprocess imc images --hpf 50
steinbock segment deepcell --minmax
steinbock measure intensities
steinbock measure regionprops
steinbock measure neighbors --type centroids --dmax 15
```

## Critical setup gotchas

Read these before trying to run anything — the workflow will **not** run from a fresh clone as-is:

- **Hardcoded `base_path`.** The Rmd sets `base_path` to `~/Desktop/XGBoost_classifier/Panel_2_10`,
  but this repo lives elsewhere (e.g. `~/Repos/genomic/...`). Point it at the actual
  `XGBoost_classifier/Panel_2_10` directory.
- **Output dir casing is inconsistent.** The code does `dir.create(.../"ROut")` and saves to `ROut`,
  while the tracked placeholder is `Rout/.gitkeep` and a commented line references `ROut_rds_all`.
  Don't assume a single canonical name.
- **Data is gitignored.** A clone contains only `panel.csv`, `images.csv`, `gates/*.rds`, and
  `.gitkeep` placeholders. The actual `img/`, `masks/`, `raw/`, `intensities/`, `regionprops/`,
  `neighbors/`, the `*_compensation/` dir, and `ROut/` outputs must be supplied separately.
- **The README's `02_`–`06_` Rmds do not exist.** The README describes a 6-stage pipeline, but all
  logic (read → compensate → label → split → train) is contained in the single `01_read_data.Rmd`.

## Architecture & data model

The central object is a `SpatialExperiment`/`SingleCellExperiment` ("SPE"): **columns = cells,
rows = IMC markers**. Key `colData` columns: `sample_id` (image), `ObjectNumber` (segmentation
label), and `cell_labels` (final cell type). The pipeline stages, all in `01_read_data.Rmd`:

1. **Read** — `read_steinbock()` builds the SPE from steinbock output.
2. **Spillover compensation** — `CATALYST` estimates a spillover matrix from the compensation TIFFs;
   `compCytof()` then **overwrites the `counts`/`exprs` assays in place** with compensated values.
   The `exprs` assay is `asinh(counts / cofactor)`-transformed.
3. **Label import** — manual gate `.rds` files in `gates/` are read and matched to cells by
   `cell_id = paste0(sample_id, "_", ObjectNumber)`. Raw gate label strings are normalized
   (trailing `_NN` suffixes stripped, `Cd4`→`CD4`, etc.).
4. **Label consolidation** — same-label duplicates collapsed; `Tumor + other` doublets resolved to
   the non-Tumor label; remaining true multi-label doublets dropped. One blank image is excluded.
5. **Biological cleanup rules** (before training): drop `Igg` cells with `IgG < 1`, drop `NK` cells
   with `CD7 < 0.5`, relabel `mregDC` with high epithelial signal (`panCK > 0.5` or `Ecad > 0.8`)
   to `Tumor`.
6. **Train/test split** — at the **image level** (all cells of an image go to one side). Uses a
   custom `split_test_rows()` optimizer that targets per-cell-type test fractions (up-weighting rare
   types like MDSC/BnT/Igg/mregDC). `caret::groupKFold` on `sample_id` is also used so cells from one
   image never leak across CV folds.
7. **Feature matrix** — DNA/Histone channels excluded via `!grepl("DNA|Histone", rownames(...))`;
   the assay is **transposed** (`t(assay(spe, "exprs")[marker_keep, ])`) because XGBoost expects
   cells-as-rows. Labels become 0-based integers; `label_mapping <- levels(factor(cell_labels))`
   defines the numeric↔name mapping. Inverse-frequency class weights are applied, with `undefined`
   and `Igg` manually down-weighted.
8. **Tuning & final model** — `xgb.cv` grid search (`multi:softprob` / `mlogloss`) over the
   image-grouped folds selects the best params and `nrounds`; `xgb.train` then retrains on all
   training cells. The 30-image held-out set is used for rare-class evaluation.

## Model artifacts contract

When the model is exported (described in README; output files not yet committed), predictions must:
- use markers in the exact order stored in `feature_names`, and
- decode predicted numeric classes via the saved `label_mapping`.

Prefer the portable `classifier_xgboost_final.json` + `classifier_xgboost_final_meta.rds` over the
R-native `.rds` model for downstream prediction.

## Directory layout

```
XGBoost_classifier/
  README.md                      # workflow narrative & modeling notes
  code/01_read_data.Rmd          # the entire pipeline
  Panel_2_10/                    # steinbock output layout for one panel/batch
    panel.csv                    # 44-marker panel (channel, name, keep flags) — tracked
    images.csv                   # per-image acquisition metadata — tracked
    gates/                       # manual gate label .rds files — tracked
    img/ masks/ raw/ intensities/ regionprops/ neighbors/   # gitignored data
    *_compensation/              # spillover TIFFs + panel/images csv — gitignored
    ROut/ (a.k.a. Rout)          # saved SPE/model outputs — gitignored
```

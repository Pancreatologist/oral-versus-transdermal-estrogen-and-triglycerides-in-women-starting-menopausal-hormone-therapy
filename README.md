# Visualization code

This repository contains the R code used to produce the figures for the article:

> **Oral versus transdermal estrogen and triglycerides in women starting menopausal hormone therapy: a target trial emulation study**

All inputs to this code are **simulated data** generated inside the script from
a fixed random seed. The code demonstrates the complete figure *design* of the
article — participant flow, forest plots, covariate-balance diagnostics, and
multiple-imputation diagnostics — but it does not contain, read, or reproduce
any real data or study results.

## Figures generated

| File | Content |
| --- | --- |
| `figure_participant_flow` | Participant flow diagram |
| `figure_forest_primary` | Primary forest plot (grouped analyses) |
| `figure_forest_full` | Extended forest plot (additional sensitivity analyses) |
| `figure_tg_baseline_followup` | Median (IQR) triglycerides, HDL-cholesterol, and LDL-cholesterol at baseline and follow-up by treatment arm |
| `figure_balance_diagnostics` | Propensity-score overlap, Love plot, weight distribution |
| `figure_mice_diagnostics` | Multiple-imputation trace plots and observed vs imputed densities |
| `figure_supp_diagnostics` | Combined diagnostics supplement (panels A–E) |
| `figure_supp_additional_results` | Combined results supplement (panels A–D) |

Each figure is written as a vector PDF and a 600–1000 dpi LZW-compressed TIFF
into the `figures/` output directory.

## Requirements

- R (>= 4.1)
- R packages from CRAN:

```r
install.packages(c("data.table", "ggplot2", "patchwork", "ragg"))
```

The Arial font family is expected for publication-style rendering; substitute
another font in the script if Arial is unavailable on your system.

## Usage

From the repository root:

```sh
Rscript figures.R
```

The script is fully self-contained and deterministic (fixed seed): it
simulates the analytic inputs, builds every figure, and writes the output
files to `figures/`. Total run time is well under a minute.

## Design notes

- Sizes follow journal double-column conventions (90 / 140 / 190 mm widths).
- Treatment arms use a fixed colour-blind-safe palette: oral = blue,
  transdermal = orange; points additionally differ in shape.
- Forest plots use a log-ratio x axis labelled in percentage points, with a
  dashed reference line at no difference and direction annotations.
- Small-cell suppression is applied wherever counts are displayed, as a
  generic disclosure-avoidance measure.

## License

This code is released under the MIT License — see [LICENSE](LICENSE).

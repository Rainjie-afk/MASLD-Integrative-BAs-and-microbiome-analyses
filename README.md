# MASLD-Integrative-BAs-and-microbiome-analyses

Reproducible R helpers for integrative bile-acid and microbiome analyses in
NAFLD/MASLD. The implementation is in `R/integrative_analysis.R`; focused
checks are in `tests/test_integrative_analysis.R`.

## Workflow

### Bile acids

Provide one matrix per compartment (liver, plasma, ileum, cecum, and feces),
with samples in rows and bile acids in columns, plus a common group vector.
`run_diablo()` fits a mixOmics DIABLO model and uses a design matrix that
connects all compartments. The resulting model can be used for selected-feature
plots, correlation/Circos visualizations, and group discrimination. Use
`run_piecewise_sem()` with explicitly specified directional component models
when evaluating the hypothesized liver → plasma → intestinal relationships.

The repository intentionally leaves plotting and SEM formulas explicit in the
analysis notebook/script: those choices depend on the study's measurement
units, covariates, and biological hypotheses rather than being silently
invented by a helper function.

### Microbiome machine learning

`relative_abundance()` normalizes each sample independently. For every outer
cross-validation split, call `prepare_fold_data()` to:

1. fit prevalence/abundance filtering on the training samples only;
2. apply that fixed feature set to the held-out samples (missing taxa become
   zero); and
3. apply the arcsine square-root transform to relative abundance proportions.

`repeated_outer_cv()` supports repeated stratified outer folds. Pass a
parameter-grid data frame to run `inner_cv_tune()` inside each outer training
split; both loops refit filtering using training samples only. The supplied
`fit_model(x, y, parameters)` and `predict_model(model, x)` callbacks keep the
classifier choice explicit. Evaluate out-of-fold predictions with
`multiclass_metrics()` (accuracy, balanced accuracy, and macro-F1). This keeps
preprocessing out of the held-out data and avoids information leakage.

The supplied helpers are model-agnostic: use the transformed matrices with
the selected classifier and tune it inside the inner loop. A two-hidden-layer
MLP/DNN may be added as an exploratory classifier, but the dataset size
(`n = 24`) is too small for confident neural-network inference; report it as
exploratory and retain the repeated out-of-fold estimates.

## Running the checks

With R installed, run:

```sh
Rscript tests/test_integrative_analysis.R
```

`mixOmics` and `piecewiseSEM` are optional dependencies and are required only
when their corresponding model helpers are called.

source(file.path("R", "integrative_analysis.R"))

counts <- matrix(c(1, 3, 0, 1, 1, 2, 4, 0, 2, 2, 1, 1),
                 nrow = 4, byrow = TRUE,
                 dimnames = list(paste0("sample", 1:4), paste0("taxon", 1:3)))

ra <- relative_abundance(counts)
stopifnot(all(abs(rowSums(ra) - 100) < 1e-10))

fold_data <- prepare_fold_data(counts, train = 1:3, test = 4,
                               min_prevalence = 0.5, min_abundance = 20)
stopifnot(identical(colnames(fold_data$train), colnames(fold_data$test)))
stopifnot(identical(fold_data$filter$features, colnames(fold_data$train)))
stopifnot(all(fold_data$train >= 0 & fold_data$train <= pi / 2))

folds <- stratified_folds(c("A", "A", "B", "B", "C", "C"), folds = 2,
                          repeats = 2, seed = 42)
stopifnot(length(folds) == 4)
stopifnot(all(vapply(folds, function(x) {
  identical(sort(c(x$train, x$test)), 1:6) && !anyDuplicated(c(x$train, x$test))
}, logical(1))))

metrics <- multiclass_metrics(c("A", "B", "C"), c("A", "C", "C"))
stopifnot(all(names(metrics) %in% c("accuracy", "balanced_accuracy", "macro_f1")))

correlations <- cross_compartment_correlations(
  list(liver = counts[, 1:2, drop = FALSE],
       plasma = counts[, 2:3, drop = FALSE]))
stopifnot(nrow(correlations) == 4)
stopifnot(nrow(circos_association_data(correlations, threshold = 0)) == 4)

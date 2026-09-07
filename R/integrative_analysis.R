# Reusable functions for the bile-acid and microbiome analyses.

relative_abundance <- function(counts) {
  counts <- as.matrix(counts)
  if (!is.numeric(counts) || any(!is.finite(counts)) || any(counts < 0)) {
    stop("counts must be a finite, non-negative numeric matrix")
  }
  totals <- rowSums(counts)
  if (any(totals == 0)) stop("samples with zero total counts cannot be normalized")
  sweep(counts, 1, totals, "/") * 100
}

fit_microbiome_filter <- function(train_counts, min_prevalence = 0.10,
                                  min_abundance = 0.01) {
  if (!is.numeric(min_prevalence) || min_prevalence < 0 || min_prevalence > 1) {
    stop("min_prevalence must be between 0 and 1")
  }
  if (!is.numeric(min_abundance) || min_abundance < 0) {
    stop("min_abundance must be non-negative")
  }
  train_ra <- relative_abundance(train_counts)
  keep <- colMeans(train_ra >= min_abundance) >= min_prevalence
  kept_features <- colnames(train_ra)[keep]
  if (is.null(kept_features)) kept_features <- character()
  list(features = kept_features, min_abundance = min_abundance)
}

apply_microbiome_filter <- function(counts, filter) {
  ra <- relative_abundance(counts)
  missing <- setdiff(filter$features, colnames(ra))
  if (length(missing)) {
    ra <- cbind(ra, matrix(0, nrow(ra), length(missing),
                            dimnames = list(rownames(ra), missing)))
  }
  ra[, filter$features, drop = FALSE]
}

arcsine_sqrt_transform <- function(relative_abundance_percent) {
  x <- as.matrix(relative_abundance_percent) / 100
  if (any(!is.finite(x)) || any(x < 0 | x > 1)) {
    stop("relative abundances must be percentages in [0, 100]")
  }
  asin(sqrt(x))
}

stratified_folds <- function(labels, folds = 5, repeats = 1, seed = 1) {
  if (folds < 2 || repeats < 1) stop("folds must be >= 2 and repeats must be >= 1")
  labels <- as.factor(labels)
  set.seed(seed)
  result <- vector("list", folds * repeats)
  k <- 1
  for (repeat_id in seq_len(repeats)) {
    fold_id <- integer(length(labels))
    for (group in levels(labels)) {
      indices <- sample(which(labels == group))
      fold_id[indices] <- rep(seq_len(folds), length.out = length(indices))
    }
    for (fold in seq_len(folds)) {
      result[[k]] <- list(repeat_id = repeat_id, fold = fold,
                          train = which(fold_id != fold),
                          test = which(fold_id == fold))
      k <- k + 1
    }
  }
  result
}

macro_f1 <- function(actual, predicted) {
  classes <- union(levels(as.factor(actual)), levels(as.factor(predicted)))
  scores <- vapply(classes, function(class) {
    tp <- sum(actual == class & predicted == class)
    fp <- sum(actual != class & predicted == class)
    fn <- sum(actual == class & predicted != class)
    if (2 * tp + fp + fn == 0) 0 else 2 * tp / (2 * tp + fp + fn)
  }, numeric(1))
  mean(scores)
}

multiclass_metrics <- function(actual, predicted) {
  actual <- as.factor(actual)
  predicted <- factor(predicted, levels = levels(actual))
  recalls <- vapply(levels(actual), function(class) {
    denominator <- sum(actual == class)
    if (denominator == 0) NA_real_ else sum(actual == class & predicted == class) / denominator
  }, numeric(1))
  c(accuracy = mean(actual == predicted, na.rm = TRUE),
    balanced_accuracy = mean(recalls, na.rm = TRUE),
    macro_f1 = macro_f1(actual, predicted))
}

# This helper deliberately receives a training-only filter. Callers should fit
# it inside every outer and inner split, never once on the complete data set.
inner_cv_tune <- function(counts, labels, fit_model, predict_model,
                          parameter_grid, folds = 3, seed = 1,
                          min_prevalence = 0.10, min_abundance = 0.01) {
  splits <- stratified_folds(labels, folds, repeats = 1, seed)
  scores <- vapply(seq_len(nrow(parameter_grid)), function(row) {
    fold_scores <- vapply(splits, function(split) {
      prepared <- prepare_fold_data(counts, split$train, split$test,
                                    min_prevalence, min_abundance)
      parameters <- as.list(parameter_grid[row, , drop = FALSE])
      model <- fit_model(prepared$train, labels[split$train], parameters)
      predictions <- predict_model(model, prepared$test)
      multiclass_metrics(labels[split$test], predictions)[["macro_f1"]]
    }, numeric(1))
    mean(fold_scores)
  }, numeric(1))
  best <- which.max(scores)
  list(parameters = as.list(parameter_grid[best, , drop = FALSE]),
       scores = scores, parameter_grid = parameter_grid)
}

repeated_outer_cv <- function(counts, labels, fit_model, predict_model,
                              folds = 5, repeats = 5, seed = 1,
                              min_prevalence = 0.10, min_abundance = 0.01,
                              parameter_grid = NULL, inner_folds = 3) {
  if (length(labels) != nrow(counts)) stop("labels must match count rows")
  splits <- stratified_folds(labels, folds, repeats, seed)
  predictions <- vector("list", length(splits))
  for (i in seq_along(splits)) {
    split <- splits[[i]]
    prepared <- prepare_fold_data(counts, split$train, split$test,
                                  min_prevalence, min_abundance)
    parameters <- list()
    if (!is.null(parameter_grid)) {
      tuning <- inner_cv_tune(counts[split$train, , drop = FALSE],
                              labels[split$train], fit_model, predict_model,
                              parameter_grid, inner_folds, seed + i,
                              min_prevalence, min_abundance)
      parameters <- tuning$parameters
    }
    model <- fit_model(prepared$train, labels[split$train], parameters)
    predictions[[i]] <- data.frame(
      repeat_id = split$repeat_id, fold = split$fold,
      sample = split$test, truth = labels[split$test],
      prediction = predict_model(model, prepared$test),
      stringsAsFactors = FALSE
    )
  }
  out_of_fold <- do.call(rbind, predictions)
  list(predictions = out_of_fold,
       metrics = multiclass_metrics(out_of_fold$truth, out_of_fold$prediction))
}

prepare_fold_data <- function(counts, train, test, min_prevalence = 0.10,
                              min_abundance = 0.01) {
  filter <- fit_microbiome_filter(counts[train, , drop = FALSE],
                                  min_prevalence, min_abundance)
  list(train = arcsine_sqrt_transform(
    apply_microbiome_filter(counts[train, , drop = FALSE], filter)),
       test = arcsine_sqrt_transform(
         apply_microbiome_filter(counts[test, , drop = FALSE], filter)),
       filter = filter)
}

run_diablo <- function(blocks, outcome, ncomp = 2, keepX = NULL) {
  if (!requireNamespace("mixOmics", quietly = TRUE)) {
    stop("The mixOmics package is required for DIABLO")
  }
  if (is.null(names(blocks)) || any(!nzchar(names(blocks)))) {
    stop("blocks must be a named list of compartment matrices")
  }
  design <- matrix(0.1, length(blocks), length(blocks),
                   dimnames = list(names(blocks), names(blocks)))
  diag(design) <- 0
  if (is.null(keepX)) keepX <- lapply(blocks, function(x) rep(min(10, ncol(x)), ncomp))
  model <- mixOmics::block.splsda(X = blocks, Y = as.factor(outcome),
                                  ncomp = ncomp, keepX = keepX, design = design)
  model
}

run_splsda <- function(x, outcome, ncomp = 2, keepX = NULL) {
  if (!requireNamespace("mixOmics", quietly = TRUE)) {
    stop("The mixOmics package is required for sPLS-DA")
  }
  if (is.null(keepX)) keepX <- rep(min(10, ncol(x)), ncomp)
  mixOmics::splsda(x, as.factor(outcome), ncomp = ncomp, keepX = keepX)
}

cross_compartment_correlations <- function(blocks, selected = NULL,
                                           method = "spearman") {
  if (is.null(names(blocks)) || length(blocks) < 2) {
    stop("blocks must contain at least two named compartments")
  }
  if (!method %in% c("pearson", "spearman", "kendall")) stop("unsupported method")
  pairs <- combn(names(blocks), 2, simplify = FALSE)
  do.call(rbind, lapply(pairs, function(pair) {
    left <- blocks[[pair[1]]]
    right <- blocks[[pair[2]]]
    if (!is.null(selected)) {
      left <- left[, intersect(colnames(left), selected), drop = FALSE]
      right <- right[, intersect(colnames(right), selected), drop = FALSE]
    }
    if (!ncol(left) || !ncol(right)) return(NULL)
    values <- cor(left, right, method = method, use = "pairwise.complete.obs")
    indices <- expand.grid(left = seq_len(ncol(left)), right = seq_len(ncol(right)))
    data.frame(compartment_1 = pair[1],
               feature_1 = colnames(left)[indices$left],
               compartment_2 = pair[2],
               feature_2 = colnames(right)[indices$right],
               correlation = as.vector(values), stringsAsFactors = FALSE)
  }))
}

circos_association_data <- function(correlations, threshold = 0.5) {
  required <- c("compartment_1", "feature_1", "compartment_2",
                "feature_2", "correlation")
  if (!all(required %in% names(correlations))) {
    stop("correlations is missing required columns")
  }
  correlations[abs(correlations$correlation) >= threshold, required]
}

fit_exploratory_dnn <- function(x, y, hidden = c(16, 8), epochs = 100,
                                seed = 1) {
  if (!requireNamespace("keras3", quietly = TRUE)) {
    stop("The keras3 package is required for the exploratory DNN")
  }
  if (length(hidden) != 2 || any(hidden < 1)) stop("hidden must have two positive sizes")
  set.seed(seed)
  model <- keras3::keras_model_sequential() |>
    keras3::layer_dense(units = hidden[1], activation = "relu",
                        input_shape = ncol(x)) |>
    keras3::layer_dense(units = hidden[2], activation = "relu") |>
    keras3::layer_dense(units = length(unique(y)), activation = "softmax")
  model |> keras3::compile(optimizer = "adam",
                           loss = "sparse_categorical_crossentropy",
                           metrics = "accuracy")
  model |> keras3::fit(x, as.integer(as.factor(y)) - 1, epochs = epochs,
                        verbose = 0)
  model
}

run_piecewise_sem <- function(data, models) {
  if (!requireNamespace("piecewiseSEM", quietly = TRUE)) {
    stop("The piecewiseSEM package is required for piecewise SEM")
  }
  if (!is.list(models) || !length(models)) stop("models must be a non-empty list")
  piecewiseSEM::psem(models, data = data)
}

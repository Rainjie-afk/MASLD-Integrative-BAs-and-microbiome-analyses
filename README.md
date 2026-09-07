<<<<<<< HEAD
# Integrative Bile Acid and Microbiome Analyses in MASLD

This repository provides a workflow for **integrative multi-compartment bile acid and microbiome analyses in MASLD**, combining bile acid profiles from the **liver, plasma, ileum, cecum, and feces**.

The workflow uses multivariate methods implemented in the **[mixOmics](https://mixomics.org/)** framework together with complementary statistical approaches to characterize coordinated alterations in bile acid metabolism across anatomical compartments and investigate potential inter-compartment relationships in MASLD.

## Bile Acid Analysis

The workflow includes:

* **DIABLO** for multiblock integration of bile acid profiles across tissues and biological compartments.
* **Circos plots** for visualization of cross-compartment associations among selected bile acids.
* **Correlation analysis** for identifying relationships among bile acids across compartments.
* **Piecewise structural equation modeling (piecewise SEM)** for evaluating potential directional relationships among bile acid profiles in the liver, plasma, ileum, cecum, and feces.
* **Optional sPLS-DA** for supervised group discrimination and feature selection.

Together, these analyses provide complementary approaches for identifying bile acids that contribute to coordinated metabolic changes across multiple anatomical compartments.

## Microbiome Machine-Learning Analysis

An optional machine-learning pipeline is included for **three-class MASLD microbiome classification**.

The analysis follows the workflow below:

1. Raw microbial counts are converted to **relative abundance (%)** independently for each sample.

2. Relative-abundance filtering is determined **exclusively from the training data within each cross-validation fold** to minimize information leakage.

3. Relative abundances are transformed using the **arcsine square-root transformation**:


4. **Repeated stratified outer cross-validation** is used to generate out-of-fold predictions and estimate model performance.

5. **Hyperparameter tuning** is performed within an **inner stratified cross-validation loop**.

6. Model performance is evaluated using **multiclass metrics** appropriate for the three-group outcome.

7. A **deep neural network (DNN)** is included as an exploratory model using a two-hidden-layer multilayer perceptron (MLP). Because the dataset contains only **n = 24 samples**, results from the DNN should be interpreted cautiously due to the limited sample size available for neural-network training.

## Overview

The overall workflow is designed to integrate complementary information from bile acid profiling and microbiome analyses in order to:

* identify coordinated bile acid alterations across anatomical compartments;
* characterize cross-compartment associations;
* evaluate potential directional relationships among bile acid profiles;
* identify features contributing to MASLD group discrimination; and
* explore the predictive potential of microbiome profiles using machine-learning approaches.

The emphasis is on combining **multivariate integration, statistical modeling, visualization, and predictive analysis** to provide a systems-level view of bile acid–microbiome relationships in MASLD.
=======
# MASLD-Integrative-BAs-and-microbiome-analyses
Integrative bile acids and microbiome analyses
>>>>>>> af3a3aeca85041c5a41e659e9f359eb2705641a7

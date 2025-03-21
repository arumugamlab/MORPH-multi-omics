# MORPH-multi-omics
Multi-Omics for Relationships in Phenotypes and Health (MORPH)

# XGBoost AUC Function: 

## Overview

This R function implements an XGBoost-based classification pipeline, including feature reduction and performance evaluation via cross-validation. It supports multiple omics datasets and outputs a range of performance metrics and model details.

## Inputs

| Parameter        | Type              | Description |
|-----------------|------------------|-------------|
| `omics`         | Character Vector  | Types of omics data (e.g., "Metabolomics", "Metagenomics", "Proteomics"). |
| `feature`       | DataFrame         | Target labels for classification. |
| `seed`          | Integer           | Random seed for reproducibility. |
| `top`           | Integer           | Number of top features to retain in feature reduction. |
| `df.matrix`     | DataFrame         | Data matrix with omics features (columns) and samples (rows). **Note:** each feature must be prefixed with its corresponding omic type (e.g.,"Proteomics.P10643")|
| `panel`         | String            | Column name on the metadata of the feature that is going to be analyzed (e.g., "PreACLF"). |
| `nfolds`        | Integer           | Number of folds for cross-validation. |
| `fecal.features` | Character Vector  | List of microbiome-related features. |
| `plasma.features` | Character Vector  | List of host-related features. |
| `metadata`      | DataFrame         | Metadata containing sample IDs and feature of interest. |

## Outputs

The function returns a list of outputs covering model performance, predictions, and feature importance:

| Output                         | Type            | Description |
|--------------------------------|----------------|-------------|
| `loss.list`                    | List (Plots)   | Loss curves for each omic type, showing training vs validation loss over one fold. |
| `cvAUC.list`                   | List           | Cross-validated AUC scores for each omic type. |
| `cvAUC.list.top`               | List           | Cross-validated AUC scores using only the top X most important features. |
| `feature.importance.list`       | List (DataFrames) | Feature importance for each fold and omic type from the full dataset model. |
| `feature.importance.list.top`   | List (DataFrames) | Feature importance using only the reduced top features. |
| `cvAUC.value`                  | List (Numeric) | Final AUC values after cross-validation for each omic type. |
| `cvAUC.se`                     | List (Numeric) | Standard error of AUC across cross-validation folds. |
| `model_list.all`                | List (Models)  | Trained XGBoost models for each fold and omic type. |
| `model_list.all.top`            | List (Models)  | XGBoost models trained using only the top features. |
| `predictions.list`              | List (Vectors) | Predicted probabilities for each fold’s validation set. |
| `predictions.list.top`          | List (Vectors) | Predicted probabilities from models trained on top features only. |
| `labelsValidation.list`         | List (Vectors) | True validation set labels for each fold. |
| `folds`                         | List (Vectors) | Indices of training/validation samples for each cross-validation fold. |
| `out`                           | List           | AUC and confidence intervals for the complete model. |
| `out.top`                       | List           | AUC and confidence intervals for the top feature model. |
| `score_list.all`                | List (DataFrames) | AUC scores, rounds, and evaluation logs per fold. |
| `confusionMatrix.list.all`      | List (Matrices) | Confusion matrices for each fold’s predictions. |
| `confusionMatrix.lit.top.all`   | List (Matrices) | Confusion matrices from top feature models. |
| `train.labels.ID.out`           | List (Vectors) | Training set sample IDs per fold. |
| `val.labels.ID.out`             | List (Vectors) | Validation set sample IDs per fold. |

## Example Usage

```r
result <- xgboost_AUC(
  omics =  c("Proteomics","Metabolomics","Microbiome", "Host"),
  feature = feature_df,
  seed = 123,
  top = 15,
  df.matrix = df_matrix,
  panel = "PreACLF",
  nfolds = 10,
  fecal.features = fecal_features,
  plasma.features = plasma_features,
  metadata = metadata_df
)

# Access results:
result$cvAUC.list          # Full model AUCs
result$cvAUC.list.top      # Top feature model AUCs
result$feature.importance.list # Feature importance for full models
result$feature.importance.list.top # Feature importance for top models
```

## Dependencies

Ensure these packages are installed:

```r
install.packages(c("dplyr", "xgboost", "caret", "pROC", "Metrics", "reshape2", "ggplot2", "xpectr"))
```

# XGBoost summary Function

## Overview

This R function  `xgboost.plots`, = generates performance plots for a trained XGBoost model and outputs relevant performance metrics and visualizations.

## Inputs

| Parameter       | Type                  | Description |
|-----------------|-----------------------|-------------|
| `xgboost.object` | Output of XGBoost_AUC | A trained XGBoost model with AUC values and feature importance information. |
| `top`           | Integer               | Number of top features to include in visualizations. |

## Outputs

| Output                          | Type            | Description |
|---------------------------------|-----------------|-------------|
| `plot.complete`                 | List (Plots)    | Complete set of plots combining AUC curves and feature upset plots. |
| `AUC.top.feature.curve`         | Plot            | Line plot showing AUC scores for top features across models. |
| `plot.AUC.all.feature.reduced`  | List (Plots)    | AUC bar plots for selected top features. |
| `AUC.barplot.all.omics`         | Plot            | Bar plot summarizing AUC values for all omics data. |
| `confusionMatrix.best.model`    | List (Plots)    | Confusion matrix for the best-performing model across feature reduced models. |
| `bestModel`                     | DataFrame       | Summary table of the best-performing models per omic. |
| `topFeatures.list`              | List (Vectors)  | List of top features from the best models. |
| `AUC.barplot.omics.all.bestmodels` | Plot         | AUC bar plot for best models across all omics. |
| `bestModel.top`                 | DataFrame       | Table summarizing the minimum number of features in the best-performing models. |
| `AUC.all.top`                   | DataFrame       | Melted data frame of AUC values across top models and features. |
| `cvAUC.list.melt`               | DataFrame       | Refined data frame of AUC values with features reordered for visualization. |

## Summary Code Breakdown

The function performs the following steps:

1. **Generate AUC Bar Plots:**
   - Calculates cross-validation AUC values and their standard deviations.
   - Creates bar plots for each omic, displaying AUC scores and error bars.

2. **Plot Top Features:**
   - Extracts performance metrics for models with varying top features.
   - Generates AUC plots with error bars and labels showing top feature performance.

3. **Best Model Identification:**
   - Identifies the best model per omic based on maximum AUC and minimun number of features.
   - Extracts feature importance tables for each best-performing model.

4. **Confusion Matrix Visualization:**
   - Creates confusion matrix plots with color-coded correct/incorrect predictions.

5. **Upset Plot Generation:**
   - Uses the `ComplexUpset` package to generate upset plots for feature overlaps across feature reduced models.



## Example Usage

```r
# Load trained XGBoost model
result <- xgboost_AUC(
  omics = c("Proteomics","Metabolomics","Microbiome", "Host"),
  feature = feature_df,
  seed = 123,
  top = 15,
  df.matrix = df_matrix,
  panel = "PreACLF",
  nfolds = 10,
  fecal.features = fecal_features,
  plasma.features = plasma_features,
  metadata = metadata_df
)

# Generate performance plots for top 10 features
Result.plots <- xgboost.plots(xgboost.object = result, top = 10)

# Access AUC plot
result$AUC.top.feature.curve
```

## Dependencies

Ensure these packages are installed:

```r
install.packages(c("ggplot2", "reshape2", "dplyr", "ComplexUpset", "stringr", "ggrepel"))
```




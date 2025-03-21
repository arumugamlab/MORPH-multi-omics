
library(dplyr)
library(xgboost)
library(caret)
library(pROC)
library(Metrics)
library(reshape2)
library(ggplot2)
library(groupdata2)
library(xpectr)



xgboost_AUC <- function(omics, feature, seed = seed, top = top, df.matrix = df.matrix, panel=panel, nfolds = nfolds, fecal.features = fecal.features, plasma.features = plasma.features, metadata = metadata) {
  
  
  # Function to generate performance plots for an XGBoost model
  # INPUT:
  #     - omics: (Character Vector)	Types of omics data (e.g., "Microbiome", "Host", "All.Clinics").
  #.    - seed:	(Integer)	Random seed for reproducibility.
  #     - top	(Integer)	Number of top features to retain in feature reduction.
  #     - df.matrix: (DataFrame) Data matrix with omics features (columns) and samples (rows).
  #     - panel: (String)	Name of feature to of intered on the medatada	
  #     - nfolds:(Integer)	Number of folds for cross-validation.
  #     - fecal.features: (Character Vector) List of microbiome-related features.
  #     -  plasma.features: (Character Vector) 	List of host-related features.
  #     -  metadata: (DataFrame)	Metadata containing sample IDs and clinical details.
  
  
  
  # Check if essential inputs exist
  if (missing(df.matrix) || missing(metadata) || missing(feature) || missing(nfolds) || missing(fecal.features)
      || missing(plasma.features) || missing(panel)) {
    stop("Missing required input")
  }

  
  # Initialize lists to store global results
  #------------------------------------------------
  
  loss.list.all <- list()
  cvAUC.list.all <- list()
  cvAUC.list.top.all<- list()
  
  cvAUC.value <- list()
  cvAUC.se <- list()
  score_list.all<-list()
  predictions.list.top <- list()
  labelsValidation.list <- list()
  folds.list <- list()
  predictions.list.all <- list()
  predictions.list.all.top<-list()
  model_list.all <- list()
  model_list.all.top <- list()
  out.list <- list()
  out.list.top <- list()
  
  feature.importance.list.top.all <- list()
  feature.importance.list.all <- list()

  confusionMatrix.list.all <- list()
  confusionMatrix.lit.top.all <- list()
  
  train.labels.list.out <- list()
  val.labels.list.out <- list()
  
  loss.list <- list()
  cvAUC.list <- list()
  cvAUC.list.top <- list()
  
  #------------------------------------------------
  
  
  # Loop through each omic type
  for (j in seq_along(omics)) {
  
    # Initialize lists to store results
    #------------------------------------------------
    score_list <- list()
    predictions.list <- list()
    predictions.list.top <- list()
    model_list <- list()
    model_list.top <- list()
    feature.importance.list <- list()
    feature.importance.list.top <- list()
    confusionMatrix.list <- list()
    confusionMatrix.list.top <- list()
    train.labels.list <- list()
    val.labels.list <- list()
   
    #------------------------------------------------
    
    omi <- omics[j]

    cat("Processing Omic: ", omi, "\n")
    
    
    # Subset features based on omic type
    if (!omi %in% c("all", "Microbiome", "Host", "All.Clinics")) {
      omics.features <- colnames(df.matrix)[startsWith(colnames(df.matrix), omi)]
      df.matrix.omics <- df.matrix %>% dplyr::select(all_of(omics.features))
    } else if (omi == "all") {
      omics.features.fecal <- which(colnames(df.matrix) %in% fecal.features)
      omics.features.host <- which(colnames(df.matrix) %in% plasma.features)
      df.matrix.omics <- df.matrix[c(omics.features.fecal, omics.features.host)]
    } else if (omi == "Microbiome") {
      omics.features <- which(colnames(df.matrix) %in% fecal.features)
      df.matrix.omics <- df.matrix[omics.features]
    } else if (omi == "Host") {
      omics.features <- which(colnames(df.matrix) %in% plasma.features)
      df.matrix.omics <- df.matrix[omics.features]
    } else if (omi == "All.Clinics") {
      df.matrix.omics <- df.matrix
    }
    
    # Remove rows with all NA values
    df.matrix.omics <- df.matrix.omics[rowSums(is.na(df.matrix.omics)) < ncol(df.matrix.omics), ]
    

    # Add target label to the data
    df.matrix.omics <- merge(feature, df.matrix.omics, by = 0)
    row.names(df.matrix.omics) <- df.matrix.omics$Row.names
    df.matrix.omics$Row.names <- NULL
    
    # Create training data and training labels and Convert data to matrix for XGBoost
    train_data <- df.matrix.omics %>% dplyr::select(-feature) %>% as.matrix()
    train_labels <- df.matrix.omics$feature
    
    # Create stratified folds for cross-validation
    xpectr::set_test_seed(1)
    train_set <- metadata %>% filter(sampleID %in% row.names(df.matrix.omics))
    train_set$idpatient <- as.factor(train_set$idpatient)
    train_set <- groupdata2::fold(train_set, k = nfolds, cat_col = panel, id_col = "idpatient")
    
    # Prepare folds

    train_set <- train_set %>% arrange(.folds)
    train_set$folds.name <- paste0("Fold", train_set$.folds)
    fold.name <- unique(train_set$folds.name)
    folds <- list()
    for (f in fold.name) {
      df.fold <- filter(train_set, folds.name == f)
      fold.id <- df.fold$sampleID
      fold <- which(row.names(train_data) %in% fold.id)
      names(fold)<-fold.id
      folds[[f]] <- fold
    }
    
    
    # Evaluate fold 1
    x <- folds$Fold1
    training_fold <- train_data[-x, ]
    val_fold <- train_data[x, ]
    train_labels.fold <- train_labels[-x]
    val_labels.fold <- train_labels[x]
    
    dtrain.fold <- xgb.DMatrix(data = training_fold, label = train_labels.fold)
    dval.fold <- xgb.DMatrix(data = val_fold, label = val_labels.fold)
    
    # Tuning in fold 1 and generate loss curve
    watchlist <- list(train = dtrain.fold, test = dval.fold)
    model <- xgb.train(
      data = dtrain.fold,
      watchlist = watchlist,
      nrounds = 20,
      eval_metric = "logloss",
      alpha = 10,
      min_child_weight = 1,
      lambda = 1,
      objective = "binary:logistic"
    )
    out.f1 <- model$evaluation_log
    out.f1.long <- melt(out.f1, id.vars = "iter")
    loss.curve <- ggplot(out.f1.long, aes(x = iter, y = value, color = variable)) +
      geom_line() +
      theme_bw()
    loss.list[[omi]] <- loss.curve
    
    # Initialize score tracking
    score_list <- data.frame(
      folds = seq(1, length(fold.name)),
      scores = rep(0, length(fold.name)),
      AUC = rep(0, length(fold.name)),
      AUC.top10 = rep(0, length(fold.name)),
      rounds = rep(0, length(fold.name))
    )
    
    
    preds_list <- data.frame(matrix(ncol = 5, nrow = nrow(val_fold)))
    colnames(preds_list) <- paste("Repeat", 1:5)
    preds_list.top <- data.frame(matrix(ncol = 5, nrow = nrow(val_fold)))
    colnames(preds_list.top) <- paste("Repeat", 1:5)
    folds.list[[omi]] <- folds
    
    # Perform cross-validation
    
    for (i in names(folds)) {
     
      x <- folds[[i]]
      training_fold <- train_data[-x, ]
      
      val_fold <- train_data[x, ]
      train_labels.fold <- train_labels[-x]
      val_labels.fold <- train_labels[x]
      

      # Skip fold if validation is all zeros
      val_labels.fold.na <- ifelse(val_labels.fold == 0, NA, val_labels.fold)
      
      if (!all(is.na(val_labels.fold.na))) {
        
        labels <- row.names(train_data)
        train.ID.fold <- labels[-x]
        val.ID.fold <- labels[x]
        train.labels.list[[i]] <- train.ID.fold
        val.labels.list[[i]] <- val.ID.fold
        names(val_labels.fold)<-val.ID.fold
        
        # Convert to XGBoost DMatrix
        training_xgb <- xgb.DMatrix(data = training_fold, label = train_labels.fold)
        testing_xgb <- xgb.DMatrix(data = val_fold, label = val_labels.fold)
        
        set.seed(seed)
        
        negative_cases <- sum(train_labels.fold == FALSE)
        positive_cases <- sum(train_labels.fold == TRUE)
        scale_pos_weight <- negative_cases / positive_cases
        
        model <- xgboost(
          data = training_xgb,
          scale_pos_weight = scale_pos_weight,
          nround = 500,
          alpha = 10,
          lambda = 1,
          eta = 0.3,
          gamma = 0,
          max_depth = 6,
          early_stopping_rounds = 20,
          eval_metric = "auc",
          objective = "binary:logistic",
          nthread = 15
        )
        
        
        # Store results
        model_list[[i]] <- model
        
        best_rounds <- model$best_iteration
        labelsValidation.list[[i]] <- val_labels.fold
        score_list[which(i == names(folds)), "rounds"] <- best_rounds
        score_list[which(i == names(folds)), "scores"] <- model$evaluation_log[[2]][score_list[which(i == names(folds)), "rounds"]]
        
        
        # Make predictions
        pred <- predict(model, testing_xgb)
        names(pred)<-val.ID.fold
        preds_list[folds[[i]], (((which(i == names(folds)) - 1) %/% nfolds) + 1)] <- pred
        predictions.list[[i]] <- pred
        
        # Calculate confusion matrix
        
        pred.fact <- ifelse(pred >= 0.5, 1, 0)
        
        val_labels.fold.factor<-factor(val_labels.fold, levels = c(0,1))
        pred.fact<-factor(pred.fact, levels = c(0,1))
        
        conf.mx <- confusionMatrix(pred.fact, val_labels.fold.factor)
        confusionMatrix.list[[i]] <- conf.mx
        
        roc.fold <- roc(val_labels.fold, pred)
        auc.fold <- Metrics::auc(actual = val_labels.fold, predicted = pred)
        score_list[which(i == names(folds)), "AUC"] <- auc.fold
        
        importance_matrix <- xgb.importance(model = model)
        feature.importance.list[[i]] <- importance_matrix
        
        # Top feature reduction model
        if(ncol(train_data)>=top){top.n=top}else{top.n=ncol(train_data)}
        
        for (u in 2:top.n) {
          
            if (nrow(importance_matrix) >= u) {
              
              
              features <- paste0("top.", u)
              
              top_features <- importance_matrix$Feature[1:u]
              training_fold.fs <- training_fold[, top_features]
              val_fold.fs <- val_fold[, top_features]
              
              # Convert to XGBoost DMatrix
              dtrain.fold.fs <- xgb.DMatrix(data = training_fold.fs, label = train_labels.fold)
              dval.fold.fs <- xgb.DMatrix(data = val_fold.fs, label = val_labels.fold)
              
              model_tuned.fs <- xgboost(
                data = dtrain.fold.fs,
                scale_pos_weight = scale_pos_weight,
                nround = 500,
                alpha = 10,
                lambda = 1,
                eta = 0.3,
                gamma = 0,
                max_depth = 6,
                early_stopping_rounds = 20,
                objective = "binary:logistic",
                nthread = 15
              )
              
              # Store top feature results
              model_list.top[[features]][[i]] <- model_tuned.fs
              best_rounds <- which.min(model_tuned.fs$evaluation_log[[2]])
              
              pred.fs <- predict(model_tuned.fs, dval.fold.fs)
              names(pred.fs)<-val.ID.fold
              preds_list.top[folds[[i]], (((which(i == names(folds)) - 1) %/% nfolds) + 1)] <- pred.fs
              predictions.list.top[[features]][[i]] <- pred.fs
              
              # Calculate confusion matrix
              pred.fact.fs <- ifelse(pred.fs >= 0.5, 1, 0)
              
              val_labels.fold<-factor(val_labels.fold, levels = c(0,1))
              pred.fact.fs<-factor(pred.fact.fs, levels = c(0,1))
              
              conf.mx.fs <- confusionMatrix(factor(pred.fact.fs), factor(val_labels.fold))
              confusionMatrix.list.top[[features]][[i]] <- conf.mx.fs
              
              
              auc.fold.fs <- Metrics::auc(actual = val_labels.fold, predicted = pred.fs)
              score_list[which(i == names(folds)), "AUC.top10"] <- auc.fold.fs
              
              importance_matrix.fs <- xgb.importance(model = model_tuned.fs)
              feature.importance.list.top[[features]][[i]] <- importance_matrix.fs
            } else {
              model_list.top[[features]][[i]] <- NULL
              feature.importance.list.top[[features]][[i]] <- NULL
              predictions.list.top[[features]][[i]] <- NULL
              confusionMatrix.list.top[[features]][[i]] <- NULL
            }
         
        }
      }
    }
    
    score_list.all[[omi]] <- score_list
    predictions.list.all[[omi]] <- predictions.list
    predictions.list.all.top[[omi]] <- predictions.list.top
    model_list.all[[omi]] <- model_list
    model_list.all.top[[omi]] <- model_list.top
    feature.importance.list.all[[omi]] <- feature.importance.list
    feature.importance.list.top.all[[omi]] <- feature.importance.list.top
    confusionMatrix.list.all[[omi]] <- confusionMatrix.list
    confusionMatrix.lit.top.all[[omi]] <- confusionMatrix.list.top
    train.labels.list.out[[omi]] <- train.labels.list
    val.labels.list.out[[omi]] <- val.labels.list
    
    
    
    # Calculate average AUC
    out <- cvAUC(predictions.list, labelsValidation.list)
    
    id.list <- lapply(predictions.list, names)
    out.ci <- ci.pooled.cvAUC(predictions.list, labelsValidation.list, label.ordering = NULL,  ids=id.list, confidence = 0.95)
    out.list[[omi]] <- out
    
    cvAUC.list[[omi]] <- list(cvAUC = out, out.ci = out.ci)
    cvAUC.value[[omi]] <- out$cvAUC
    cvAUC.se[[omi]] <- out.ci$se
    
    # Calculate average AUC for top features
    
    if (!is.null(names(predictions.list.top))) {
      
      # Remove empty folds
      predictions.list.top <- predictions.list.top[sapply(predictions.list.top, function(x) !(is.list(x) && length(x) == 0))]
      
      #Ensure predictions and lavidadtion folds have the same length
    for (u in names(predictions.list.top)) {

        predictions.list.top.u <- predictions.list.top[[u]]
    
        labelsValidation.list.fs<-labelsValidation.list[names(predictions.list.top.u)]
        
        labelsValidation.list.fs2 <- lapply(names(labelsValidation.list.fs), function(l) {
          labelsValidation.list.fs[[l]][names(predictions.list.top.u[[l]])]
        })
        
        names(labelsValidation.list.fs2) <- names(labelsValidation.list.fs)
        
        # Calculate average AUC for top features
        
        
         out.fs <- cvAUC(predictions.list.top.u, labelsValidation.list.fs2)
        
        id.list.fs <- lapply(predictions.list.top.u, names)
        out.ci.fs <- ci.pooled.cvAUC(predictions.list.top.u, labelsValidation.list.fs2, ids=id.list.fs, confidence = 0.95)
        
        out.list.top[[omi]][[u]] <- out.fs
        cvAUC.list.top[[omi]][[u]] <- list(cvAUC = out.fs, out.ci = out.ci.fs)
  
    }
      
      }
    
    
  }
  
  list.return <- list(
    loss.list.all = loss.list,
    cvAUC.list.all = cvAUC.list,
    cvAUC.list.top.all= cvAUC.list.top,
    feature.importance.list = feature.importance.list.all,
    feature.importance.list.top = feature.importance.list.top.all,
    cvAUC.value = cvAUC.value,
    cvAUC.se = cvAUC.se,
    model_list.all = model_list.all,
    model_list.all.top = model_list.all.top,
    predictions.list = predictions.list.all,
    predictions.list.top = predictions.list.all.top,
    labelsValidation.list = labelsValidation.list,
    folds = folds.list,
    out = out.list,
    out.top = out.list.top,
    score_list.all = score_list.all,
    confusionMatrix.list.all = confusionMatrix.list.all,
    confusionMatrix.lit.top.all = confusionMatrix.lit.top.all,
    train.labels.ID.out = train.labels.list.out,
    val.labels.ID.out = val.labels.list.out
  )
  
  return(list.return)
}

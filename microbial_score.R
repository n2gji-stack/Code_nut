#---------------Feature selection----------------------#
library(caret)
library(Boruta)
library(randomForest)
library(gridExtra)
library(glmnet) 
library(dplyr)  
set.seed(123)
Microbiome <- Annot_Taxon$taxon_new
Microbiome <- Microbiome[grep("^s_", Microbiome)]
species_data <- ToBeUsed.Taxon[, Microbiome]
min_abundance <- 0.0001
min_prevalence <- 0.1
prevalence <- colMeans(species_data > 0) 
mean_abundance <- colMeans(species_data) 
species_to_keep <- names(prevalence)[prevalence >= min_prevalence & mean_abundance >= min_abundance]
filtered_species_data <- species_data[, species_to_keep] #385
ast_transformed_data <- asin(sqrt(filtered_species_data))
######################################################################################
######################################################################################
run_boruta_rf_pipeline <- function(target_var_name, data_source, predictors_df) {
  
  cat("\n####################################################################\n")
  cat("STARTING ANALYSIS FOR:", target_var_name, "\n")
  cat("####################################################################\n")
  if(is.null(data_source[[target_var_name]])) {
    stop(paste("Error: Variable", target_var_name, "not found in ToBeUsed.Taxon"))
  }
  
  tmp_data <- data.frame(
    target = data_source[[target_var_name]], 
    predictors_df
  )
  tmp_data <- na.omit(tmp_data)
  
  set.seed(123)
  folds <- createFolds(tmp_data$target, k = 5, returnTrain = TRUE)
  
  all_predictions <- data.frame()
  feature_frequency <- c() 
  
  for(i in 1:5){
    cat("  >> Running Fold", i, "/ 5 ...\n")
    
    train_idx <- folds[[i]]
    train_data <- tmp_data[train_idx, ]
    test_data  <- tmp_data[-train_idx, ] 
    
    set.seed(123 + i) 
    fs_boruta <- Boruta(x = train_data[, -1], y = train_data$target, 
                        doTrace = 0, maxRuns = 500, 
                        getImp = getImpLegacyRfGini)
    selected_feats <- getSelectedAttributes(fs_boruta, withTentative = FALSE)
    
    if(length(selected_feats) == 0){
      selected_feats <- getSelectedAttributes(fs_boruta, withTentative = TRUE)
    }
    
    if(length(selected_feats) == 0){
      mean_prediction <- mean(train_data$target, na.rm = TRUE)
      preds <- rep(mean_prediction, nrow(test_data))
    } else {
      feature_frequency <- c(feature_frequency, selected_feats)
      
      x_train <- train_data[, selected_feats, drop = FALSE]
      x_test  <- test_data[, selected_feats, drop = FALSE]
      
      my_mtry <- min(floor(sqrt(ncol(x_train))), ncol(x_train))
      if(my_mtry < 1) my_mtry <- 1
      
      rf_model <- randomForest(x = x_train, 
                               y = train_data$target,
                               ntree = 500,
                               mtry = my_mtry)
      
      preds <- predict(rf_model, newdata = x_test)
    }
    
    fold_res <- data.frame(Observed = test_data$target, 
                           Predicted = preds, 
                           Fold = i)
    all_predictions <- rbind(all_predictions, fold_res)
  }
  
  cat("  >> CV Loop Finished. Calculating metrics...\n")
  
  cv_pearson <- cor.test(all_predictions$Observed, all_predictions$Predicted, method = "pearson")
  cv_spearman <- cor.test(all_predictions$Observed, all_predictions$Predicted, method = "spearman")
  
  feat_table <- sort(table(feature_frequency), decreasing = TRUE)
  stability_threshold <- 2
  final_features <- names(feat_table)[feat_table >= stability_threshold]
  
  cat("  >> Stable features selected:", length(final_features), "\n")
  if(length(final_features) > 0) {
    cat("     ", paste(head(final_features, 5), collapse=", "), "...\n")
  }
  x_final <- tmp_data[, final_features, drop = FALSE]
  y_final <- tmp_data$target
  best_mtry_final <- 1
  if(length(final_features) > 1) {
    tune_final_mtry <- function(feature_cols, target) {
      n_feat <- ncol(feature_cols)
      mtry_grid <- 1:min(n_feat, 15) 
      errors <- sapply(mtry_grid, function(m) {
        set.seed(1234)
        rf <- randomForest(x = feature_cols, y = target, 
                           mtry = m, ntree = 500, importance = FALSE)
        return(rf$mse[length(rf$mse)])
      })
      return(mtry_grid[which.min(errors)])
    }
    best_mtry_final <- tune_final_mtry(x_final, y_final)
  }
  
  cat("  >> Training Final Model with mtry =", best_mtry_final, "\n")
  
  set.seed(1234)
  final_rf_model <- randomForest(x = x_final, 
                                 y = y_final,
                                 mtry = best_mtry_final,
                                 ntree = 1000,
                                 importance = TRUE)
  
  print(cv_pearson)
  return(list(
    Model = final_rf_model,
    CV_Predictions = all_predictions,
    Stable_Features = final_features,
    CV_Pearson = cv_pearson,
    CV_Spearman = cv_spearman,
    Data_Used = tmp_data # 如果后续需要用到清洗后的数据
  ))
}
# ==============================================================================
# 
# ==============================================================================
# all_nut_r_log_zscore
# ------------------------------------------------------------------------------
result_all <- run_boruta_rf_pipeline("all_nut_r_log_zscore", ToBeUsed.Taxon, ast_transformed_data)
Final_RF_Model_all <- result_all$Model
gc() 
cor.test(Final_RF_Model_all$predicted,ToBeUsed.Taxon$all_nut)
cor.test(Final_RF_Model_all$predicted,ToBeUsed.Taxon$all_nut,method = "spearman")
#  treenut_r_log_zscore
# ------------------------------------------------------------------------------
result_tree <- run_boruta_rf_pipeline("treenut_r_log_zscore", ToBeUsed.Taxon, ast_transformed_data)
Final_RF_Model_tree <- result_tree$Model
gc()
cor.test(Final_RF_Model_tree$predicted,ToBeUsed.Taxon$treenut_r_log_zscore,method = "spearman")
cor.test(Final_RF_Model_tree$predicted,ToBeUsed.Taxon$treenut_r_log_zscore)
#  peanut_r_logzscore

# ------------------------------------------------------------------------------
#result_pea <- run_boruta_rf_pipeline("peanut_r_log_zscore", ToBeUsed.Taxon, ast_transformed_data)
#NULL
#  walnut_r_logzscore
# ------------------------------------------------------------------------------
result_wn <- run_boruta_rf_pipeline("walnut_r_log_zscore", ToBeUsed.Taxon, ast_transformed_data)
Final_RF_Model_wn <- result_wn$Model
gc()
cor.test(Final_RF_Model_wn$predicted,ToBeUsed.Taxon$walnut)
cor.test(Final_RF_Model_wn$predicted,ToBeUsed.Taxon$walnut,method = "spearman")
# --- Helper to safely extract predictions ---
get_preds_safely <- function(model_obj, original_data) {
  if (!is.null(model_obj) && !is.null(model_obj$predicted)) {
    return(model_obj$predicted)
  } else {
    warning("Model is NULL, returning NAs.")
    return(rep(NA, nrow(original_data)))
  }
}
cor.test(ToBeUsed.Taxon$all_nut, cv_df$Predicted, method = "pearson")
cor.test(cv_df$Observed, cv_df$Predicted, method = "spearman")
ToBeUsed.Taxon$pred_wn <- get_preds_safely(Final_RF_Model_wn, ToBeUsed.Taxon)
ToBeUsed.Taxon$pred_tn <- get_preds_safely(Final_RF_Model_tree, ToBeUsed.Taxon)
ToBeUsed.Taxon$pred_an <- get_preds_safely(Final_RF_Model_all, ToBeUsed.Taxon)
save(Final_RF_Model_wn,Final_RF_Model_tree,Final_RF_Model_all,file = "/udd/n2gji/micro/nut_final/nut_Dec/data/RF.RData")
gut_an<-rownames(importance(Final_RF_Model_all))
gut_tn<-rownames(importance(Final_RF_Model_tree))
gut_wn<-rownames(importance(Final_RF_Model_wn))

############################################################################################################
############################################################################################################
#######################################SHAP############################################################
############################################################################################################
library(shapviz)
library(kernelshap)
#Treenut
target_variable_name <- "treenut_r_log_zscore" 
data_tn<- cbind(target = ToBeUsed.Taxon[[target_variable_name]],ast_transformed_data[rownames(importance(Final_RF_Model_tree))] )
pfun <- function(object, newdata) {
  predict(object, newdata = newdata)
}

X_data <- data_tn[, -c(1)] 
set.seed(1234)
if(nrow(X_data) > 100){
  bg_X <- X_data[sample(nrow(X_data), 100), ]
} else {
  bg_X <- X_data
}
shap_values_tn <- kernelshap(
  object = Final_RF_Model_tree,
  X = X_data,
  pred_fun = pfun,
  bg_X = bg_X,
  exact = FALSE  
)
sv_tn <- shapviz(shap_values_tn)
print(sv_tn)
name_mapping <- c(
  "s_0564" = "Blautia Obeum",
  "s_0636" = "UBA11774 sp003507655",
  "s_0864" = "Roseburia Faecis",
  "s_1158" = "Lachnospira sp NSJ 43",
  "s_1862" = "Adlercreutzia Equolifaciens")
sv_original_names <- colnames(sv_tn$S)
sv_new_names <- ifelse(sv_original_names %in% names(name_mapping), name_mapping[sv_original_names], sv_original_names)
sv_new_names
colnames(sv_tn$S) <- sv_new_names
colnames(sv_tn$X) <- sv_new_names
imp_tn<-sv_importance(sv_tn, kind = "bar")
imp_tn<-imp_tn+ 
  ggplot2::theme(
    axis.text.y = ggplot2::element_text(face = "italic")
  )
imp_tn
ggsave(
  filename = "/udd/n2gji/micro/nut_final/nut_Dec/output/shap_bar_tn_1215.pdf",
  plot = imp_tn,
  width = 6,
  height = 3,
  units = "in"
)
bee_tn<-sv_importance(sv_tn,kind = "beeswarm")
bee_tn<-bee_tn+ 
  ggplot2::theme(
    axis.text.y = ggplot2::element_text(face = "italic")
  )
bee_tn
ggsave(
  filename = "/udd/n2gji/micro/nut_final/nut_Dec/output/shap_bee_tn_1215.pdf",
  plot = bee_tn,
  width = 6,
  height = 3,
  units = "in"
)
#####################################WALNUT##########################################################
target_variable_name <- "walnut_r_log_zscore" 
data_wn<- cbind(target = ToBeUsed.Taxon[[target_variable_name]],ast_transformed_data[rownames(importance(Final_RF_Model_wn))] )
X_data <- data_wn[, -c(1)] 
set.seed(1234)
if(nrow(X_data) > 100){
  bg_X <- X_data[sample(nrow(X_data), 100), ]
} else {
  bg_X <- X_data
}
shap_values <- kernelshap(
  object = Final_RF_Model_wn,
  X = X_data,
  pred_fun = pfun,
  bg_X = bg_X
)
sv_wn <- shapviz(shap_values)
print(sv_wn)
name_mapping_wn <- c(
  "s_0636" = "UBA11774 sp003507655",
  "s_1312" = "Dysosmobacter sp BX15",
  "s_0789" = "Dysosmobacter sp900544615",
  "s_2130" = "Eubacteriales SGB15145")
sv_original_names <- colnames(sv_wn$S)
sv_new_names <- ifelse(sv_original_names %in% names(name_mapping_wn), name_mapping_wn[sv_original_names], sv_original_names)
colnames(sv_wn$S) <- sv_new_names
colnames(sv_wn$X) <- sv_new_names
imp_wn<-sv_importance(sv_wn, kind = "bar")
imp_wn<-imp_wn+ 
  ggplot2::theme(
    axis.text.y = ggplot2::element_text(face = "italic")
  )
imp_wn
ggsave(
  filename = "/udd/n2gji/micro/nut_final/nut_Dec/output/shap_bar_wn_1225.pdf",
  plot = imp_wn,
  width = 6,
  height = 3,
  units = "in"
)
bee_wn<-sv_importance(sv_wn,kind = "beeswarm")
bee_wn<-bee_wn+ 
  ggplot2::theme(
    axis.text.y = ggplot2::element_text(face = "italic")
  )
bee_wn
ggsave(
  filename = "/udd/n2gji/micro/nut_final/nut_Dec/output/shap_bee_wn_1225.pdf",
  plot = bee_wn,
  width = 6,
  height = 3,
  units = "in"
)
#####################################Total nuts##########################################################
target_variable_name <- "all_nut_r_log_zscore" 
data_an<- cbind(target = ToBeUsed.Taxon[[target_variable_name]],ast_transformed_data[rownames(importance(Final_RF_Model_all))] )
X_data <- data_an[, -c(1)] 
set.seed(1234)
if(nrow(X_data) > 100){
  bg_X <- X_data[sample(nrow(X_data), 100), ]
} else {
  bg_X <- X_data
}
shap_values_an <- kernelshap(
  object = Final_RF_Model_all,
  X =  X_data ,
  pred_fun = pfun,
  bg_X = bg_X
)
sv_an <- shapviz(shap_values_an)
sv_an
sv_original_names <- colnames(sv_an$S)
name_mapping_an <- c(
  "s_0564"="Blautia Obeum",
  "s_0636" = "UBA11774 sp003507655",
  "s_1862" = "Adlercreutzia Equolifaciens"
)
sv_new_names <- ifelse(sv_original_names %in% names(name_mapping_an), name_mapping_an[sv_original_names], sv_original_names)
colnames(sv_an$S) <- sv_new_names
colnames(sv_an$X) <- sv_new_names
imp_an<-sv_importance(sv_an, kind = "bar")
imp_an<-imp_an+ 
  ggplot2::theme(
    axis.text.y = ggplot2::element_text(face = "italic")
  )
imp_an
ggsave(
  filename = "/udd/n2gji/micro/nut_final/nut_Dec/output/shap_bar_an_1215.pdf",
  plot = imp_an,
  width = 6,
  height = 3,
  units = "in"
)
bee_an<-sv_importance(sv_an,kind = "beeswarm")
bee_an<-bee_an+ 
  ggplot2::theme(
    axis.text.y = ggplot2::element_text(face = "italic")
  )
bee_an
ggsave(
  filename = "/udd/n2gji/micro/nut_final/nut_Dec/output/shap_bee_an_1215.pdf",
  plot = bee_an,
  width = 6,
  height = 3,
  units = "in"
)
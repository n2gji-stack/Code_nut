#install.packages('ggalt')
#install.packages('ggpubr')
#install.packages("ggnewscale")
#BiocManager::install("treeio")
#BiocManager::install("ggtree")
#BiocManager::install("ggtreeExtra")
library(chanmetab)
library(Biobase)
library(dplyr)
library(stringr)
library(naniar)
library(ggplot2)
library(naniar)
library(forcats)
library(ggsci)
library(glmnet)
library(caret)
library(ggvenn)
library(ggrepel)
library(tidyverse)
library(ggtree)
library(treeio)
library(ape)
library(ggnewscale)
library(ggtreeExtra)
library(MetBrewer)
library(viridis)
library(RColorBrewer)
library(patchwork)
library(ggpubr)
library(tableone)
library(survival)
library(survminer)
library(ggalt)
library(treeio)
library(BiocManager)
# 重新安装
# 再次尝试加载
library(ggnewscale)
library(ggtreeExtra)
library(dplyr)
library(data.table)
library(tidyverse)
library(circlize)
library(ggsci)
library(ComplexHeatmap)
##########################################Using predicted micro bial score to select metabolites########################################
AdjVars<-c("ageyr","bmi_bld","totMETs_paq", "ahei_g_nnut",
           "probio_2m_fec","antibio_12m_fec","acid_2m_fec","colsc_2m_fec",
           "bs_1","bs_2","bs_3","calor_fo_dr_wtavg")

hmdb_vars <- grep("^HMDB", colnames(ToBeUsed.Taxon),value = TRUE)
results_met=data.frame(met=hmdb_vars)
met_aligned=results_met %>% pull(met)

for (i in 1:295) {
  fit=lm(as.formula(paste0(met_aligned[i], "~ walnut_r_log_zscore +",paste(AdjVars,collapse = " + "))), data=ToBeUsed.Taxon)
  summ=summary(fit)
  results_met[i, "beta_wnut_adj"]=round(summ$coefficients[2,1],4)
  results_met[i, "se_wnut_adj"]=summ$coefficients[2,2]
  results_met[i, "p_wnut_adj"]=summ$coefficients[2,4]
  results_met[i, "n_wnut_adj"]=nrow(model.frame(fit))
  
  fit=lm(as.formula(paste0(met_aligned[i], "~ all_nut_r_log_zscore +",paste(AdjVars,collapse = " + "))), data=ToBeUsed.Taxon)
  summ=summary(fit)
  results_met[i, "beta_anut_adj"]=round(summ$coefficients[2,1],4)
  results_met[i, "se_anut_adj"]=summ$coefficients[2,2]
  results_met[i, "p_anut_adj"]=summ$coefficients[2,4]
  results_met[i, "n_anut_adj"]=nrow(model.frame(fit))
  
  fit=lm(as.formula(paste0(met_aligned[i], "~ treenut_r_log_zscore  +",paste(AdjVars,collapse = " + "))), data=ToBeUsed.Taxon)
  summ=summary(fit)
  results_met[i, "beta_tnut_adj"]=round(summ$coefficients[2,1],4)
  results_met[i, "se_tnut_adj"]=summ$coefficients[2,2]
  results_met[i, "p_tnut_adj"]=summ$coefficients[2,4]
  results_met[i, "n_tnut_adj"]=nrow(model.frame(fit))
  
  fit=lm(as.formula(paste0(met_aligned[i], "~ peanut_r_log_zscore  +",paste(AdjVars,collapse = " + "))), data=ToBeUsed.Taxon)
  summ=summary(fit)
  results_met[i, "beta_pnut_adj"]=round(summ$coefficients[2,1],4)
  results_met[i, "se_pnut_adj"]=summ$coefficients[2,2]
  results_met[i, "p_pnut_adj"]=summ$coefficients[2,4]
  results_met[i, "n_pnut_adj"]=nrow(model.frame(fit))
  
  fit=lm(as.formula(paste0(met_aligned[i], "~ peanut_exbut_r_log_zscore  +",paste(AdjVars,collapse = " + "))), data=ToBeUsed.Taxon)
  summ=summary(fit)
  results_met[i, "beta_pnutex_adj"]=round(summ$coefficients[2,1],4)
  results_met[i, "se_pnutex_adj"]=summ$coefficients[2,2]
  results_met[i, "p_pnutex_adj"]=summ$coefficients[2,4]
  results_met[i, "n_pnutex_adj"]=nrow(model.frame(fit))
  
}
results_met$p_wnut_adj_fdr <- p.adjust(results_met$p_wnut_adj, method = "fdr")
results_met$p_anut_adj_fdr <- p.adjust(results_met$p_anut_adj, method = "fdr")
results_met$p_tnut_adj_fdr <- p.adjust(results_met$p_tnut_adj, method = "fdr")
results_met$p_pnut_adj_fdr <- p.adjust(results_met$p_pnut_adj, method = "fdr")
results_met$p_pnutex_adj_fdr <- p.adjust(results_met$p_pnutex_adj, method = "fdr")
results_met<- results_met %>% mutate(
  wnut_intake=case_when(beta_wnut_adj>0 & p_wnut_adj_fdr<0.15  ~ 'positive',
                        beta_wnut_adj<0 & p_wnut_adj_fdr<0.15  ~ 'negative',
                        T ~ 'insignificant'),
  anut_intake=case_when(beta_anut_adj>0 & p_anut_adj_fdr<0.15  ~ 'positive',
                        beta_anut_adj<0 & p_anut_adj_fdr<0.15  ~ 'negative',
                        T ~ 'insignificant'),
  tnut_intake=case_when(beta_tnut_adj>0 & p_tnut_adj_fdr<0.15  ~ 'positive',
                        beta_tnut_adj<0 & p_tnut_adj_fdr<0.15  ~ 'negative',
                        T ~ 'insignificant'),
  pnut_intake=case_when(beta_pnut_adj>0 & p_pnut_adj_fdr<0.15  ~ 'positive',
                        beta_pnut_adj<0 & p_pnut_adj_fdr<0.15  ~ 'negative',
                        T ~ 'insignificant'),
  pnutex_intake=case_when(beta_pnutex_adj>0 & p_pnutex_adj_fdr<0.15  ~ 'positive',
                          beta_pnutex_adj<0 & p_pnutex_adj_fdr<0.15  ~ 'negative',
                          T ~ 'insignificant')
)

f_lvs_match<-f_lvs %>% mutate(met=rownames(f_lvs)) %>% 
  dplyr::select("met",'method',"metabolite_name",'biochemical_name',"class_metabolon","sub_class_metabolon","super_class_metabolon")
results_met<-merge(results_met, f_lvs_match, by = 'met', all = FALSE, sort = TRUE)
names(results_met)
results_mlvs<-results_met
write.csv(results_mlvs,"/udd/n2gji/micro/nut_final/nut_Dec/output/lm_nut_met_mlvs.csv")
##############micro score & metabolites#################
ToBeUsed.Taxon$pred_wn_s<-scale(ToBeUsed.Taxon$pred_wn)
ToBeUsed.Taxon$pred_tn_s<-scale(ToBeUsed.Taxon$pred_tn)
ToBeUsed.Taxon$pred_an_s<-scale(ToBeUsed.Taxon$pred_an)
ToBeUsed.Taxon$ahei_g_nnut<-as.factor(ToBeUsed.Taxon$ahei_g_nnut)
hmdb_vars <- grep("^HMDB", colnames(ToBeUsed.Taxon),value = TRUE)
results_met=data.frame(met=hmdb_vars)
met_aligned=results_met %>% pull(met)
for (i in 1:295) {
  fit=lm(as.formula(paste0(met_aligned[i], "~ pred_wn_s +",paste(AdjVars,collapse = " + "))), data=ToBeUsed.Taxon)
  summ=summary(fit)
  results_met[i, "beta_wnut_adj"]=round(summ$coefficients[2,1],4)
  results_met[i, "se_wnut_adj"]=summ$coefficients[2,2]
  results_met[i, "p_wnut_adj"]=summ$coefficients[2,4]
  results_met[i, "n_wnut_adj"]=nrow(model.frame(fit))
  
  fit=lm(as.formula(paste0(met_aligned[i], "~ pred_an_s +",paste(AdjVars,collapse = " + "))), data=ToBeUsed.Taxon)
  summ=summary(fit)
  results_met[i, "beta_anut_adj"]=round(summ$coefficients[2,1],4)
  results_met[i, "se_anut_adj"]=summ$coefficients[2,2]
  results_met[i, "p_anut_adj"]=summ$coefficients[2,4]
  results_met[i, "n_anut_adj"]=nrow(model.frame(fit))
  
  fit=lm(as.formula(paste0(met_aligned[i], "~ pred_tn_s  +",paste(AdjVars,collapse = " + "))), data=ToBeUsed.Taxon)
  summ=summary(fit)
  results_met[i, "beta_tnut_adj"]=round(summ$coefficients[2,1],4)
  results_met[i, "se_tnut_adj"]=summ$coefficients[2,2]
  results_met[i, "p_tnut_adj"]=summ$coefficients[2,4]
  results_met[i, "n_tnut_adj"]=nrow(model.frame(fit))
}
results_met$p_wnut_adj_fdr <- p.adjust(results_met$p_wnut_adj, method = "fdr")
results_met$p_anut_adj_fdr <- p.adjust(results_met$p_anut_adj, method = "fdr")
results_met$p_tnut_adj_fdr <- p.adjust(results_met$p_tnut_adj, method = "fdr")
results_met<- results_met %>% mutate(
  wnut_intake=case_when(beta_wnut_adj>0 & p_wnut_adj_fdr<0.15  ~ 'positive',
                        beta_wnut_adj<0 & p_wnut_adj_fdr<0.15  ~ 'negative',
                        T ~ 'insignificant'),
  anut_intake=case_when(beta_anut_adj>0 & p_anut_adj_fdr<0.15  ~ 'positive',
                        beta_anut_adj<0 & p_anut_adj_fdr<0.15  ~ 'negative',
                        T ~ 'insignificant'),
  tnut_intake=case_when(beta_tnut_adj>0 & p_tnut_adj_fdr<0.15  ~ 'positive',
                        beta_tnut_adj<0 & p_tnut_adj_fdr<0.15  ~ 'negative',
                        T ~ 'insignificant')
  
)
results_met<-merge(results_met, f_lvs_match, by = 'met', all = FALSE, sort = TRUE)
results_micro<-results_met
write.csv(results_micro,"/udd/n2gji/micro/nut_final/nut_Dec/output/lm_mrs_met_mlvs.csv")
##############
table(
  results_micro$met[results_micro$anut_intake != "insignificant"] %in%
    results_mlvs$met[results_mlvs$anut_intake != "insignificant"]
)#0
table(
  results_micro$met[results_micro$wnut_intake != "insignificant"] %in%
    results_mlvs$met[results_mlvs$wnut_intake != "insignificant"]
)#FALSE  TRUE #16    42 
table(
  results_micro$met[results_micro$tnut_intake != "insignificant"] %in%
    results_mlvs$met[results_mlvs$tnut_intake != "insignificant"]
)# 39    38 

distict_wn <- subset(
  results_micro,
  wnut_intake != "insignificant" &
    met %in% results_mlvs$met[results_mlvs$wnut_intake == "insignificant"]
)
distict_tn <- subset(
  results_micro,
  tnut_intake != "insignificant" &
    met %in% results_mlvs$met[results_mlvs$tnut_intake == "insignificant"]
)
##############Build metabolite score##################

metabolite_cols <- grep("^HMDB", names(ToBeUsed.Taxon), value = TRUE)


covariate_cols_species <- c("ageyr","bmi_bld","totMETs_paq", "ahei_g_nnut",
                            "probio_2m_fec","antibio_12m_fec","acid_2m_fec","colsc_2m_fec",
                            "bs_1","bs_2","bs_3","calor_fo_dr_wtavg")
run_loocv_analysis <- function(data, 
                               species_col,    # e.g., "pred_an"
                               nut_col,        # e.g., "all_nut_r_log_zscore"
                               metab_cols, 
                               covar_cols, 
                               fdr_thresh = 0.15) {
  
  cat(paste0("\n=== 开始分析组合: Species=[", species_col, "] -> Nut=[", nut_col, "] ===\n"))
  cat(paste0("时间: ", Sys.time(), "\n"))
  
  # 初始化存储容器
  n_samples <- nrow(data)
  res <- list(
    preds_strat1 = numeric(n_samples),
    preds_strat2 = numeric(n_samples),
    preds_strat3 = numeric(n_samples),
    preds_strat4 = numeric(n_samples),
    
    coefs_strat1 = vector("list", n_samples),
    coefs_strat2 = vector("list", n_samples),
    coefs_strat3 = vector("list", n_samples),
    coefs_strat4 = vector("list", n_samples),
    
    selected_metabs_strat1 = vector("list", n_samples),
    selected_metabs_strat2 = vector("list", n_samples),
    selected_metabs_strat3 = vector("list", n_samples),
    selected_metabs_strat4 = vector("list", n_samples)
  )
  
  # LOOCV 循环
  for (i in 1:n_samples) {
    if (i %% 50 == 0) cat(paste("  ...处理折叠", i, "/", n_samples, "\n"))
    
    train_data <- data[-i, ]
    test_data <- data[i, , drop = FALSE]
    
    X_train_metab <- as.matrix(train_data[, metab_cols])
    X_test_metab <- as.matrix(test_data[, metab_cols])
    y_train <- train_data[[nut_col]]
    
    # --- 策略 1: ENet ---
    cv_fit_enet <- tryCatch({
      cv.glmnet(X_train_metab, y_train, alpha = 0.5, nfolds = 10, family = "gaussian")
    }, error = function(e) NULL)
    
    if (!is.null(cv_fit_enet)) {
      res$preds_strat1[i] <- predict(cv_fit_enet, newx = X_test_metab, s = "lambda.min")
      
      c_raw <- coef(cv_fit_enet, s = "lambda.min")
      valid_idx <- which(as.numeric(c_raw) != 0 & rownames(c_raw) != "(Intercept)")
      valid_names <- rownames(c_raw)[valid_idx]
      
      coef_vec <- rep(0, length(metab_cols))
      names(coef_vec) <- metab_cols
      coef_vec[valid_names] <- as.numeric(c_raw)[valid_idx]
      
      res$coefs_strat1[[i]] <- coef_vec
      res$selected_metabs_strat1[[i]] <- valid_names
    } else {
      res$preds_strat1[i] <- mean(y_train, na.rm=TRUE)
      res$coefs_strat1[[i]] <- setNames(rep(0, length(metab_cols)), metab_cols)
    }
    
    # --- 策略 2: Species 关联筛选 (修正版: 先收集P值，后做FDR) ---
    p_vals_s2 <- numeric(length(metab_cols))
    names(p_vals_s2) <- metab_cols
    
    for (nm in metab_cols) {
      fmla <- paste0(nm, " ~ ", species_col, " + ", paste(covar_cols, collapse = " + "))
      try({
        # 仅取交集列防止报错
        cols_needed <- c(nm, species_col, covar_cols)
        cols_valid <- intersect(cols_needed, names(train_data))
        m <- lm(as.formula(fmla), data = train_data[, cols_valid], na.action = na.omit)
        if (species_col %in% rownames(summary(m)$coefficients)) {
          p_vals_s2[nm] <- summary(m)$coefficients[species_col, "Pr(>|t|)"]
        } else {
          p_vals_s2[nm] <- NA
        }
      }, silent = TRUE)
    }
    
    sel_metabs_s2 <- character(0)
    valid_p_s2 <- p_vals_s2[!is.na(p_vals_s2)]
    if (length(valid_p_s2) > 0) {
      # 真正的 FDR 校正
      adj_p_s2 <- p.adjust(valid_p_s2, method = "fdr")
      sel_metabs_s2 <- names(adj_p_s2[adj_p_s2 < fdr_thresh])
    }
    res$selected_metabs_strat2[[i]] <- sel_metabs_s2
    
    # 策略 2 Ridge
    if (length(sel_metabs_s2) > 0) {
      fit_s2 <- tryCatch(cv.glmnet(as.matrix(train_data[, sel_metabs_s2, drop=F]), y_train, alpha=0, nfolds=10), error=function(e) NULL)
      if(!is.null(fit_s2)) {
        res$preds_strat2[i] <- predict(fit_s2, newx=as.matrix(test_data[, sel_metabs_s2, drop=F]), s="lambda.min")
        # 存系数逻辑同上(略简写)
        c_raw <- coef(fit_s2, s="lambda.min"); v_nms <- rownames(c_raw)[-1]
        c_vec <- setNames(rep(0, length(metab_cols)), metab_cols)
        c_vec[v_nms] <- as.numeric(c_raw)[-1]
        res$coefs_strat2[[i]] <- c_vec
      } else { res$preds_strat2[i] <- mean(y_train, na.rm=T) }
    } else { res$preds_strat2[i] <- mean(y_train, na.rm=T) }
    
    # --- 策略 3: Nut 关联筛选 (修正版: 保持逻辑，应用FDR < 0.2) ---
    p_vals_s3 <- numeric(length(metab_cols))
    names(p_vals_s3) <- metab_cols
    
    for (nm in metab_cols) {
      fmla <- paste0(nm, " ~ ", nut_col, " + ", paste(covar_cols, collapse = " + "))
      try({
        cols_needed <- c(nm, nut_col, covar_cols)
        cols_valid <- intersect(cols_needed, names(train_data))
        m <- lm(as.formula(fmla), data = train_data[, cols_valid], na.action = na.omit)
        if (nut_col %in% rownames(summary(m)$coefficients)) {
          p_vals_s3[nm] <- summary(m)$coefficients[nut_col, "Pr(>|t|)"]
        } else {
          p_vals_s3[nm] <- NA
        }
      }, silent = TRUE)
    }
    
    sel_metabs_s3 <- character(0)
    valid_p_s3 <- p_vals_s3[!is.na(p_vals_s3)]
    if (length(valid_p_s3) > 0) {
      adj_p_s3 <- p.adjust(valid_p_s3, method = "fdr")
      sel_metabs_s3 <- names(adj_p_s3[adj_p_s3 < fdr_thresh])
    }
    res$selected_metabs_strat3[[i]] <- sel_metabs_s3
    
    # 策略 3 Ridge
    if (length(sel_metabs_s3) > 0) {
      fit_s3 <- tryCatch(cv.glmnet(as.matrix(train_data[, sel_metabs_s3, drop=F]), y_train, alpha=0, nfolds=10), error=function(e) NULL)
      if(!is.null(fit_s3)) {
        res$preds_strat3[i] <- predict(fit_s3, newx=as.matrix(test_data[, sel_metabs_s3, drop=F]), s="lambda.min")
        c_raw <- coef(fit_s3, s="lambda.min"); v_nms <- rownames(c_raw)[-1]
        c_vec <- setNames(rep(0, length(metab_cols)), metab_cols)
        c_vec[v_nms] <- as.numeric(c_raw)[-1]
        res$coefs_strat3[[i]] <- c_vec
      } else { res$preds_strat3[i] <- mean(y_train, na.rm=T) }
    } else { res$preds_strat3[i] <- mean(y_train, na.rm=T) }
    
    # --- 策略 4: 并集 ---
    union_metabs <- unique(c(sel_metabs_s2, sel_metabs_s3))
    res$selected_metabs_strat4[[i]] <- union_metabs
    
    if (length(union_metabs) > 0) {
      fit_s4 <- tryCatch(cv.glmnet(as.matrix(train_data[, union_metabs, drop=F]), y_train, alpha=0, nfolds=10), error=function(e) NULL)
      if(!is.null(fit_s4)) {
        res$preds_strat4[i] <- predict(fit_s4, newx=as.matrix(test_data[, union_metabs, drop=F]), s="lambda.min")
        c_raw <- coef(fit_s4, s="lambda.min"); v_nms <- rownames(c_raw)[-1]
        c_vec <- setNames(rep(0, length(metab_cols)), metab_cols)
        c_vec[v_nms] <- as.numeric(c_raw)[-1]
        res$coefs_strat4[[i]] <- c_vec
      } else { res$preds_strat4[i] <- mean(y_train, na.rm=T) }
    } else { res$preds_strat4[i] <- mean(y_train, na.rm=T) }
    
  } # End LOOCV loop
  
  return(res)
}
tasks <- list(
  list(sp = "pred_an", nut = "all_nut_r_log_zscore"),
  list(sp = "pred_wn", nut = "walnut_r_log_zscore"),
  list(sp = "pred_tn", nut = "treenut_r_log_zscore")
)
out_dir <- "LOOCV_Results_FDR02"
if(!dir.exists(out_dir)) dir.create(out_dir)

for (k in 1:length(tasks)) {
  task <- tasks[[k]]
  
  file_name <- paste0(out_dir, "/res_", task$sp, "_", gsub("_r_log_zscore", "", task$nut), ".rds")
  
  if (file.exists(file_name)) {
    cat(paste0("跳过任务 ", k, ": ", file_name, " (文件已存在)\n"))
    next
  }
  
  tryCatch({
    
    result_list <- run_loocv_analysis(
      data = ToBeUsed.Taxon,
      species_col = task$sp,
      nut_col = task$nut,
      metab_cols = metabolite_cols,
      covar_cols = covariate_cols_species,
      fdr_thresh = 0.15 
    )
    
    # 立即保存结果到硬盘
    saveRDS(result_list, file = file_name)
    cat(paste0(">>> Success ", k, " is saved in ", file_name, "\n"))
    
  }, error = function(e) {
    # 如果报错，打印错误但不停止脚本，继续下一个任务
    cat(paste0("!!! ERROR:  ", k, " (", task$sp, ") FAILED !!!\n"))
    cat("ERROR INFO: ", e$message, "\n")
  })
  
  gc()
  cat("\n------------------------------------------\n")
}
summary_df <- data.frame()

for (k in 1:length(tasks)) {
  task <- tasks[[k]]
  
  nut_short <- gsub("_r_log_zscore", "", task$nut)
  file_name <- paste0(out_dir, "/res_", task$sp, "_", nut_short, ".rds")
  
  if (file.exists(file_name)) {
    res <- readRDS(file_name)
    
    y_true <- ToBeUsed.Taxon[[task$nut]]
    
    r_s1 <- cor(res$preds_strat1, y_true, use = "complete.obs") # ENet Only
    r_s2 <- cor(res$preds_strat2, y_true, use = "complete.obs") # Species Screen
    r_s3 <- cor(res$preds_strat3, y_true, use = "complete.obs") # Nut Screen
    r_s4 <- cor(res$preds_strat4, y_true, use = "complete.obs") # Union (主要关注这个)
    tmp <- data.frame(
      Species = task$sp,
      Nutrient = task$nut,
      R_Strat1_ENet = round(r_s1, 4),
      R_Strat2_SpScr = round(r_s2, 4),
      R_Strat3_NutScr = round(r_s3, 4),
      R_Strat4_Union = round(r_s4, 4)
    )
    summary_df <- rbind(summary_df, tmp)
    
  } else {
    cat(paste0("文件不存在: ", file_name, "\n"))
  }
}

# 4. 打印最终结果表
print(summary_df)


peanut_nut_col <- "peanut_r_log_zscore"
peanut_out_file <- paste0(out_dir, "/udd/n2gji/micro/data/res_no_species_peanut.rds")

if (file.exists(peanut_out_file)) {
  cat("跳过 Peanut 任务: 文件已存在\n")
} else {
  
  # 初始化容器
  n_samples <- nrow(ToBeUsed.Taxon)
  res_peanut <- list(
    preds_strat1 = numeric(n_samples), # ENet all
    preds_strat3 = numeric(n_samples), # Nut screen
    
    coefs_strat1 = vector("list", n_samples),
    coefs_strat3 = vector("list", n_samples),
    
    selected_metabs_strat1 = vector("list", n_samples),
    selected_metabs_strat3 = vector("list", n_samples)
  )
  
  # Peanut 专用 LOOCV 循环
  for (i in 1:n_samples) {
    if (i %% 50 == 0) cat(paste("  [Peanut] 处理折叠", i, "/", n_samples, "\n"))
    
    train_data <- ToBeUsed.Taxon[-i, ]
    test_data  <- ToBeUsed.Taxon[i, , drop = FALSE]
    
    X_train_metab <- as.matrix(train_data[, metabolite_cols])
    X_test_metab  <- as.matrix(test_data[, metabolite_cols])
    y_train       <- train_data[[peanut_nut_col]]
    
    cv_fit_enet <- tryCatch({
      cv.glmnet(X_train_metab, y_train, alpha = 0.5, nfolds = 10, family = "gaussian")
    }, error = function(e) NULL)
    
    if (!is.null(cv_fit_enet)) {
      res_peanut$preds_strat1[i] <- predict(cv_fit_enet, newx = X_test_metab, s = "lambda.min")
      
      c_raw <- coef(cv_fit_enet, s = "lambda.min")
      valid_idx <- which(as.numeric(c_raw) != 0 & rownames(c_raw) != "(Intercept)")
      valid_names <- rownames(c_raw)[valid_idx]
      
      coef_vec <- rep(0, length(metabolite_cols))
      names(coef_vec) <- metabolite_cols
      coef_vec[valid_names] <- as.numeric(c_raw)[valid_idx]
      
      res_peanut$coefs_strat1[[i]] <- coef_vec
      res_peanut$selected_metabs_strat1[[i]] <- valid_names
    } else {
      res_peanut$preds_strat1[i] <- mean(y_train, na.rm=TRUE)
    }
    p_vals_s3 <- numeric(length(metabolite_cols))
    names(p_vals_s3) <- metabolite_cols
    
    for (nm in metabolite_cols) {
      fmla <- paste0(nm, " ~ ", peanut_nut_col, " + ", paste(covariate_cols_species, collapse = " + "))
      try({
        cols_needed <- c(nm, peanut_nut_col, covariate_cols_species)
        cols_valid <- intersect(cols_needed, names(train_data))
        m <- lm(as.formula(fmla), data = train_data[, cols_valid], na.action = na.omit)
        if (peanut_nut_col %in% rownames(summary(m)$coefficients)) {
          p_vals_s3[nm] <- summary(m)$coefficients[peanut_nut_col, "Pr(>|t|)"]
        } else {
          p_vals_s3[nm] <- NA
        }
      }, silent = TRUE)
    }
    
    sel_metabs_s3 <- character(0)
    valid_p_s3 <- p_vals_s3[!is.na(p_vals_s3)]
    if (length(valid_p_s3) > 0) {
      adj_p_s3 <- p.adjust(valid_p_s3, method = "fdr")
      sel_metabs_s3 <- names(adj_p_s3[adj_p_s3 < 0.15]) # 保持与你主函数一致的 0.15
    }
    res_peanut$selected_metabs_strat3[[i]] <- sel_metabs_s3
    
    # 策略 3 Ridge 回归
    if (length(sel_metabs_s3) > 0) {
      fit_s3 <- tryCatch(cv.glmnet(as.matrix(train_data[, sel_metabs_s3, drop=F]), y_train, alpha=0, nfolds=10), error=function(e) NULL)
      if(!is.null(fit_s3)) {
        res_peanut$preds_strat3[i] <- predict(fit_s3, newx=as.matrix(test_data[, sel_metabs_s3, drop=F]), s="lambda.min")
        
        c_raw <- coef(fit_s3, s="lambda.min"); v_nms <- rownames(c_raw)[-1]
        c_vec <- setNames(rep(0, length(metabolite_cols)), metabolite_cols)
        c_vec[v_nms] <- as.numeric(c_raw)[-1]
        res_peanut$coefs_strat3[[i]] <- c_vec
      } else { res_peanut$preds_strat3[i] <- mean(y_train, na.rm=T) }
    } else { res_peanut$preds_strat3[i] <- mean(y_train, na.rm=T) }
    
  } # End Peanut Loop
  
  saveRDS(res_peanut, file = peanut_out_file)
  cat(paste0(">>> Peanut 分析完成，结果保存至: ", peanut_out_file, "\n"))
}

summary_df <- data.frame()


for (k in 1:length(tasks)) {
  task <- tasks[[k]]
  nut_short <- gsub("_r_log_zscore", "", task$nut)
  file_name <- paste0(out_dir, "/res_", task$sp, "_", nut_short, ".rds")
  
  if (file.exists(file_name)) {
    res <- readRDS(file_name)
    y_true <- ToBeUsed.Taxon[[task$nut]]
    
    tmp <- data.frame(
      Species = task$sp,
      Nutrient = task$nut,
      R_Strat1_ENet = round(cor(res$preds_strat1, y_true, use = "complete.obs"), 4),
      R_Strat2_SpScr = round(cor(res$preds_strat2, y_true, use = "complete.obs"), 4),
      R_Strat3_NutScr = round(cor(res$preds_strat3, y_true, use = "complete.obs"), 4),
      R_Strat4_Union = round(cor(res$preds_strat4, y_true, use = "complete.obs"), 4)
    )
    summary_df <- rbind(summary_df, tmp)
  }
}

# B. 处理 Peanut 任务
peanut_file <- paste0(out_dir, "/res_no_species_peanut.rds")
if (file.exists(peanut_file)) {
  res_pn <- readRDS(peanut_file)
  y_true_pn <- ToBeUsed.Taxon[["peanut_r_log_zscore"]]
  
  tmp_pn <- data.frame(
    Species = "NA (No pred_pn)",
    Nutrient = "peanut_r_log_zscore",
    R_Strat1_ENet = round(cor(res_pn$preds_strat1, y_true_pn, use = "complete.obs"), 4),
    R_Strat2_SpScr = NA, # 不适用
    R_Strat3_NutScr = round(cor(res_pn$preds_strat3, y_true_pn, use = "complete.obs"), 4),
    R_Strat4_Union = NA  # 不适用
  )
  summary_df <- rbind(summary_df, tmp_pn)
}

print(summary_df)


task_sp <- "pred_an"
task_nut <- "all_nut_r_log_zscore"
nut_col <- task_nut # 你后续代码用了这个变量名
min_selection_proportion<-0.8
rds_file <- paste0("/udd/n2gji/micro/data/LOOCV_Results_FDR02/res_", task_sp, "_", gsub("_r_log_zscore", "", task_nut), ".rds")

if(file.exists(rds_file)) {
  cat("正在加载结果文件:", rds_file, "\n")
  res_obj <- readRDS(rds_file)
  
  # === 关键：将 list 中的结果提取出来赋值给你的变量 ===
  
  # 策略 1
  predictions_strat1_nut <- res_obj$preds_strat1
  all_coefs_strat1       <- res_obj$coefs_strat1
  selected_metabs_per_fold_strat1 <- res_obj$selected_metabs_strat1
  
  # 策略 2
  predictions_strat2_nut <- res_obj$preds_strat2
  all_coefs_strat2       <- res_obj$coefs_strat2
  selected_metabs_per_fold_strat2 <- res_obj$selected_metabs_strat2
  
  # 策略 3
  predictions_strat3_nut <- res_obj$preds_strat3
  all_coefs_strat3       <- res_obj$coefs_strat3
  selected_metabs_per_fold_strat3 <- res_obj$selected_metabs_strat3
  
  # 策略 4
  predictions_strat4_nut <- res_obj$preds_strat4
  all_coefs_strat4       <- res_obj$coefs_strat4
  selected_metabs_per_fold_strat4 <- res_obj$selected_metabs_strat4
  
  cat("数据加载并解包完成！可以运行后续评分代码了。\n")
  
} else {
  stop(paste("找不到文件:", rds_file))
}
final_results_summary <- list()
num_folds <- nrow(ToBeUsed.Taxon) 

calculate_metrics <- function(actual, predicted) {
  rmse <- sqrt(mean((actual - predicted)^2, na.rm = TRUE))
  rss <- sum((actual - predicted)^2, na.rm = TRUE)
  tss <- sum((actual - mean(actual, na.rm = TRUE))^2, na.rm = TRUE)
  r2 <- 1 - (rss / tss)
  return(list(RMSE = rmse, R2 = r2))
}

process_strategy_results <- function(strategy_name, predictions, all_coefs, selected_metabs_per_fold, actual_nut_values, total_metabolite_cols, num_folds, min_prop) {
  metrics <- calculate_metrics(actual_nut_values, predictions)
  
  if (length(all_coefs) > 0 && !is.null(all_coefs[[1]])) {
    raw_avg_coefs <- colMeans(do.call(rbind, all_coefs), na.rm = TRUE)
  } else {
    raw_avg_coefs <- rep(0, length(total_metabolite_cols))
    names(raw_avg_coefs) <- total_metabolite_cols
  }
  
  all_selected_metabs_flat <- unlist(selected_metabs_per_fold)
  if (length(all_selected_metabs_flat) > 0) {
    selection_counts <- table(all_selected_metabs_flat)
    selection_frequencies <- selection_counts / num_folds
    selection_frequencies_df <- as.data.frame(selection_frequencies)
    colnames(selection_frequencies_df) <- c("Metabolite", "Frequency")
    selection_frequencies_df <- selection_frequencies_df[order(selection_frequencies_df$Frequency, decreasing = TRUE), ]
  } else {
    selection_frequencies_df <- data.frame(Metabolite = character(0), Frequency = numeric(0))
  }
  
  filtered_avg_coefs <- raw_avg_coefs
  if (nrow(selection_frequencies_df) > 0) {
    # 注意：这里使用传入的 min_prop 参数
    metabs_to_zero_out <- selection_frequencies_df$Metabolite[selection_frequencies_df$Frequency < min_prop]
    filtered_avg_coefs[metabs_to_zero_out] <- 0
  }
  
  missing_metabs <- setdiff(total_metabolite_cols, names(filtered_avg_coefs))
  if (length(missing_metabs) > 0) {
    temp_coefs <- rep(0, length(total_metabolite_cols))
    names(temp_coefs) <- total_metabolite_cols
    temp_coefs[names(filtered_avg_coefs)] <- filtered_avg_coefs
    filtered_avg_coefs <- temp_coefs
  }
  
  return(list(
    Description = paste0(strategy_name, ": "),
    RMSE = metrics$RMSE,
    R2 = metrics$R2,
    Raw_Average_Metabolite_Coefficients = raw_avg_coefs,
    Selection_Frequencies = selection_frequencies_df,
    Filtered_Average_Metabolite_Coefficients = filtered_avg_coefs
  ))
}

# Strategy 1
final_results_summary$Strategy1 <- process_strategy_results(
  "策略 1: ENet", predictions_strat1_nut, all_coefs_strat1, selected_metabs_per_fold_strat1, 
  ToBeUsed.Taxon[[nut_col]], metabolite_cols, num_folds, min_selection_proportion
)

# Strategy 2
final_results_summary$Strategy2 <- process_strategy_results(
  "策略 2: Species-Screen", predictions_strat2_nut, all_coefs_strat2, selected_metabs_per_fold_strat2, 
  ToBeUsed.Taxon[[nut_col]], metabolite_cols, num_folds, min_selection_proportion
)

# Strategy 3 (修正了你的原始代码中的变量名错误)
final_results_summary$Strategy3 <- process_strategy_results(
  "策略 3: Nut-Screen", predictions_strat3_nut, all_coefs_strat3, selected_metabs_per_fold_strat3, 
  ToBeUsed.Taxon[[nut_col]], metabolite_cols, num_folds, min_selection_proportion
)

# Strategy 4
final_results_summary$Strategy4 <- process_strategy_results(
  "策略 4: Union", predictions_strat4_nut, all_coefs_strat4, selected_metabs_per_fold_strat4, 
  ToBeUsed.Taxon[[nut_col]], metabolite_cols, num_folds, min_selection_proportion
)
table(final_results_summary$Strategy4$Selection_Frequencies)

# -------------------------------------------------------------------------------
# 3. 结果保存与筛选 (生成 merged_strategyX_data)
# -------------------------------------------------------------------------------
save(final_results_summary, file = "met_score_an_select.RData")

generate_merged_data <- function(strat_res, strategy_num, f_lvs_data) {
  freq_df <- strat_res$Selection_Frequencies
  coef_vec <- strat_res$Raw_Average_Metabolite_Coefficients
  
  coef_df <- data.frame(
    Metabolite = names(coef_vec),
    Coefficient = as.numeric(coef_vec),
    stringsAsFactors = FALSE
  )
  
  merged <- left_join(coef_df, freq_df, by = "Metabolite") %>%
    arrange(desc(Frequency), desc(abs(Coefficient))) %>%
    subset(Frequency >= 0.8) # 硬编码筛选 0.8
  
  # 合并注释信息
  if(exists("f_lvs_data") && nrow(merged) > 0) {
    f_lvs_match <- f_lvs_data %>% 
      mutate(met = rownames(.)) %>% 
      dplyr::select(any_of(c("met","method","metabolite_name","biochemical_name",
                             "class_metabolon","sub_class_metabolon","super_class_metabolon")))
    merged <- merge(merged, f_lvs_match, by.x = 'Metabolite', by.y = 'met', all = FALSE, sort = TRUE)
  }
  return(merged)
}

merged_strategy1_data <- generate_merged_data(final_results_summary$Strategy1, 1, f_lvs)
merged_strategy2_data <- generate_merged_data(final_results_summary$Strategy2, 2, f_lvs)
merged_strategy3_data <- generate_merged_data(final_results_summary$Strategy3, 3, f_lvs)

merged_strategy4_data <- generate_merged_data(final_results_summary$Strategy4, 4, f_lvs)
valid_mets <- unique(c(merged_strategy2_data$Metabolite, merged_strategy3_data$Metabolite))
merged_strategy4_data <- merged_strategy4_data %>% filter(Metabolite %in% valid_mets)

merged_data_dfs <- list(
  "Strategy1" = merged_strategy1_data,
  "Strategy2" = merged_strategy2_data,
  "Strategy3" = merged_strategy3_data,
  "Strategy4" = merged_strategy4_data
)
save(merged_data_dfs, file = "met_score_an.RData")

for (strat_key in names(merged_data_dfs)) {
  cat(paste0("计算 ", strat_key, " 的代谢物评分...\n"))
  current_merged_df <- merged_data_dfs[[strat_key]]
  
  # 初始化全零系数
  score_coefs <- rep(0, length(metabolite_cols))
  names(score_coefs) <- metabolite_cols
  
  # 填充筛选后的系数
  valid_metabs <- intersect(current_merged_df$Metabolite, metabolite_cols)
  if (length(valid_metabs) > 0) {
    coefs_to_use <- current_merged_df %>% 
      filter(Metabolite %in% valid_metabs) %>% 
      dplyr::select(Metabolite, Coefficient)
    score_coefs[coefs_to_use$Metabolite] <- coefs_to_use$Coefficient
  }
  
  # 矩阵乘法计算得分
  X_metabolites <- as.matrix(ToBeUsed.Taxon[, names(score_coefs)])
  B_coefs <- as.matrix(score_coefs, ncol = 1)
  metabolite_scores <- X_metabolites %*% B_coefs
  
  # 添加新列 (strategy1_met_score_an)
  score_col_name <- paste0(tolower(strat_key), "_met_score_an")
  ToBeUsed.Taxon[[score_col_name]] <- as.numeric(metabolite_scores) 
}

cat("评分计算完成！\n")

# -------------------------------------------------------------------------------
# 5. 验证结果
# -------------------------------------------------------------------------------
cat("\n真实营养素与各策略评分的相关性:\n")
print(cor(ToBeUsed.Taxon[[nut_col]], ToBeUsed.Taxon$strategy1_met_score_an, use="complete.obs"))
print(cor(ToBeUsed.Taxon[[nut_col]], ToBeUsed.Taxon$strategy2_met_score_an, use="complete.obs"))
print(cor(ToBeUsed.Taxon[[nut_col]], ToBeUsed.Taxon$strategy3_met_score_an, use="complete.obs"))
print(cor(ToBeUsed.Taxon[[nut_col]], ToBeUsed.Taxon$strategy4_met_score_an, use="complete.obs"))


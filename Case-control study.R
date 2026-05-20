#################NHS t2D case-control study#############
library(mediation)
library(naniar)
#----------Diabetes case-control---------#
diabetes_met <- merge_metab_data(endpoints = c("diabetes"),
                                 cohorts = "nhs1",
                                 transformation = "transform_ln_z_score",
                                 keep_unknown_metabolites = F,
                                 keep_failed_pm_metabolites = F, 
                                 impute_missing_function = "impute_one_half_min",
                                 impute_cutoff=0.3,
                                 combine_cohorts = TRUE)

dim(diabetes_met$expr_set$all_cohorts) # Features 351   Samples 1456
names(p_db)
p_db=pData(diabetes_met$expr_set$all_cohorts)
f_db=fData(diabetes_met$expr_set$all_cohorts)
e_db=data.frame(exprs(diabetes_met$expr_set$all_cohorts))

db_cv <- f_db %>% filter(mean_cv<30 & mean_icc>0.4) # 351 to 310 metabolites
e_db_cv <- subset(e_db, (rownames(e_db) %in% rownames(db_cv))) # 310 metabolites
a_db=data.frame(t(e_db_cv))
rownames(p_db) <- paste("X", rownames(p_db), sep="")
identical(rownames(a_db), rownames(p_db)) # should be TRUE
z_db=cbind(a_db,p_db) 
non_na_counts <- colMeans(is.na(z_db[,1:310]))
non_na_counts 
z_db_part1 <- z_db[,1:310][, non_na_counts <0.3] # 310 to 302 metabolites
z_db=cbind(z_db_part1, z_db[, 311:332])
names(z_db)
missing <- miss_case_summary(z_db[,1:302])#all 0
z_db$idsix<-as.numeric(substr(z_db$id, 1, 6))
z_db$id<-z_db$idsix
#names(z_db)
#write.csv(z_db,paste0(path,'db_cc_nhs.csv'),row.names = FALSE, na='.')
#z_db<-read.csv('/udd/n2gji/micro/db_cc_nhs.csv')
# Start here
nhs_db_cohort<-read.csv('/udd/n2gji/data/nhs1_db_cc.csv')
nhs_db_cohort<- nhs_db_cohort %>%
  group_by(id) %>%
  filter(row_number() == 1) %>%
  ungroup()

nhs_db_cohort<-nhs_db_cohort%>%
  mutate(pnut_all=rowSums(cbind(pnuta,pbuta),na.rm=F))

# 150+ missing value of nut,nearly 12%, I don't know wether it is correct to impute#
#Impute with median#
nhs_db_cohort$nuta<-ifelse(is.na(nhs_db_cohort$nuta),median(nhs_db_cohort$nuta,na.rm = T),nhs_db_cohort$nuta)
nhs_db_cohort$onuta<-ifelse(is.na(nhs_db_cohort$onuta),median(nhs_db_cohort$onuta,na.rm = T),nhs_db_cohort$onuta)
nhs_db_cohort$pnut_all<-ifelse(is.na(nhs_db_cohort$pnut_all),median(nhs_db_cohort$pnut_all,na.rm = T),nhs_db_cohort$pnut_all)
summary(nhs_db_cohort$fish)
names(z_db)
z_db$idsix<-as.numeric(substr(z_db$id, 1, 6))
z_db$id<-z_db$idsix
nhs_db_cohort<-merge(nhs_db_cohort,z_db[,c(303,307,309:310)],by="id",all.x = T)

nhs_db_cohort$db<-ifelse(nhs_db_cohort$caco=="case",1,0)
table(nhs_db_cohort$db)
#load("/udd/n2gji/micro/nut_final/t2d/data/nhs_db_cc.RData")
#############################Calculate score###########################
metabolite_cols <- grep("^HMDB", colnames(nhs_db_cohort), value = TRUE)
#-----------------------T2D case-control---------------#
for (strat_key in names(merged_data_dfs_tn_eln)) {
  cat(paste0("计算 ", strat_key, " 的代谢物评分...\n"))
  
  # 获取当前策略的 merged_data 数据框
  current_merged_df <- merged_data_dfs_tn_eln[[strat_key]]
  
  # 1. 准备 score_coefs 向量
  # 初始化一个包含所有代谢物（来自 metabolite_cols）的系数向量，所有系数都为 0。
  score_coefs <- rep(0, length(metabolite_cols))
  names(score_coefs) <- metabolite_cols
  
  # 从 current_merged_df 中提取代谢物名称和它们的 Coefficient
  # 并只考虑那些在总 metabolite_cols 中存在的代谢物
  metabs_in_merged_and_full_list <- intersect(current_merged_df$Metabolite, metabolite_cols)
  
  if (length(metabs_in_merged_and_full_list) > 0) {
    coefficients_from_merged <- current_merged_df %>%
      filter(Metabolite %in% metabs_in_merged_and_full_list) %>%
      dplyr::select(Metabolite, Coefficient)
    
    # 将从 merged_data 中获取的系数填充到 score_coefs 向量中
    # 这样，未在 merged_data 中列出但在 metabolite_cols 中的代谢物，其系数将保持为 0
    score_coefs[coefficients_from_merged$Metabolite] <- coefficients_from_merged$Coefficient
  } else {
    cat(paste0("  注意：策略 ", strat_key, " 的 merged_data 中没有找到与总代谢物列表匹配的代谢物。该策略在 nhs_db_cohort 上所有评分将为0。\n"))
  }
  model_expected_metabolites <- names(score_coefs) 
  
  X_cohort_metabolites <- matrix(0, nrow = nrow(nhs_db_cohort), ncol = length(model_expected_metabolites))
  colnames(X_cohort_metabolites) <- model_expected_metabolites
  
  # 识别 nhs_db_cohort 中实际存在的模型期望代谢物
  present_metabs_in_current_cohort <- intersect(model_expected_metabolites, names(nhs_db_cohort))
  
  # 将 nhs_db_cohort 中存在的代谢物数据复制到 X_cohort_metabolites 中
  if (length(present_metabs_in_current_cohort) > 0) {
    # 确保提取的数据是数值型，并处理潜在的NA
    temp_data <- as.matrix(nhs_db_cohort[, present_metabs_in_current_cohort])
    temp_data[is.na(temp_data)] <- 0 # 将 nhs_db_cohort 中代谢物的NA值也视为0
    X_cohort_metabolites[, present_metabs_in_current_cohort] <- temp_data
  }
  
  # 3. 将系数向量转换为列矩阵，以便进行矩阵乘法
  B_coefs <- as.matrix(score_coefs, ncol = 1)
  
  # 4. 计算代谢物评分
  metabolite_scores <- X_cohort_metabolites %*% B_coefs
  
  # 5. 将计算出的评分添加为 nhs_db_cohort 的新列
  score_col_name <- paste0("eln_tn_score_", tolower(strat_key), "")
  nhs_db_cohort[[score_col_name]] <- as.numeric(metabolite_scores) 
}
for (strat_key in names(merged_data_dfs_wn_eln)) {
  cat(paste0("计算 ", strat_key, " 的代谢物评分...\n"))
  
  # 获取当前策略的 merged_data 数据框
  current_merged_df <- merged_data_dfs_wn_eln[[strat_key]]
  
  # 1. 准备 score_coefs 向量
  # 初始化一个包含所有代谢物（来自 metabolite_cols）的系数向量，所有系数都为 0。
  score_coefs <- rep(0, length(metabolite_cols))
  names(score_coefs) <- metabolite_cols
  
  # 从 current_merged_df 中提取代谢物名称和它们的 Coefficient
  # 并只考虑那些在总 metabolite_cols 中存在的代谢物
  metabs_in_merged_and_full_list <- intersect(current_merged_df$Metabolite, metabolite_cols)
  
  if (length(metabs_in_merged_and_full_list) > 0) {
    coefficients_from_merged <- current_merged_df %>%
      filter(Metabolite %in% metabs_in_merged_and_full_list) %>%
      dplyr :: select(Metabolite, Coefficient)
    
    # 将从 merged_data 中获取的系数填充到 score_coefs 向量中
    # 这样，未在 merged_data 中列出但在 metabolite_cols 中的代谢物，其系数将保持为 0
    score_coefs[coefficients_from_merged$Metabolite] <- coefficients_from_merged$Coefficient
  } else {
    cat(paste0("  注意：策略 ", strat_key, " 的 merged_data 中没有找到与总代谢物列表匹配的代谢物。该策略在 nhs_db_cohort 上所有评分将为0。\n"))
  }
  
  # 2. 准备 nhs_db_cohort 的代谢物矩阵 (X_cohort_metabolites)
  #    这个矩阵将包含所有 score_coefs 期望的代谢物列，即使 nhs_db_cohort 中没有。
  #    缺失的列将自动填充为0。
  
  # 获取模型期望的所有代谢物列名，并确保它们的顺序与 score_coefs 一致
  model_expected_metabolites <- names(score_coefs) 
  
  # 创建一个与模型期望列数一致的空矩阵，并用0填充
  X_cohort_metabolites <- matrix(0, nrow = nrow(nhs_db_cohort), ncol = length(model_expected_metabolites))
  colnames(X_cohort_metabolites) <- model_expected_metabolites
  
  # 识别 nhs_db_cohort 中实际存在的模型期望代谢物
  present_metabs_in_current_cohort <- intersect(model_expected_metabolites, names(nhs_db_cohort))
  
  # 将 nhs_db_cohort 中存在的代谢物数据复制到 X_cohort_metabolites 中
  if (length(present_metabs_in_current_cohort) > 0) {
    # 确保提取的数据是数值型，并处理潜在的NA
    temp_data <- as.matrix(nhs_db_cohort[, present_metabs_in_current_cohort])
    temp_data[is.na(temp_data)] <- 0 # 将 nhs_db_cohort 中代谢物的NA值也视为0
    X_cohort_metabolites[, present_metabs_in_current_cohort] <- temp_data
  }
  
  # 3. 将系数向量转换为列矩阵，以便进行矩阵乘法
  B_coefs <- as.matrix(score_coefs, ncol = 1)
  
  # 4. 计算代谢物评分
  metabolite_scores <- X_cohort_metabolites %*% B_coefs
  
  # 5. 将计算出的评分添加为 nhs_db_cohort 的新列
  score_col_name <- paste0("eln_wn_score_", tolower(strat_key), "")
  nhs_db_cohort[[score_col_name]] <- as.numeric(metabolite_scores) 
}
for (strat_key in names(merged_data_dfs_an_eln)) {
  cat(paste0("计算 ", strat_key, " 的代谢物评分...\n"))
  
  # 获取当前策略的 merged_data 数据框
  current_merged_df <- merged_data_dfs_an_eln[[strat_key]]
  
  # 1. 准备 score_coefs 向量
  # 初始化一个包含所有代谢物（来自 metabolite_cols）的系数向量，所有系数都为 0。
  score_coefs <- rep(0, length(metabolite_cols))
  names(score_coefs) <- metabolite_cols
  
  # 从 current_merged_df 中提取代谢物名称和它们的 Coefficient
  # 并只考虑那些在总 metabolite_cols 中存在的代谢物
  metabs_in_merged_and_full_list <- intersect(current_merged_df$Metabolite, metabolite_cols)
  
  if (length(metabs_in_merged_and_full_list) > 0) {
    coefficients_from_merged <- current_merged_df %>%
      filter(Metabolite %in% metabs_in_merged_and_full_list) %>%
      dplyr:: select(Metabolite, Coefficient)
    
    # 将从 merged_data 中获取的系数填充到 score_coefs 向量中
    # 这样，未在 merged_data 中列出但在 metabolite_cols 中的代谢物，其系数将保持为 0
    score_coefs[coefficients_from_merged$Metabolite] <- coefficients_from_merged$Coefficient
  } else {
    cat(paste0("  注意：策略 ", strat_key, " 的 merged_data 中没有找到与总代谢物列表匹配的代谢物。该策略在 nhs_db_cohort 上所有评分将为0。\n"))
  }
  
  # 2. 准备 nhs_db_cohort 的代谢物矩阵 (X_cohort_metabolites)
  #    这个矩阵将包含所有 score_coefs 期望的代谢物列，即使 nhs_db_cohort 中没有。
  #    缺失的列将自动填充为0。
  
  # 获取模型期望的所有代谢物列名，并确保它们的顺序与 score_coefs 一致
  model_expected_metabolites <- names(score_coefs) 
  
  # 创建一个与模型期望列数一致的空矩阵，并用0填充
  X_cohort_metabolites <- matrix(0, nrow = nrow(nhs_db_cohort), ncol = length(model_expected_metabolites))
  colnames(X_cohort_metabolites) <- model_expected_metabolites
  
  # 识别 nhs_db_cohort 中实际存在的模型期望代谢物
  present_metabs_in_current_cohort <- intersect(model_expected_metabolites, names(nhs_db_cohort))
  
  # 将 nhs_db_cohort 中存在的代谢物数据复制到 X_cohort_metabolites 中
  if (length(present_metabs_in_current_cohort) > 0) {
    # 确保提取的数据是数值型，并处理潜在的NA
    temp_data <- as.matrix(nhs_db_cohort[, present_metabs_in_current_cohort])
    temp_data[is.na(temp_data)] <- 0 # 将 nhs_db_cohort 中代谢物的NA值也视为0
    X_cohort_metabolites[, present_metabs_in_current_cohort] <- temp_data
  }
  
  # 3. 将系数向量转换为列矩阵，以便进行矩阵乘法
  B_coefs <- as.matrix(score_coefs, ncol = 1)
  
  # 4. 计算代谢物评分
  metabolite_scores <- X_cohort_metabolites %*% B_coefs
  
  # 5. 将计算出的评分添加为 nhs_db_cohort 的新列
  score_col_name <- paste0("eln_an_score_", tolower(strat_key), "")
  nhs_db_cohort[[score_col_name]] <- as.numeric(metabolite_scores) 
}
for (strat_key in names(merged_data_dfs_pn_eln)) {
  cat(paste0("计算 ", strat_key, " 的代谢物评分...\n"))
  
  # 获取当前策略的 merged_data 数据框
  current_merged_df <- merged_data_dfs_pn_eln[[strat_key]]
  
  # 1. 准备 score_coefs 向量
  # 初始化一个包含所有代谢物（来自 metabolite_cols）的系数向量，所有系数都为 0。
  score_coefs <- rep(0, length(metabolite_cols))
  names(score_coefs) <- metabolite_cols
  
  # 从 current_merged_df 中提取代谢物名称和它们的 Coefficient
  # 并只考虑那些在总 metabolite_cols 中存在的代谢物
  metabs_in_merged_and_full_list <- intersect(current_merged_df$Metabolite, metabolite_cols)
  
  if (length(metabs_in_merged_and_full_list) > 0) {
    coefficients_from_merged <- current_merged_df %>%
      filter(Metabolite %in% metabs_in_merged_and_full_list) %>%
      dplyr:: select(Metabolite, Coefficient)
    
    # 将从 merged_data 中获取的系数填充到 score_coefs 向量中
    # 这样，未在 merged_data 中列出但在 metabolite_cols 中的代谢物，其系数将保持为 0
    score_coefs[coefficients_from_merged$Metabolite] <- coefficients_from_merged$Coefficient
  } else {
    cat(paste0("  注意：策略 ", strat_key, " 的 merged_data 中没有找到与总代谢物列表匹配的代谢物。该策略在 nhs_db_cohort 上所有评分将为0。\n"))
  }
  
  # 2. 准备 nhs_db_cohort 的代谢物矩阵 (X_cohort_metabolites)
  #    这个矩阵将包含所有 score_coefs 期望的代谢物列，即使 nhs_db_cohort 中没有。
  #    缺失的列将自动填充为0。
  
  # 获取模型期望的所有代谢物列名，并确保它们的顺序与 score_coefs 一致
  model_expected_metabolites <- names(score_coefs) 
  
  # 创建一个与模型期望列数一致的空矩阵，并用0填充
  X_cohort_metabolites <- matrix(0, nrow = nrow(nhs_db_cohort), ncol = length(model_expected_metabolites))
  colnames(X_cohort_metabolites) <- model_expected_metabolites
  
  # 识别 nhs_db_cohort 中实际存在的模型期望代谢物
  present_metabs_in_current_cohort <- intersect(model_expected_metabolites, names(nhs_db_cohort))
  
  # 将 nhs_db_cohort 中存在的代谢物数据复制到 X_cohort_metabolites 中
  if (length(present_metabs_in_current_cohort) > 0) {
    # 确保提取的数据是数值型，并处理潜在的NA
    temp_data <- as.matrix(nhs_db_cohort[, present_metabs_in_current_cohort])
    temp_data[is.na(temp_data)] <- 0 # 将 nhs_db_cohort 中代谢物的NA值也视为0
    X_cohort_metabolites[, present_metabs_in_current_cohort] <- temp_data
  }
  
  # 3. 将系数向量转换为列矩阵，以便进行矩阵乘法
  B_coefs <- as.matrix(score_coefs, ncol = 1)
  
  # 4. 计算代谢物评分
  metabolite_scores <- X_cohort_metabolites %*% B_coefs
  
  # 5. 将计算出的评分添加为 nhs_db_cohort 的新列
  score_col_name <- paste0("eln_pn_score_", tolower(strat_key), "")
  nhs_db_cohort[[score_col_name]] <- as.numeric(metabolite_scores) 
}
nhs_ahei <-read.csv("/udd/n2gji/data/nhs_ahei.csv")#Go to nut_0228 calclate nnut_ahei
nhs_ahei<-nhs_ahei %>% mutate(ahei_nnut86=rowSums(cbind(ahei2010_vegI86,ahei2010_frtI86,
                                                        ahei2010_ssbI86,ahei2010_whgrnI86,
                                                        ahei2010_nutI86,ahei2010_rmtI86,
                                                        ahei2010_etohI86,ahei2010_ptranI86,
                                                        ahei2010_omegaI86,ahei2010_polyI86,ahei2010_naI86),na.rm = F),
                              ahei_nnut90=rowSums(cbind(ahei2010_vegI90,ahei2010_frtI90,
                                                        ahei2010_ssbI90,ahei2010_whgrnI90,
                                                        ahei2010_nutI90,ahei2010_rmtI90,
                                                        ahei2010_etohI90,ahei2010_ptranI90,
                                                        ahei2010_omegaI90,ahei2010_polyI90,ahei2010_naI90),na.rm = F),
                              ahei_nnut94=rowSums(cbind(ahei2010_vegI94,ahei2010_frtI94,
                                                        ahei2010_ssbI94,ahei2010_whgrnI94,
                                                        ahei2010_nutI94,ahei2010_rmtI94,
                                                        ahei2010_etohI94,ahei2010_ptranI94,
                                                        ahei2010_omegaI94,ahei2010_polyI94,ahei2010_naI94),na.rm = F),
                              ahei_nnut98=rowSums(cbind(ahei2010_vegI98,ahei2010_frtI98,
                                                        ahei2010_ssbI98,ahei2010_whgrnI98,
                                                        ahei2010_nutI98,ahei2010_rmtI98,
                                                        ahei2010_etohI98,ahei2010_ptranI98,
                                                        ahei2010_omegaI98,ahei2010_polyI98,ahei2010_naI98),na.rm = F),
                              ahei_nnut02=rowSums(cbind(ahei2010_vegI02,ahei2010_frtI02,
                                                        ahei2010_ssbI02,ahei2010_whgrnI02,
                                                        ahei2010_nutI02,ahei2010_rmtI02,
                                                        ahei2010_etohI02,ahei2010_ptranI02,
                                                        ahei2010_omegaI02,ahei2010_polyI02,ahei2010_naI02),na.rm = F),
                              ahei_nnut06=rowSums(cbind(ahei2010_vegI06,ahei2010_frtI06,
                                                        ahei2010_ssbI06,ahei2010_whgrnI06,
                                                        ahei2010_nutI06,ahei2010_rmtI06,
                                                        ahei2010_etohI06,ahei2010_ptranI06,
                                                        ahei2010_omegaI06,ahei2010_polyI06,ahei2010_naI06),na.rm = F),
                              ahei_nnut10=rowSums(cbind(ahei2010_vegI10,ahei2010_frtI10,
                                                        ahei2010_ssbI10,ahei2010_whgrnI10,
                                                        ahei2010_nutI10,ahei2010_rmtI10,
                                                        ahei2010_etohI10,ahei2010_ptranI10,
                                                        ahei2010_omegaI10,ahei2010_polyI10,ahei2010_naI10),na.rm = F))
names(nhs_ahei)
nhs_db_cohort<-merge(nhs_db_cohort,nhs_ahei[,c(1,211)],by="id",all.x = T)
summary(nhs_db_cohort$ahei_nnut90)

nhs_db_cohort <-nhs_db_cohort %>% 
  #---------------------------------------------------------------------------
#---------------- missing values -------------------------------------------
mutate(pa_metbase=case_when(
  is.na(pa_metbase)~act,
  TRUE ~ pa_metbase
)) %>%  # still NA: 51
  mutate(act=case_when(
    is.na(act)~pa_metbase,
    TRUE ~ act
  )) %>% # still NA: 51
  mutate(pa_metbase=case_when(
    is.na(pa_metbase)~median(pa_metbase,na.rm=TRUE),
    TRUE ~ pa_metbase
  )) %>% 
  mutate(act=case_when(
    is.na(act)~median(act,na.rm=TRUE),
    TRUE ~ act
  )) %>% 
  mutate(aheibase=case_when(
    is.na(aheibase)~median(aheibase,na.rm=TRUE),
    TRUE ~ aheibase
  )) %>%  # No missing after carry on from bqx
  mutate(ahei_nnut90=case_when(
    is.na(ahei_nnut90)~median(ahei_nnut90,na.rm=TRUE),
    TRUE ~ ahei_nnut90
  )) %>%  # No missing after carry on from bqx
  #mutate(smoke=case_when(
  #  is.na(smoke)~ as.character(smk), #  smoke with NA also presents with NA in smk
  #  TRUE ~ smoke
  #))  %>%
  mutate(smoke=as.factor(smoke)) %>% 
  mutate(smoke=case_when(
    smoke=='1' ~'never',
    smoke=='2' ~'past',
    smoke=='3' ~'current',
    TRUE ~ smoke
  )) %>%
  # mutate(smoke=case_when(
  #   is.na(smoke)~ 'missing', # still 92 missing
  # TRUE ~ smoke
  # )) %>%
  mutate(bmi=case_when(
    is.na(bmi) ~ bmicontbase, # still 88 NAs
    TRUE ~ bmi 
  )) %>% 
  mutate(bmicontbase=case_when(
    is.na(bmicontbase) ~ bmi, # also 88 NAs
    TRUE ~ bmicontbase
  )) %>% 
  mutate(bmi=case_when(
    is.na(bmi) ~ median(bmi,na.rm=TRUE), # impute 88 missing with median
    TRUE ~ bmi
  )) %>% 
  mutate(bmicontbase=case_when(
    is.na(bmicontbase) ~ median(bmicontbase,na.rm=TRUE),
    TRUE ~ bmicontbase
  )) %>% 
  mutate(bmibase3cat=case_when(
    bmicontbase<25 ~ 1,
    25<=bmicontbase&bmicontbase<30 ~ 2,
    bmicontbase>=30 ~ 3
  )) %>% 
  #---------------------------------------------------------------------------
mutate(calorbaseq=ntile(calorbase,5)) %>% 
  mutate(pa_metbaseq=ntile(pa_metbase,5)) %>% 
  mutate(aheibaseq=ntile(aheibase,5)) %>% 
  mutate(aheinnutq=ntile(ahei_nnut90,5)) %>% 
  mutate(anveprorbaseq=ntile(anvepror,5)) %>% 
  mutate(alcoq=ntile(alco,5)) %>% 
  #---------------------------------------------------------------------------
mutate(sex=factor(sex,levels=c(0,1),labels=c('Female','Male'))) %>% 
  mutate(calorbaseq=factor(calorbaseq,levels=c(1:5))) %>% 
  mutate(smk=factor(smk,levels=c(1:3),labels=c('never','past','current'))) %>% 
  mutate(pa_metbaseq=factor(pa_metbaseq,levels=c(1:5))) %>%
  mutate(alcoq=factor(alcoq,levels=c(1:5))) %>%
  mutate(aheibaseq=factor(aheibaseq,levels=c(1:5))) %>% 
  mutate(aheinnutq=factor(aheinnutq,levels=c(1:5))) %>% 
  mutate(anveprorbaseq=factor(anveprorbaseq,levels=c(1:5))) %>% 
  mutate(bmibase3cat=factor(bmibase3cat,levels=c(1:3),labels=c('<25','<30','>=30'))) %>% 
  mutate(bmi3cat=factor(bmi3cat,levels=c(1:3),labels=c('<25','<30','>=30'))) 
nhs_db_cohort$alcoq[is.na(nhs_db_cohort$alcoq)]<-1
nhs_db_cohort$smoke[nhs_db_cohort$smoke!="never"&nhs_db_cohort$smoke!="past"&nhs_db_cohort$smoke!="current"]<-"never"
###Merge ahei without nut####
# 根据 db=0 的组计算五分位分箱标准
db0_quantiles <- nhs_db_cohort %>% 
  filter(db == 0) %>% 
  summarise(
    eln1_an_cutoff = quantile(eln_an_score_strategy1, probs = c(0, 0.20, 0.40, 0.60, 0.80, 1), na.rm = TRUE),
    # eln2_an_cutoff= quantile(eln_an_score_strategy2, probs = c(0, 0.20, 0.40, 0.60, 0.80, 1), na.rm = TRUE),
    eln3_an_cutoff = quantile(eln_an_score_strategy3, probs = c(0, 0.20, 0.40, 0.60, 0.80, 1), na.rm = TRUE),
    eln4_an_cutoff= quantile(eln_an_score_strategy4, probs = c(0, 0.20, 0.40, 0.60, 0.80, 1), na.rm = TRUE),
    eln1_tn_cutoff = quantile(eln_tn_score_strategy1, probs = c(0, 0.20, 0.40, 0.60, 0.80, 1), na.rm = TRUE),
    eln2_tn_cutoff= quantile(eln_tn_score_strategy2, probs = c(0, 0.20, 0.40, 0.60, 0.80, 1), na.rm = TRUE),
    eln3_tn_cutoff = quantile(eln_tn_score_strategy3, probs = c(0, 0.20, 0.40, 0.60, 0.80, 1), na.rm = TRUE),
    eln4_tn_cutoff= quantile(eln_tn_score_strategy4, probs = c(0, 0.20, 0.40, 0.60, 0.80, 1), na.rm = TRUE),
    eln1_wn_cutoff = quantile(eln_wn_score_strategy1, probs = c(0, 0.20, 0.40, 0.60, 0.80, 1), na.rm = TRUE),
    eln2_wn_cutoff= quantile(eln_wn_score_strategy2, probs = c(0, 0.20, 0.40, 0.60, 0.80, 1), na.rm = TRUE),
    eln3_wn_cutoff = quantile(eln_wn_score_strategy3, probs = c(0, 0.20, 0.40, 0.60, 0.80, 1), na.rm = TRUE),
    eln4_wn_cutoff= quantile(eln_wn_score_strategy4, probs = c(0, 0.20, 0.40, 0.60, 0.80, 1), na.rm = TRUE),
    eln1_pn_cutoff = quantile(eln_pn_score_strategy1, probs = c(0, 0.20, 0.40, 0.60, 0.80, 1), na.rm = TRUE),
    #eln2_pn_cutoff= quantile(eln_pn_score_strategy2, probs = c(0, 0.20, 0.40, 0.60, 0.80, 1), na.rm = TRUE),
    eln3_pn_cutoff = quantile(eln_pn_score_strategy3, probs = c(0, 0.20, 0.40, 0.60, 0.80, 1), na.rm = TRUE),
    #eln4_pn_cutoff= quantile(eln_pn_score_strategy4, probs = c(0, 0.20, 0.40, 0.60, 0.80, 1), na.rm = TRUE)
    
  )
View(db0_quantiles)
# 使用 db=0 的五分位标准进行全体数据分箱
nhs_db_cohort <- nhs_db_cohort %>%
  mutate(
    score1_an_elnq = cut(eln_an_score_strategy1,
                         breaks = unlist(db0_quantiles$eln1_an_cutoff),
                         labels = c("Q1", "Q2", "Q3", "Q4", "Q5"),
                         include.lowest = TRUE),
    # score2_an_elnq= cut(eln_an_score_strategy2,
    #                    breaks = unlist(db0_quantiles$eln2_an_cutoff),
    #                 labels = c("Q1", "Q2", "Q3", "Q4", "Q5"),
    #               include.lowest = TRUE),
    score3_an_elnq = cut(eln_an_score_strategy3,
                         breaks = unlist(db0_quantiles$eln3_an_cutoff),
                         labels = c("Q1", "Q2", "Q3", "Q4", "Q5"),
                         include.lowest = TRUE),
    score4_an_elnq= cut(eln_an_score_strategy4,
                        breaks = unlist(db0_quantiles$eln4_an_cutoff),
                        labels = c("Q1", "Q2", "Q3", "Q4", "Q5"),
                        include.lowest = TRUE),
    score1_tn_elnq = cut(eln_tn_score_strategy1,
                         breaks = unlist(db0_quantiles$eln1_tn_cutoff),
                         labels = c("Q1", "Q2", "Q3", "Q4", "Q5"),
                         include.lowest = TRUE),
    score2_tn_elnq= cut(eln_tn_score_strategy2,
                        breaks = unlist(db0_quantiles$eln2_tn_cutoff),
                        labels = c("Q1", "Q2", "Q3", "Q4", "Q5"),
                        include.lowest = TRUE),
    score3_tn_elnq = cut(eln_tn_score_strategy3,
                         breaks = unlist(db0_quantiles$eln3_tn_cutoff),
                         labels = c("Q1", "Q2", "Q3", "Q4", "Q5"),
                         include.lowest = TRUE),
    score4_tn_elnq= cut(eln_tn_score_strategy4,
                        breaks = unlist(db0_quantiles$eln4_tn_cutoff),
                        labels = c("Q1", "Q2", "Q3", "Q4", "Q5"),
                        include.lowest = TRUE),
    score1_wn_elnq = cut(eln_wn_score_strategy1,
                         breaks = unlist(db0_quantiles$eln1_wn_cutoff),
                         labels = c("Q1", "Q2", "Q3", "Q4", "Q5"),
                         include.lowest = TRUE),
    score2_wn_elnq= cut(eln_wn_score_strategy2,
                        breaks = unlist(db0_quantiles$eln2_wn_cutoff),
                        labels = c("Q1", "Q2", "Q3", "Q4", "Q5"),
                        include.lowest = TRUE),
    score3_wn_elnq = cut(eln_wn_score_strategy3,
                         breaks = unlist(db0_quantiles$eln3_wn_cutoff),
                         labels = c("Q1", "Q2", "Q3", "Q4", "Q5"),
                         include.lowest = TRUE),
    score4_wn_elnq= cut(eln_wn_score_strategy4,
                        breaks = unlist(db0_quantiles$eln4_wn_cutoff),
                        labels = c("Q1", "Q2", "Q3", "Q4", "Q5"),
                        include.lowest = TRUE),
    score1_pn_elnq = cut(eln_pn_score_strategy1,
                         breaks = unlist(db0_quantiles$eln1_pn_cutoff),
                         labels = c("Q1", "Q2", "Q3", "Q4", "Q5"),
                         include.lowest = TRUE),
    score3_pn_elnq = cut(eln_pn_score_strategy3,
                         breaks = unlist(db0_quantiles$eln3_pn_cutoff),
                         labels = c("Q1", "Q2", "Q3", "Q4", "Q5"),
                         include.lowest = TRUE),
  )
dim(nhs_db_cohort)
nhs_db_cohort %>%
  group_by(caco) %>%
  summarise(
    Q1     = quantile(ahei_nnut90, 0.25, na.rm = TRUE),
    median = median(ahei_nnut90, na.rm = TRUE),
    Q3     = quantile(ahei_nnut90, 0.75, na.rm = TRUE)
  )

#############################Start case-control###########################
run_clogit <- function(data, score_vars, covariates = NULL) {
  results <- data.frame()
  
  for (score_var in score_vars) {
    # 构造公式
    formula_str <- paste("db ~", score_var)
    if (!is.null(covariates)) {
      formula_str <- paste(formula_str, "+", paste(covariates, collapse = " + "))
    }
    formula_str <- paste(formula_str, "+ strata(matchid) ")
    formula <- as.formula(formula_str)
    
    # 拟合模型并提取结果
    model_summary <- clogit(formula, data = data) %>%
      broom::tidy() %>%
      filter(grepl(score_var, term)) %>%
      transmute(
        Variable = score_var,
        Level = term,
        OR = exp(estimate),
        Lower_95CI = exp(estimate - 1.96 * std.error),
        Upper_95CI = exp(estimate + 1.96 * std.error),
        P_Value = p.value
      )
    
    results <- bind_rows(results, model_summary)
  }
  
  return(results)
}

# --- 定义变量和协变量 ---

score_vars_q <- c("score4_an_elnq", "score4_tn_elnq", "score4_wn_elnq","score3_pn_elnq")
score_vars_cont<-c("eln_tn_score_strategy4",
                   "eln_wn_score_strategy4","eln_an_score_strategy4",
                   "eln_pn_score_strategy3")
nhs_db_cohort<- nhs_db_cohort %>%
  mutate(across(all_of(score_vars_cont), ~ as.numeric(scale(.))))
covariates_basic <- c("pa_metbaseq", "alcoq", "aheinnutq", "smoke", "dbfh","factor(hbcbase)", "factor(htnbase)")
covariates_full <- c(covariates_basic, "factor(bmibase3cat)")
results <- list(
  "q_basic" = run_clogit(nhs_db_cohort, score_vars_q),
  "q_adj1" = run_clogit(nhs_db_cohort, score_vars_q, covariates_basic),
  "q_adj2" = run_clogit(nhs_db_cohort, score_vars_q, covariates_full),
  "cont_basic" = run_clogit(nhs_db_cohort, score_vars_cont),
  "cont_adj1" = run_clogit(nhs_db_cohort, score_vars_cont, covariates_basic),
  "cont_adj2" = run_clogit(nhs_db_cohort, score_vars_cont, covariates_full)
)
print(results)
df <- as.data.frame(do.call(rbind, results))
setwd("/udd/n2gji/micro/nut_final/nut_Dec/output")
write.csv(df, "nhs_db2_results.csv")
score_variables <- c(
  "eln_tn_score_strategy4",
  "eln_wn_score_strategy4",
  "eln_an_score_strategy4",
  "eln_pn_score_strategy3"
)

# 将 nhs_db_cohort 从宽格式转换为长格式
nhs_db_long <- nhs_db_cohort %>%
  pivot_longer(
    cols = all_of(score_variables), # 指定要转换的列
    names_to = "score_type",        # 新列，用于存放原来的列名（即分数类型）
    values_to = "score_value"       # 新列，用于存放原来的值（即分数）
  )
nhs_db_long$db<-as.factor(nhs_db_long$db)
# 查看转换后的长数据
head(nhs_db_long)
p<-ggplot(nhs_db_long, aes(x = score_value, fill = factor(db))) +
  # 使用 geom_histogram，设置透明度 alpha 使其可见
  # position="identity" 确保直方图在原位置绘制，而不是堆叠
  geom_histogram(alpha = 0.6, position = "identity", bins = 30) + 
  
  # 按 score_type 创建分面，scales="free" 让每个子图有独立的x,y轴
  facet_wrap(~ score_type, scales = "free") +
  
  # 自定义颜色，使其更美观
  scale_fill_manual(values = c("0" = "#0072B2", "1" = "#D55E00")) +
  
  # 添加标题和标签
  labs(
    title = "Distribution of Scores for DB=0 and DB=1",
    subtitle = "Overlapping histograms with transparency",
    x = "Score Value",
    y = "Frequency (Count)",
    fill = "DB Group" # 图例标题
  ) +
  theme_minimal() +
  theme(legend.position = "top")
p
##########Mediation#############
load("/udd/n2gji/micro/nut_final/t2d/dbcc_data.RData")
library(mediation)
merged_db$db<-ifelse(merged_db$db==1,0,1)
covariates_basic <- c("pa_metbaseq", "alcoq", "aheinnutq", "smoke", "dbfh", 
                      "factor(hbcbase)", "factor(htnbase)")

# 2. 构建模型公式
# 中介模型公式 (M ~ X + C)
# 解释：bmicontbase 如何受 eln_tn_score_strategy1 和其他协变量的影响
mediator_formula_str <- paste("bmicontbase ~ eln_tn_score_strategy4 +", 
                              paste(covariates_basic, collapse = " + "))
mediator_formula <- as.formula(mediator_formula_str)

# 结局模型公式 (Y ~ X + M + C + MatchID)
# 解释：db 如何受 eln_tn_score_strategy1, bmicontbase, 其他协变量和匹配ID的影响
# 注意：我们将中介变量bmicontbase和factor(matchid)也加入模型
outcome_formula_str <- paste("db ~ eln_tn_score_strategy4 + bmicontbase +", 
                             paste(covariates_basic, collapse = " + "), 
                             "+ factor(matchid)") # 调整匹配ID
outcome_formula <- as.formula(outcome_formula_str)

# 我们可以打印出来检查一下公式是否正确
print("中介模型公式:")
print(mediator_formula)
print("结局模型公式:")
print(outcome_formula)
model.m <- lm(mediator_formula, data = merged_db)
summary(model.m) # 查看模型结果

# 2. 拟合结局模型 (逻辑回归)
# Y = db, X = eln_tn_score_strategy1, M = bmicontbase
# family = "binomial" 表示这是逻辑回归
# factor(matchid) 用于调整病例对照的匹配设计
model.y <- glm(outcome_formula, data = merged_db, family = "binomial")
summary(model.y) # 查看模型结果

# 3. 运行中介分析
# 使用 mediate() 函数
# treat = 暴露变量的名称
# mediator = 中介变量的名称
# sims = 模拟次数，用于计算置信区间，建议1000次或更多
set.seed(123) # 设置随机种子以保证结果可重复
mediation_results <- mediate(model.m = model.m,
                             model.y = model.y,
                             treat = "eln_tn_score_strategy4",
                             mediator = "bmicontbase",
                             sims = 1000) # 建议1000次模拟

# 4. 查看并解释中介分析结果
summary(mediation_results)



# 2. 构建模型公式
# 中介模型公式 (M ~ X + C)
# 解释：bmicontbase 如何受 eln_tn_score_strategy1 和其他协变量的影响
mediator_formula_str <- paste("bmicontbase ~ eln_wn_score_strategy4 +", 
                              paste(covariates_basic, collapse = " + "))
mediator_formula <- as.formula(mediator_formula_str)

# 结局模型公式 (Y ~ X + M + C + MatchID)
# 解释：db 如何受 eln_tn_score_strategy1, bmicontbase, 其他协变量和匹配ID的影响
# 注意：我们将中介变量bmicontbase和factor(matchid)也加入模型
outcome_formula_str <- paste("db ~ eln_wn_score_strategy4 + bmicontbase +", 
                             paste(covariates_basic, collapse = " + "), 
                             "+ factor(matchid)") # 调整匹配ID
outcome_formula <- as.formula(outcome_formula_str)

# 我们可以打印出来检查一下公式是否正确
print("中介模型公式:")
print(mediator_formula)
print("结局模型公式:")
print(outcome_formula)
model.m <- lm(mediator_formula, data = merged_db)
summary(model.m) # 查看模型结果

# 2. 拟合结局模型 (逻辑回归)
# Y = db, X = eln_tn_score_strategy1, M = bmicontbase
# family = "binomial" 表示这是逻辑回归
# factor(matchid) 用于调整病例对照的匹配设计
model.y <- glm(outcome_formula, data = merged_db, family = "binomial")
summary(model.y) # 查看模型结果

# 3. 运行中介分析
# 使用 mediate() 函数
# treat = 暴露变量的名称
# mediator = 中介变量的名称
# sims = 模拟次数，用于计算置信区间，建议1000次或更多
set.seed(123) # 设置随机种子以保证结果可重复
mediation_results_wn <- mediate(model.m = model.m,
                                model.y = model.y,
                                treat = "eln_wn_score_strategy4",
                                mediator = "bmicontbase",
                                sims = 1000) # 建议1000次模拟

# 4. 查看并解释中介分析结果
summary(mediation_results_wn)






mediator_formula_str <- paste("bmicontbase ~ eln_an_score_strategy4 +", 
                              paste(covariates_basic, collapse = " + "))
mediator_formula <- as.formula(mediator_formula_str)

# 结局模型公式 (Y ~ X + M + C + MatchID)
# 解释：db 如何受 eln_tn_score_strategy1, bmicontbase, 其他协变量和匹配ID的影响
# 注意：我们将中介变量bmicontbase和factor(matchid)也加入模型
outcome_formula_str <- paste("db ~ eln_an_score_strategy4 + bmicontbase +", 
                             paste(covariates_basic, collapse = " + "), 
                             "+ factor(matchid)") # 调整匹配ID
outcome_formula <- as.formula(outcome_formula_str)

# 我们可以打印出来检查一下公式是否正确
print("中介模型公式:")
print(mediator_formula)
print("结局模型公式:")
print(outcome_formula)
model.m <- lm(mediator_formula, data = merged_db)
summary(model.m) # 查看模型结果

# 2. 拟合结局模型 (逻辑回归)
# Y = db, X = eln_tn_score_strategy1, M = bmicontbase
# family = "binomial" 表示这是逻辑回归
# factor(matchid) 用于调整病例对照的匹配设计
model.y <- glm(outcome_formula, data = merged_db, family = "binomial")
summary(model.y) # 查看模型结果

# 3. 运行中介分析
# 使用 mediate() 函数
# treat = 暴露变量的名称
# mediator = 中介变量的名称
# sims = 模拟次数，用于计算置信区间，建议1000次或更多
set.seed(123) # 设置随机种子以保证结果可重复
mediation_results_an <- mediate(model.m = model.m,
                                model.y = model.y,
                                treat = "eln_an_score_strategy4",
                                mediator = "bmicontbase",
                                sims = 1000) # 建议1000次模拟

# 4. 查看并解释中介分析结果
summary(mediation_results_an)






mediator_formula_str <- paste("bmicontbase ~ eln_pn_score_strategy3 +", 
                              paste(covariates_basic, collapse = " + "))
mediator_formula <- as.formula(mediator_formula_str)

# 结局模型公式 (Y ~ X + M + C + MatchID)
# 解释：db 如何受 eln_tn_score_strategy1, bmicontbase, 其他协变量和匹配ID的影响
# 注意：我们将中介变量bmicontbase和factor(matchid)也加入模型
outcome_formula_str <- paste("db ~ eln_pn_score_strategy3 + bmicontbase +", 
                             paste(covariates_basic, collapse = " + "), 
                             "+ factor(matchid)") # 调整匹配ID
outcome_formula <- as.formula(outcome_formula_str)

# 我们可以打印出来检查一下公式是否正确
print("中介模型公式:")
print(mediator_formula)
print("结局模型公式:")
print(outcome_formula)
model.m <- lm(mediator_formula, data = merged_db)
summary(model.m) # 查看模型结果

# 2. 拟合结局模型 (逻辑回归)
# Y = db, X = eln_tn_score_strategy1, M = bmicontbase
# family = "binomial" 表示这是逻辑回归
# factor(matchid) 用于调整病例对照的匹配设计
model.y <- glm(outcome_formula, data = merged_db, family = "binomial")
summary(model.y) # 查看模型结果

# 3. 运行中介分析
# 使用 mediate() 函数
# treat = 暴露变量的名称
# mediator = 中介变量的名称
# sims = 模拟次数，用于计算置信区间，建议1000次或更多
set.seed(123) # 设置随机种子以保证结果可重复
mediation_results_pn <- mediate(model.m = model.m,
                                model.y = model.y,
                                treat = "eln_pn_score_strategy3",
                                mediator = "bmicontbase",
                                sims = 1000) # 建议1000次模拟

# 4. 查看并解释中介分析结果
summary(mediation_results_pn)

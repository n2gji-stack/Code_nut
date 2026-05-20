# ---------------------------------------------------------------------------------------------------------

# Purpose:           1. Get data get from the sas (diet,biomarkers,basic) and several R (different micorbiome set data) together;
#					           2. Processing and Unifying all the varaible and make it into several longdata seperatedly for 
#                       Taxon, DNA enzyme, RNA enzyme, DNA pathways, RNA pathways, each with annotations too
#                       Format will be row-sample*time, column-metadata+microbiome 

#                    The file will be ready used for IPA analysis in the MLVS

#Study period:       2012-2013                                                       
#Var of interest:    merge and organize of all datasets
#Programmer:         Guliyeerke Jigeer (Based on Yang Hu)                                                                                                          

# ---------------------------------------------------------------------------------------------------------;
rm(list=ls())
setwd("/udd/n2gji/micro/data")
source("MLVSfunctions.R")
##Read in MLVS metabolites data
library(chanmetab)

lvs.lnz.new <- merge_metab_data(endpoints = c("lvs"),
                                collection_to_use = "substudy",
                                cohorts="hpfs",
                                # methods = c("C8-pos", "HILIC-neg", "HILIC-pos"),
                                transformation = "transform_ln_z_score",
                                keep_unknown_metabolites = F,
                                keep_failed_pm_metabolites = F,
                                impute_cutoff=0.05)
f_lvs=chanmetab::fData(lvs.lnz.new$expr_set$hpfs)
lvs.lnz.new$expr_set$hpfs

temp_met<-as.data.frame(t(lvs.lnz.new$expr_set$hpfs@assayData$exprs))
temp_met$id<-substr(rownames(temp_met),1,6)
temp_phe<-as.data.frame(lvs.lnz.new$expr_set$hpfs@phenoData@data)
temp_phe$id<-substr(temp_phe$id,1,6)
lvs_hpfs_pm<-merge(temp_phe,temp_met)
lvs_hpfs_fea<-as.data.frame(lvs.lnz.new$expr_set$hpfs@featureData@data)
# ---------------------------------------------------------------------------;
#
#                        Processing the Phenotype data
#
# ---------------------------------------------------------------------------;
#                      Varaibles for sample collections
# ---------------------------------------------------------------------------;

basic = read.csv("Basic.csv",header=T)
# keep those with at least one fecal date
basic$fecal_sample_date_1 = as.character(basic$fecal_sample_date_1)
basic$fecal_sample_date_2 = as.character(basic$fecal_sample_date_2)
basic$fecal_sample_date_3 = as.character(basic$fecal_sample_date_3)
basic$fecal_sample_date_4 = as.character(basic$fecal_sample_date_4)

basic = basic[which(basic$fecal_sample_date_1!="" | basic$fecal_sample_date_2!="" | 
                      basic$fecal_sample_date_3!="" | basic$fecal_sample_date_4!=""),]

basic[which(basic$fecal_sample_date_1==""),"fecal_sample_date_1"]<-NA
basic[which(basic$fecal_sample_date_2==""),"fecal_sample_date_2"]<-NA
basic[which(basic$fecal_sample_date_3==""),"fecal_sample_date_3"]<-NA
basic[which(basic$fecal_sample_date_4==""),"fecal_sample_date_4"]<-NA

table(rowSums(is.na(basic[,c("fecal_sample_date_1","fecal_sample_date_2","fecal_sample_date_3","fecal_sample_date_4")])))
# 0   1  2   3 
# 212 8 191  4

#head(basic)
dim(basic) # 415 306
Biomarkers = read.csv("Biomarker.csv",header=T)
names(Biomarkers)
dim(Biomarkers) # 686 23

# merge basic and biomarkers
basic = merge(basic,Biomarkers,by="id",all.x=T)
#head(basic)
dim(basic) # 415 328
#------------------------------------------------
# the dateof birth as month to easy calculate age
# also get the day for matching other measures

# check on missing
names(basic)

# calcualte dob
#basic$dob = dateproc("dateofbirth",basic,"months")
#basic$dobday = dateproc("dateofbirth",basic,"days")
#summary(basic$dob) # no missing
#summary(basic$dobday) # no missing
#------------------------------------------------
# date and time of blood draw

summary(basic$irt_blood_mon_bqu1) # missing 8
summary(basic$irt_blood_day_bqu1) # missing 10

summary(basic$irt_blood_mon_bqu2) # missing 11
summary(basic$irt_blood_day_bqu2) # missing 11

# the dates directly provided have wrong ranges - calculate the date myself
# age at bld draw
basic$mon_bqu1 = as.character(basic$mon_bqu1)
basic$mon_bqu2 = as.character(basic$mon_bqu2)

basic$datem_bld1 = date2proc("year_bqu1","mon_bqu1","day_bqu1",basic,"months")
basic$dated_bld1 = date2proc("year_bqu1","mon_bqu1","day_bqu1",basic,"days")

basic$datem_bld2 = date2proc("year_bqu2","mon_bqu2","day_bqu2",basic,"months")
basic$dated_bld2 = date2proc("year_bqu2","mon_bqu2","day_bqu2",basic,"days")

summary(basic$datem_bld1) # missing 8
summary(basic$dated_bld1) # missing 10

summary(basic$datem_bld2) # missing 11
summary(basic$dated_bld2) # missing 11
# the same missing but the ranges are right now
#------------------------------------------------
# date and time of fecal samples

basic$fecal_sample_date_1 = as.character(basic$fecal_sample_date_1)
basic$fecal_sample_date_2 = as.character(basic$fecal_sample_date_2)
basic$fecal_sample_date_3 = as.character(basic$fecal_sample_date_3)
basic$fecal_sample_date_4 = as.character(basic$fecal_sample_date_4)

basic$datem_fec1 = dateproc("fecal_sample_date_1",basic,"months")
basic$dated_fec1 = dateproc("fecal_sample_date_1",basic,"days")

basic$datem_fec2 = dateproc("fecal_sample_date_2",basic,"months")
basic$dated_fec2 = dateproc("fecal_sample_date_2",basic,"days")

basic$datem_fec3 = dateproc("fecal_sample_date_3",basic,"months")
basic$dated_fec3 = dateproc("fecal_sample_date_3",basic,"days")

basic$datem_fec4 = dateproc("fecal_sample_date_4",basic,"months")
basic$dated_fec4 = dateproc("fecal_sample_date_4",basic,"days")

summary(basic$datem_fec1) # missing 120
summary(basic$dated_fec1) # missing 120

summary(basic$datem_fec2) # missing 130
summary(basic$dated_fec2) # missing 130

summary(basic$datem_fec3) # missing 75
summary(basic$dated_fec3) # missing 75

summary(basic$datem_fec4) # missing 77
summary(basic$dated_fec4) # missing 77

#------------------------------------------------
# age at bld and fecal collection

#basic$age_bld1 = round((basic$datem_bld1-basic$dob)/12,1)
#basic$age_bld2 = round((basic$datem_bld2-basic$dob)/12,1)

#basic$age_fec1 = round((basic$datem_fec1-basic$dob)/12,1)
#basic$age_fec2 = round((basic$datem_fec2-basic$dob)/12,1)
#basic$age_fec3 = round((basic$datem_fec3-basic$dob)/12,1)
#basic$age_fec4 = round((basic$datem_fec4-basic$dob)/12,1)

#rbind(
# summary(basic$age_bld1),
#summary(basic$age_bld2))
#      Min. 1st Qu. Median     Mean 3rd Qu. Max. NA's
# [1,] 46.7    66.2  68.90 68.11990  72.600 82.6    8
# [2,] 47.2    66.8  69.65 68.75272  73.125 83.2   11

#rbind(
#  summary(basic$age_fec1),
# summary(basic$age_fec2),
#summary(basic$age_fec3),
#summary(basic$age_fec4))
#      Min. 1st Qu. Median     Mean 3rd Qu. Max. NA's
# [1,] 46.4    66.0  69.00 67.89051   72.55 82.5  120
# [2,] 46.4    66.2  69.00 68.05298   72.60 82.5  130
# [3,] 47.3    66.7  69.50 68.87235   72.90 82.9   75
# [4,] 47.3    66.7  69.45 68.85000   72.80 82.9   77

#------------------------------------------------
# time of fecal sample collection
# 1=morning, 2=afternoon, and 3=night and midnight

basic$time_fec1 = timeproc("stool_time_1",basic)
basic$time_fec2 = timeproc("stool_time_2",basic)
basic$time_fec3 = timeproc("stool_time_3",basic)
basic$time_fec4 = timeproc("stool_time_4",basic)

rbind(
  table(basic$time_fec1),
  table(basic$time_fec2),
  table(basic$time_fec3),
  table(basic$time_fec4))

#       1  2  3
#[1,] 179 39 69
#[2,] 185 24 69
#[3,] 219 32 74
#[4,] 222 34 77

#------------------------------------------------
# fasting status of blood
# 1= <2; 2= 2-4; 3= 5-7; 4= 8-11; 5= 12+

table(basic$lasteat_bqu1) # consider 6 as missing
table(basic$lasteat_bqu2)

basic$fast_bld1 = ifelse(basic$lasteat_bqu1>=4,1,0)
basic$fast_bld2 = ifelse(basic$lasteat_bqu2>=4,1,0)
table(data.frame(fast=basic$fast_bld1,eattime=basic$lasteat_bqu1),useNA="ifany")
#       eattime
# fast     2   4   5   6 <NA>
#   0      1   0   0   0    0
#   1      0 105 287  20    0
#   <NA>   0   0   0   0    2

table(data.frame(fast=basic$fast_bld2,eattime=basic$lasteat_bqu2),useNA="ifany")
#       eattime
# fast     1   3   4   5   6 <NA>
#   0      2   4   0   0   0    0
#   1      0   0 111 269  18    0
#   <NA>   0   0   0   0   0   11

# if missing, consider fasting -> will not lose power on covs
basic$fast_bld1[is.na(basic$fast_bld1)] = 1
basic$fast_bld2[is.na(basic$fast_bld2)] = 1


# ---------------------------------------------------------------------------;
#                   Calculate Basic variables for participants
# ---------------------------------------------------------------------------;

#------------------------------------------------
# smoking status

basic$smoke_bld1 = 0
basic$smoke_bld1[which(basic$smoke_bqu1=="Yes")]=1
basic$smoke_bld2 = 0
basic$smoke_bld2[which(basic$smoke_bqu2=="Yes")]=1

table(basic$smoke_bld1) 
# 7 adn 4 people smoke-> too little

#------------------------------------------------
# BMI

# at PA1
summary(basic$bmi_paq1) 
#   Min. 1st Qu.  Median    Mean 3rd Qu.    Max.    NA's 
#  12.55   23.57   25.52   26.11   28.12   41.77       5

# at PA2
summary(basic$wt_paq2) # missing 239
summary(basic$height_paq1) # missing 186
basic$bmi_paq2 = (basic$wt_paq2*0.453592) / (basic$height_paq1*0.0254) / (basic$height_paq1*0.0254)
summary(basic$bmi_paq2) # missing 13
#   Min. 1st Qu.  Median    Mean 3rd Qu.    Max.    NA's 
#  11.99   23.57   25.53   26.07   28.02   41.05      13


# first blood draw
summary(basic$weight_bqu1) # missing 231
summary(basic$weight_bqu2) # missing 255
basic$bmi_bld1 = (basic$weight_bqu1*0.453592) / (basic$height_paq1*0.0254) / (basic$height_paq1*0.0254)
basic$bmi_bld2 = (basic$weight_bqu2*0.453592) / (basic$height_paq1*0.0254) / (basic$height_paq1*0.0254)
summary(basic$bmi_bld1) # missing 239
summary(basic$bmi_bld2) # missing 262

#   Min. 1st Qu.  Median    Mean 3rd Qu.    Max.    NA's 
#  18.37   23.46   25.54   26.11   28.13   40.42      16

#   Min. 1st Qu.  Median    Mean 3rd Qu.    Max.    NA's 
#  18.47   23.60   25.54   26.17   28.06   41.50      18

# if missing can be impute by pad pa
basic[is.na(basic$bmi_bld1),"bmi_bld1"] = basic[is.na(basic$bmi_bld1),"bmi_paq1"]
basic[is.na(basic$bmi_bld2),"bmi_bld2"] = basic[is.na(basic$bmi_bld2),"bmi_paq2"]

# if still missing, use the other time to impute -> will not lose power on covs
basic[is.na(basic$bmi_bld1),"bmi_bld1"] = basic[is.na(basic$bmi_bld1),"bmi_bld2"]
basic[is.na(basic$bmi_bld2),"bmi_bld2"] = basic[is.na(basic$bmi_bld2),"bmi_bld1"]
basic[is.na(basic$bmi_bld1),"bmi_bld1"] = basic[is.na(basic$bmi_bld1),"bmi_paq2"]
basic[is.na(basic$bmi_bld2),"bmi_bld2"] = basic[is.na(basic$bmi_bld2),"bmi_paq1"]

# if still missing, use median -> will not lose power on covs, 4 still missing, it's fine
basic[is.na(basic$bmi_bld1),"bmi_bld1"] = median(basic[,"bmi_bld1"],na.rm=TRUE)
basic[is.na(basic$bmi_bld2),"bmi_bld2"] = median(basic[,"bmi_bld2"],na.rm=TRUE)

rbind(
  summary(basic$bmi_bld1),
  summary(basic$bmi_bld2)) 

#          Min.  1st Qu.   Median     Mean  3rd Qu.     Max.
# [1,] 18.36588 23.42317 25.52304 26.07795 28.05406 40.41562
# [2,] 18.47057 23.58442 25.54006 26.15082 28.09331 41.50060


#------------------------------------------------
# physical activity

rbind(
  summary(basic$totMETs_paq1),
  summary(basic$totMETs_paq2))

#          Min. 1st Qu.   Median     Mean  3rd Qu.     Max. NA's
# [1,] 11.30206 79.6240 112.9809 117.8041 151.6844 352.9783    2
# [2,] 14.87406 81.5667 115.9390 122.4292 152.6404 475.0892    8

# if missimg, use value another time and then median to impute 
basic$totMETs_paq1[is.na(basic$totMETs_paq1)] = basic$totMETs_paq2[is.na(basic$totMETs_paq1)]
basic$totMETs_paq2[is.na(basic$totMETs_paq2)] = basic$totMETs_paq1[is.na(basic$totMETs_paq2)]

rbind(
  summary(basic$totMETs_paq1),
  summary(basic$totMETs_paq2))

#         Min.  1st Qu.   Median     Mean  3rd Qu.     Max.
#[1,] 11.30206 79.29955 112.9809 117.7192 151.1911 352.9783
#[2,] 14.87406 81.56670 116.0424 122.7477 152.7237 475.0892


#------------------------------------------------
# medications at blood collections

basic$steroid48h_bld1 = ifelse(basic$ster48hr_bqu1=="",NA,ifelse(basic$ster48hr_bqu1=="no",0,1))
basic$steroidmth_bld1 = ifelse(basic$stermo_bqu1=="",NA,ifelse(basic$stermo_bqu1=="no",0,1))

basic$antibio48h_bld1 = ifelse(basic$antibi48hr_bqu1=="",NA,ifelse(basic$antibi48hr_bqu1=="no",0,1))
basic$antibiomth_bld1 = ifelse(basic$antibimo_bqu1=="",NA,ifelse(basic$antibimo_bqu1=="no",0,1))

basic$aspirin48h_bld1 = ifelse(basic$aspirin48hr_bqu1=="",NA,ifelse(basic$aspirin48hr_bqu1=="no",0,1))
basic$aspirinmth_bld1 = ifelse(basic$aspirinmo_bqu1=="",NA,ifelse(basic$aspirinmo_bqu1=="no",0,1))

basic$acetamp48h_bld1 = ifelse(basic$aceta48hr_bqu1=="",NA,ifelse(basic$aceta48hr_bqu1=="no",0,1))
basic$acetampmth_bld1 = ifelse(basic$acetamo_bqu1=="",NA,ifelse(basic$acetamo_bqu1=="no",0,1))

basic$ibuprof48h_bld1 = ifelse(basic$ibupro48hr_bqu1=="",NA,ifelse(basic$ibupro48hr_bqu1=="no",0,1))
basic$ibuprofmth_bld1 = ifelse(basic$ibupromo_bqu1=="",NA,ifelse(basic$ibupromo_bqu1=="no",0,1))

basic$h2block48h_bld1 = ifelse(basic$h2bloc48hr_bqu1=="",NA,ifelse(basic$h2bloc48hr_bqu1=="no",0,1))
basic$h2blockmth_bld1 = ifelse(basic$h2blocmo_bqu1=="",NA,ifelse(basic$h2blocmo_bqu1=="no",0,1))

rbind(
  table(basic$steroid48h_bld1,useNA="ifany"),
  table(basic$steroidmth_bld1,useNA="ifany"),
  table(basic$antibio48h_bld1,useNA="ifany"),
  table(basic$antibiomth_bld1,useNA="ifany"),
  table(basic$aspirin48h_bld1,useNA="ifany"),
  table(basic$aspirinmth_bld1,useNA="ifany"),
  table(basic$acetamp48h_bld1,useNA="ifany"),
  table(basic$acetampmth_bld1,useNA="ifany"),
  table(basic$ibuprof48h_bld1,useNA="ifany"),
  table(basic$ibuprofmth_bld1,useNA="ifany"),
  table(basic$h2block48h_bld1,useNA="ifany"),
  table(basic$h2blockmth_bld1,useNA="ifany"))

#         0   1
# [1,] 413   1 # ster
# [2,] 410   4 # ster
# [3,] 406   8 # antibiotics
# [4,] 401  13 # antibiotics
# [5,] 238 176 # asp
# [6,] 341  73 # asp
# [7,] 392  22 # acet
# [8,] 376  38 # acet
# [9,] 371  43 # ibupro
# [10,] 322  92 # ibupro
# [11,] 401  13 # h2blocker
# [12,] 407   7 # h2blocker

#         0   1
#  [1,] 689   1  
#  [2,] 684   6  
#  [3,] 673  17  
#  [4,] 672  18  
#  [5,] 414 276  
#  [6,] 570 120  
#  [7,] 655  35  
#  [8,] 625  65  
#  [9,] 606  84  
# [10,] 537 153  
# [11,] 667  23  
# [12,] 676  14  

basic$steroid48h_bld2 = ifelse(basic$ster48hr_bqu2=="",NA,ifelse(basic$ster48hr_bqu2=="no",0,1))
basic$steroidmth_bld2 = ifelse(basic$stermo_bqu2=="",NA,ifelse(basic$stermo_bqu2=="no",0,1))

basic$antibio48h_bld2 = ifelse(basic$antibi48hr_bqu2=="",NA,ifelse(basic$antibi48hr_bqu2=="no",0,1))
basic$antibiomth_bld2 = ifelse(basic$antibimo_bqu2=="",NA,ifelse(basic$antibimo_bqu2=="no",0,1))

basic$aspirin48h_bld2 = ifelse(basic$aspirin48hr_bqu2=="",NA,ifelse(basic$aspirin48hr_bqu2=="no",0,1))
basic$aspirinmth_bld2 = ifelse(basic$aspirinmo_bqu2=="",NA,ifelse(basic$aspirinmo_bqu2=="no",0,1))

basic$acetamp48h_bld2 = ifelse(basic$aceta48hr_bqu2=="",NA,ifelse(basic$aceta48hr_bqu2=="no",0,1))
basic$acetampmth_bld2 = ifelse(basic$acetamo_bqu2=="",NA,ifelse(basic$acetamo_bqu2=="no",0,1))

basic$ibuprof48h_bld2 = ifelse(basic$ibupro48hr_bqu2=="",NA,ifelse(basic$ibupro48hr_bqu2=="no",0,1))
basic$ibuprofmth_bld2 = ifelse(basic$ibupromo_bqu2=="",NA,ifelse(basic$ibupromo_bqu2=="no",0,1))

basic$h2block48h_bld2 = ifelse(basic$h2bloc48hr_bqu2=="",NA,ifelse(basic$h2bloc48hr_bqu2=="no",0,1))
basic$h2blockmth_bld2 = ifelse(basic$h2blocmo_bqu2=="",NA,ifelse(basic$h2blocmo_bqu2=="no",0,1))

rbind(
  table(basic$steroid48h_bld2),
  table(basic$steroidmth_bld2),
  table(basic$antibio48h_bld2),
  table(basic$antibiomth_bld2),
  table(basic$aspirin48h_bld2),
  table(basic$aspirinmth_bld2),
  table(basic$acetamp48h_bld2),
  table(basic$acetampmth_bld2),
  table(basic$ibuprof48h_bld2),
  table(basic$ibuprofmth_bld2),
  table(basic$h2block48h_bld2),
  table(basic$h2blockmth_bld2))

#         0   1
# [1,] 402   2 # ster
# [2,] 402   2 # ster
# [3,] 392  12 # antibiotics
# [4,] 385  19 # antibiotics
# [5,] 234 170 # asp
# [6,] 334  70 # asp
# [7,] 378  26 # acet
# [8,] 366  38 # acet
# [9,] 362  42 # ibupro
# [10,] 326  78 # ibupro
# [11,] 386  18 # h2blocker
# [12,] 397   7 # h2blocker


#------------------------------------------------
# medications and others at stool 1

basic$antibio_12m_fec12 = ifelse(basic$ant_12mo_qu1=="",NA,ifelse(basic$ant_12mo_qu1=="Checked",1,0))
basic$chm_12m_fec12 = ifelse(basic$chm_12mo_qu1=="",NA,ifelse(basic$chm_12mo_qu1=="Checked",1,0))
basic$ims_12m_fec12 = ifelse(basic$ims_12mon_qu1=="",NA,ifelse(basic$ims_12mon_qu1=="Checked",1,0))
basic$colsc_2m_fec12 = ifelse(basic$colsc_2mo_qu1=="",NA,ifelse(basic$colsc_2mo_qu1=="Yes",1,0))
basic$ctscan_2m_fec12 = ifelse(basic$ctscan_2mo_qu1=="",NA,ifelse(basic$ctscan_2mo_qu1=="Yes",1,0))
basic$cdiarr_2m_fec12 = ifelse(basic$cdiarr_2mo_qu1=="",NA,ifelse(basic$cdiarr_2mo_qu1=="Yes",1,0))
basic$adiarr_2m_fec12 = ifelse(basic$adiarr_2mo_qu1=="",NA,ifelse(basic$adiarr_2mo_qu1=="Yes",1,0))
basic$probio_2m_fec12 = ifelse(basic$probio_2mo_qu1=="",NA,ifelse(basic$probio_2mo_qu1=="Yes",1,0))
basic$acid_2m_fec12 = ifelse(basic$acid_med_2mo_qu1=="",NA,ifelse(basic$acid_med_2mo_qu1=="Yes",1,0))
basic$bile_2m_fec12 = ifelse(basic$bile_med_2mo_qu1=="",NA,ifelse(basic$bile_med_2mo_qu1=="Yes",1,0))


basic$yog_2m_fec12 = NA
basic$yog_2m_fec12[which(basic$yog_2mo_qu1=="1-6 times a week")]=3.5  # 3.5/w
basic$yog_2m_fec12[which(basic$yog_2mo_qu1=="Daily")]=7  # 7/w
basic$yog_2m_fec12[which(basic$yog_2mo_qu1=="More than daily")]=14  # 7/w assuming 2 per day
basic$yog_2m_fec12[which(basic$yog_2mo_qu1=="Never")]=0  # 0/w
basic$yog_2m_fec12[which(basic$yog_2mo_qu1=="Rarely")]=1  # 1/w

table(basic$dietpref_qu1)
basic$meatpref_fec12 = NA
basic$meatpref_fec12[which(basic$dietpref_qu1=="Standard diet")]="standard"
basic$meatpref_fec12[which(basic$dietpref_qu1=="Standard diet with poultry and/or fish (no red meat)")]="standard no red meat" 
basic$meatpref_fec12[which(basic$dietpref_qu1=="Vegan (no meat, dairy, or animal products)")]="vegan"
basic$meatpref_fec12[which(basic$dietpref_qu1=="Vegetarian (no meat)")]="vegetarian" 
table(basic$dietpref_qu1)

basic$stooltype_fec1 = basic$stool_type_batch1_qu1
basic$stooltype_fec2 = basic$stool_type_batch2_qu1

rbind(
  table(basic$stooltype_fec1,useNA="ifany"),
  table(basic$stooltype_fec2,useNA="ifany"))
#       1  2  3   4  5  6 7 <NA>
# fec1 10 27 74 131 24 21 5  123
# fec2 12 31 81 120 15 14 6  136

# dummy var on type 1-6

basic$stooltype_fec1.1 = ifelse(is.na(basic$stooltype_fec1),0,ifelse(basic$stooltype_fec1==1,1,0))
basic$stooltype_fec1.2 = ifelse(is.na(basic$stooltype_fec1),0,ifelse(basic$stooltype_fec1==2,1,0))
basic$stooltype_fec1.3 = ifelse(is.na(basic$stooltype_fec1),0,ifelse(basic$stooltype_fec1==3,1,0))
basic$stooltype_fec1.4 = ifelse(is.na(basic$stooltype_fec1),0,ifelse(basic$stooltype_fec1==4,1,0))
basic$stooltype_fec1.5 = ifelse(is.na(basic$stooltype_fec1),0,ifelse(basic$stooltype_fec1==5,1,0))
basic$stooltype_fec1.6 = ifelse(is.na(basic$stooltype_fec1),0,ifelse(basic$stooltype_fec1==6,1,0))

basic$stooltype_fec2.1 = ifelse(is.na(basic$stooltype_fec2),0,ifelse(basic$stooltype_fec2==1,1,0))
basic$stooltype_fec2.2 = ifelse(is.na(basic$stooltype_fec2),0,ifelse(basic$stooltype_fec2==2,1,0))
basic$stooltype_fec2.3 = ifelse(is.na(basic$stooltype_fec2),0,ifelse(basic$stooltype_fec2==3,1,0))
basic$stooltype_fec2.4 = ifelse(is.na(basic$stooltype_fec2),0,ifelse(basic$stooltype_fec2==4,1,0))
basic$stooltype_fec2.5 = ifelse(is.na(basic$stooltype_fec2),0,ifelse(basic$stooltype_fec2==5,1,0))
basic$stooltype_fec2.6 = ifelse(is.na(basic$stooltype_fec2),0,ifelse(basic$stooltype_fec2==6,1,0))

rbind(
  table(basic$stooltype_fec1.1),
  table(basic$stooltype_fec1.2),
  table(basic$stooltype_fec1.3),
  table(basic$stooltype_fec1.4),
  table(basic$stooltype_fec1.5),
  table(basic$stooltype_fec1.6))

#        0   1
# [1,] 405  10
# [2,] 388  27
# [3,] 341  74
# [4,] 284 131
# [5,] 391  24
# [6,] 394  21

rbind(
  table(basic$stooltype_fec2.1),
  table(basic$stooltype_fec2.2),
  table(basic$stooltype_fec2.3),
  table(basic$stooltype_fec2.4),
  table(basic$stooltype_fec2.5),
  table(basic$stooltype_fec2.6))

#       0   1
# [1,] 403  12
# [2,] 384  31
# [3,] 334  81
# [4,] 295 120
# [5,] 400  15
# [6,] 401  14

#------------------------------------------------
# medications and others at stool 2

basic$antibio_12m_fec34 = ifelse(basic$ant_12mo_qu2=="",NA,ifelse(basic$ant_12mo_qu2=="Checked",1,0))
basic$chm_12m_fec34 = ifelse(basic$chm_12mo_qu2=="",NA,ifelse(basic$chm_12mo_qu2=="Checked",1,0))
basic$ims_12m_fec34 = ifelse(basic$ims_12mon_qu2=="",NA,ifelse(basic$ims_12mon_qu2=="Checked",1,0))
basic$colsc_2m_fec34 = ifelse(basic$colsc_2mo_qu2=="",NA,ifelse(basic$colsc_2mo_qu2=="Yes",1,0))
basic$ctscan_2m_fec34 = ifelse(basic$ctscan_2mo_qu2=="",NA,ifelse(basic$ctscan_2mo_qu2=="Yes",1,0))
basic$cdiarr_2m_fec34 = ifelse(basic$cdiarr_2mo_qu2=="",NA,ifelse(basic$cdiarr_2mo_qu2=="Yes",1,0))
basic$adiarr_2m_fec34 = ifelse(basic$adiarr_2mo_qu2=="",NA,ifelse(basic$adiarr_2mo_qu2=="Yes",1,0))
basic$probio_2m_fec34 = ifelse(basic$probio_2mo_qu2=="",NA,ifelse(basic$probio_2mo_qu2=="Yes",1,0))
basic$acid_2m_fec34 = ifelse(basic$acid_med_2mo_qu2=="",NA,ifelse(basic$acid_med_2mo_qu2=="Yes",1,0))
basic$bile_2m_fec34 = ifelse(basic$bile_med_2mo_qu2=="",NA,ifelse(basic$bile_med_2mo_qu2=="Yes",1,0))


basic$yog_2m_fec34 = NA
basic$yog_2m_fec34[which(basic$yog_2mo_qu2=="1-6 times a week")]=3.5  # 3.5/w
basic$yog_2m_fec34[which(basic$yog_2mo_qu2=="Daily")]=7  # 7/w
basic$yog_2m_fec34[which(basic$yog_2mo_qu2=="More than daily")]=14  # 7/w assuming 2 per day
basic$yog_2m_fec34[which(basic$yog_2mo_qu2=="Never")]=0  # 0/w
basic$yog_2m_fec34[which(basic$yog_2mo_qu2=="Rarely")]=1  # 1/w

table(basic$dietpref_qu2)
basic$meatpref_fec34 = NA
basic$meatpref_fec34[which(basic$dietpref_qu2=="Standard diet")]="standard"
basic$meatpref_fec34[which(basic$dietpref_qu2=="Standard diet with poultry and/or fish (no red meat)")]="standard no red meat" 
basic$meatpref_fec34[which(basic$dietpref_qu2=="Vegan (no meat, dairy, or animal products)")]="vegan"
basic$meatpref_fec34[which(basic$dietpref_qu2=="Vegetarian (no meat)")]="vegetarian" 


basic$stooltype_fec3 = basic$stool_type_batch3_qu2
basic$stooltype_fec4 = basic$stool_type_batch4_qu2

rbind(
  table(basic$stooltype_fec3,useNA="ifany"),
  table(basic$stooltype_fec4,useNA="ifany"))

# dummy var on type 1-6

basic$stooltype_fec3.1 = ifelse(is.na(basic$stooltype_fec3),0,ifelse(basic$stooltype_fec3==1,1,0))
basic$stooltype_fec3.2 = ifelse(is.na(basic$stooltype_fec3),0,ifelse(basic$stooltype_fec3==2,1,0))
basic$stooltype_fec3.3 = ifelse(is.na(basic$stooltype_fec3),0,ifelse(basic$stooltype_fec3==3,1,0))
basic$stooltype_fec3.4 = ifelse(is.na(basic$stooltype_fec3),0,ifelse(basic$stooltype_fec3==4,1,0))
basic$stooltype_fec3.5 = ifelse(is.na(basic$stooltype_fec3),0,ifelse(basic$stooltype_fec3==5,1,0))
basic$stooltype_fec3.6 = ifelse(is.na(basic$stooltype_fec3),0,ifelse(basic$stooltype_fec3==6,1,0))

basic$stooltype_fec4.1 = ifelse(is.na(basic$stooltype_fec4),0,ifelse(basic$stooltype_fec4==1,1,0))
basic$stooltype_fec4.2 = ifelse(is.na(basic$stooltype_fec4),0,ifelse(basic$stooltype_fec4==2,1,0))
basic$stooltype_fec4.3 = ifelse(is.na(basic$stooltype_fec4),0,ifelse(basic$stooltype_fec4==3,1,0))
basic$stooltype_fec4.4 = ifelse(is.na(basic$stooltype_fec4),0,ifelse(basic$stooltype_fec4==4,1,0))
basic$stooltype_fec4.5 = ifelse(is.na(basic$stooltype_fec4),0,ifelse(basic$stooltype_fec4==5,1,0))
basic$stooltype_fec4.6 = ifelse(is.na(basic$stooltype_fec4),0,ifelse(basic$stooltype_fec4==6,1,0))

rbind(
  table(basic$stooltype_fec3.1),
  table(basic$stooltype_fec3.2),
  table(basic$stooltype_fec3.3),
  table(basic$stooltype_fec3.4),
  table(basic$stooltype_fec3.5),
  table(basic$stooltype_fec3.6))

#        0   1
# [1,] 402  13
# [2,] 380  35
# [3,] 332  83
# [4,] 275 140
# [5,] 387  28
# [6,] 391  24

rbind(
  table(basic$stooltype_fec4.1),
  table(basic$stooltype_fec4.2),
  table(basic$stooltype_fec4.3),
  table(basic$stooltype_fec4.4),
  table(basic$stooltype_fec4.5),
  table(basic$stooltype_fec4.6))

#        0   1
# [1,] 400  15
# [2,] 382  33
# [3,] 341  74
# [4,] 270 145
# [5,] 384  31
# [6,] 383  32



#------------------------------------------------
# long-term diatery patterns, whole grain, red meat
# calcualte cumulate average 
basic$lt_fiber = rowMeans(basic[,c("aofib86a", "aofib90a", "aofib94a", "aofib98a", "aofib02a", "aofib06a", "aofib10a")],na.rm=TRUE)
basic$lt_whgrn = rowMeans(basic[,c("whgrn86a", "whgrn90a", "whgrn94a", "whgrn98a", "whgrn02a", "whgrn06a", "whgrn10a")],na.rm=TRUE)
basic$lt_ahei = rowMeans(basic[,c("ahei86_a", "ahei90_a", "ahei94_a", "ahei98_a", "ahei02_a", "ahei06_a", "ahei10_a")],na.rm=TRUE)
basic$lt_ahei_na = rowMeans(basic[,c("ahei86_na", "ahei90_na", "ahei94_na", "ahei98_na", "ahei02_na", "ahei06_na", "ahei10_na")],na.rm=TRUE)
basic$lt_ahei_nowgr = rowMeans(basic[,c("ahei86_nowgr", "ahei90_nowgr", "ahei94_nowgr", "ahei98_nowgr", "ahei02_nowgr", "ahei06_nowgr", "ahei10_nowgr")],na.rm=TRUE)
basic$lt_readmt = rowMeans(basic[,c("rmea86w", "rmea90w", "rmea94w", "rmea98w", "rmea02w", "rmea06w", "rmea10w")],na.rm=TRUE)
basic$lt_nut = rowMeans(basic[,c("nut86", "nut90", "nut94", "nut98", "nut02", "nut06", "nut10")],na.rm=TRUE)
basic$lt_veg = rowMeans(basic[,c("tveg86", "tveg90", "tveg94", "tveg98", "tveg02", "tveg06", "tveg10")],na.rm=TRUE)
basic$lt_fru = rowMeans(basic[,c("fruit86", "fruit90", "fruit94", "fruit98", "fruit02", "fruit06", "fruit10")],na.rm=TRUE)
basic$lt_trypt = rowMeans(basic[,c("trypt86n", "trypt90n", "trypt94n", "trypt98n", "trypt02n", "trypt06n", "trypt10n")],na.rm=TRUE)

rbind(
  summary(basic$lt_fiber),
  summary(basic$lt_whgrn),
  summary(basic$lt_ahei),
  summary(basic$lt_ahei_na),
  summary(basic$lt_ahei_nowgr),
  summary(basic$lt_readmt),
  summary(basic$lt_nut),
  summary(basic$lt_veg),
  summary(basic$lt_fru),
  summary(basic$lt_trypt))
# for individuals missing diet, use median to replace it because it is only for adjustment
basic$lt_fiber[is.na(basic$lt_fiber)] = median(basic$lt_fiber,na.rm=TRUE)
basic$lt_whgrn[is.na(basic$lt_whgrn)] = median(basic$lt_whgrn,na.rm=TRUE)
basic$lt_ahei[is.na(basic$lt_ahei)] = median(basic$lt_ahei,na.rm=TRUE)
basic$lt_ahei_na[is.na(basic$lt_ahei_na)] = median(basic$lt_ahei_na,na.rm=TRUE)
basic$lt_ahei_nowgr[is.na(basic$lt_ahei_nowgr)] = median(basic$lt_ahei_nowgr,na.rm=TRUE)
basic$lt_readmt[is.na(basic$lt_readmt)] = median(basic$lt_readmt,na.rm=TRUE)
basic$lt_nut[is.na(basic$lt_nut)] = median(basic$lt_nut,na.rm=TRUE)
basic$lt_veg[is.na(basic$lt_veg)] = median(basic$lt_veg,na.rm=TRUE)
basic$lt_fru[is.na(basic$lt_fru)] = median(basic$lt_fru,na.rm=TRUE)
basic$lt_trypt[is.na(basic$lt_trypt)] = median(basic$lt_trypt,na.rm=TRUE)

rbind(
  summary(basic$lt_fiber),  
  summary(basic$lt_whgrn),
  summary(basic$lt_ahei),
  summary(basic$lt_ahei_na),
  summary(basic$lt_ahei_nowgr),
  summary(basic$lt_readmt),
  summary(basic$lt_nut),
  summary(basic$lt_veg),
  summary(basic$lt_fru),
  summary(basic$lt_trypt))
# Min.    1st Qu.     Median       Mean    3rd Qu.       Max.
# [1,]  6.2157143 28.0400000 32.3385714 34.5884464 39.5671429 141.194286
# [2,] 26.4614410 53.7588039 57.5406547 58.1610204 62.7056533  84.686264
# [3,] 23.2471553 47.5250727 51.9774852 52.0151667 55.7247181  80.043407
# [4,] 25.1277902 50.3746898 53.6339512 54.2746080 58.8201462  79.828804
# [5,]  0.0000000  0.6950000  0.9009524  0.9729771  1.1528571   4.505714
# [6,]  0.0100000  0.4164286  0.5242857  0.6383462  0.7539286   2.437143
# [7,]  0.8282429  3.0730071  3.7512571  3.7874890  4.2764500  10.673129
# [8,]  0.0000000  1.1521429  1.4442857  1.6102547  1.9121429   6.731429
# [9,]  0.4800000  0.8907143  0.9942857  1.0042679  1.1000000   1.892857

#------------------------------------------------
# main outcomes TMAO Biomarkers - use adjusted and batch corrected biomarker levels

rbind(
  summary(basic$adj_choline1),
  summary(basic$adj_tmao1),
  summary(basic$adj_carnitine1))

#            Min.   1st Qu.    Median      Mean   3rd Qu.     Max. NA's
# [1,]  7.8652108 14.185211 16.001450 16.345984 18.483532 31.08697    6
# [2,]  0.2378718  2.659617  3.766012  5.965381  6.045487 95.12899    6
# [3,] 20.8957661 33.459071 39.124598 39.516096 44.169475 79.53857    6
par(mar=c(1,1,1,1))

lineardis("adj_choline1",basic,3)
lineardis("adj_choline1",basic,4) # use log

lineardis("adj_tmao1",basic,3) 
lineardis("adj_tmao1",basic,4) # use log

lineardis("adj_carnitine1",basic,3) 
lineardis("adj_carnitine1",basic,4) # use log


rbind(
  summary(basic$adj_choline2),
  summary(basic$adj_tmao2),
  summary(basic$adj_carnitine2))

#            Min.   1st Qu.    Median      Mean   3rd Qu.     Max. NA's
# [1,]  8.6486877 14.176395 16.047315 16.490332 18.492523 29.53252   10
# [2,]  0.1845096  2.575375  3.748037  5.895166  6.205797 92.76132   10
# [3,] 22.0843911 34.032773 39.408879 39.891876 44.738170 73.11639   10

lineardis("adj_choline2",basic,3)
lineardis("adj_choline2",basic,4) # use log

lineardis("adj_tmao2",basic,3)
lineardis("adj_tmao2",basic,4) # use log

lineardis("adj_carnitine2",basic,3)
lineardis("adj_carnitine2",basic,4) # use log


# standardizaiton should be done on all values including 1 and 2

basic$choline_std_bld1 = standwinso("adj_choline1","adj_choline2",basic,"log",4)[,1]
basic$tmao_std_bld1 = standwinso("adj_tmao1","adj_tmao2",basic,"log",4)[,1]
basic$carnitine_std_bld1 = standwinso("adj_carnitine1","adj_carnitine2",basic,"log",4)[,1]

basic$choline_std_bld2 = standwinso("adj_choline1","adj_choline2",basic,"log",4)[,2]
basic$tmao_std_bld2 = standwinso("adj_tmao1","adj_tmao2",basic,"log",4)[,2]
basic$carnitine_std_bld2 = standwinso("adj_carnitine1","adj_carnitine2",basic,"log",4)[,2]

# check final distributions

par(mfrow=c(2,3))
hist(basic$choline_std_bld1)
hist(basic$tmao_std_bld1)
hist(basic$carnitine_std_bld1)
hist(basic$choline_std_bld2)
hist(basic$tmao_std_bld2)
hist(basic$carnitine_std_bld2)

# ready to use now


#------------------------------------------------
# blood lipids, CRP, and HbA1c - use orignal except for crp which should use batch corrected values

rbind(
  summary(basic$tc_plasma1),
  summary(basic$hdlc_plasma1),
  summary(basic$tg_plasma1),
  summary(basic$adj_crp1),
  summary(basic$hba1cp1))

#              Min.     1st Qu.      Median       Mean    3rd Qu.      Max. NA's
# [1,] 97.000000000 160.0000000 181.0000000 185.689487 205.000000 406.00000    6
# [2,] 27.600000000  45.5000000  55.2000000  56.622249  65.900000 110.80000    6
# [3,] 33.000000000  67.0000000  89.0000000 101.990220 122.000000 482.00000    6
# [4,]  0.003839386   0.2892978   0.6859276   1.583236   1.524472  21.77817    6
# [5,]  4.820000000   5.4900000   5.6800000   5.735966   5.890000   9.16000    6

lineardis("tc_plasma1",basic,3) # 
lineardis("tc_plasma1",basic,4) # use log

lineardis("hdlc_plasma1",basic,3) #
lineardis("hdlc_plasma1",basic,4) # use log

lineardis("tg_plasma1",basic,3) #
lineardis("tg_plasma1",basic,4) # use log

lineardis("adj_crp1",basic,3) # 
lineardis("adj_crp1",basic,4) # use log, but has two hits distribution due to very low levels of some samples

lineardis("hba1cp1",basic,3) #
lineardis("hba1cp1",basic,4) # use log

# lineardis("adj_el1",basic,3) #
# lineardis("adj_el1",basic,4) # use log

round(rbind(
  summary(basic$tc_plasma2),
  summary(basic$hdlc_plasma2),
  summary(basic$tg_plasma2),
  summary(basic$adj_crp2),
  summary(basic$hba1cp2)))

#      Min. 1st Qu. Median Mean 3rd Qu. Max. NA's
# [1,]  105     160    181  183     203  353   10
# [2,]   27      44     54   55      61   99   10
# [3,]   29      65     88  101     116  468   10
# [4,]    0       0      1    2       2   26   10
# [5,]    5       5      6    6       6    7   10

lineardis("tc_plasma2",basic,3) #
lineardis("tc_plasma2",basic,4) # use log

lineardis("hdlc_plasma2",basic,3) # 
lineardis("hdlc_plasma2",basic,4) # use log

lineardis("tg_plasma2",basic,3) # 
lineardis("tg_plasma2",basic,4) # use log

lineardis("adj_crp2",basic,3) #
lineardis("adj_crp2",basic,4) # use log, but has two hits distribution due to very low levels of some samples

lineardis("hba1cp2",basic,3) #
lineardis("hba1cp2",basic,4) # use log

lineardis("hba1cp2",basic,3) #
lineardis("hba1cp2",basic,4) # use log

# lineardis("adj_el2",basic,3) #
# lineardis("adj_el2",basic,4) # use log

basic$tc.hdl1 = basic$tc_plasma1/basic$hdlc_plasma1
basic$tc.hdl2 = basic$tc_plasma2/basic$hdlc_plasma2
lineardis("tc.hdl1",basic,3) #
lineardis("tc.hdl1",basic,4) # use log
lineardis("tc.hdl2",basic,3) # 
lineardis("tc.hdl2",basic,4) # use log


basic$tc_std_bld1 = standwinso("tc_plasma1","tc_plasma2",basic,"log",4)[,1]
basic$hdlc_std_bld1 = standwinso("hdlc_plasma1","hdlc_plasma2",basic,"log",4)[,1]
basic$tg_std_bld1 = standwinso("tg_plasma1","tg_plasma2",basic,"log",4)[,1]
basic$crp_std_bld1 = standwinso("adj_crp1","adj_crp2",basic,"log",4)[,1]
basic$hba1c_std_bld1 = standwinso("hba1cp1","hba1cp2",basic,"log",4)[,1]
#basic$el_std_bld1 = standwinso("adj_el1","adj_el2",basic,"log",4)[,1]

basic$tc_std_bld2 = standwinso("tc_plasma1","tc_plasma2",basic,"log",4)[,2]
basic$hdlc_std_bld2 = standwinso("hdlc_plasma1","hdlc_plasma2",basic,"log",4)[,2]
basic$tg_std_bld2 = standwinso("tg_plasma1","tg_plasma2",basic,"log",4)[,2]
basic$crp_std_bld2 = standwinso("adj_crp1","adj_crp2",basic,"log",4)[,2]
basic$hba1c_std_bld2 = standwinso("hba1cp1","hba1cp2",basic,"log",4)[,2]
#basic$el_std_bld2 = standwinso("adj_el1","adj_el2",basic,"log",4)[,2]

basic$tc.hdl_std_bld1 = standwinso("tc.hdl1","tc.hdl2",basic,"log",4)[,1]
basic$tc.hdl_std_bld2 = standwinso("tc.hdl1","tc.hdl2",basic,"log",4)[,2]

# check final distributions

par(mfrow=c(4,3))
hist(basic$tc_std_bld1)
hist(basic$hdlc_std_bld1)
hist(basic$tg_std_bld1)
hist(basic$crp_std_bld1)
hist(basic$hba1c_std_bld1)
# hist(basic$el_std_bld1)
hist(basic$tc.hdl_std_bld1)

hist(basic$tc_std_bld2)
hist(basic$hdlc_std_bld2)
hist(basic$tg_std_bld2)
hist(basic$crp_std_bld2)
hist(basic$hba1c_std_bld2)
# hist(basic$el_std_bld2)
hist(basic$tc.hdl_std_bld2)


# ---------------------------------------------------------------------------;
#                        Get basic file to use
# ---------------------------------------------------------------------------;

names(basic)

Basic1 = basic[,c("id","datem_bld1","dated_bld1","datem_fec1","dated_fec1","time_fec1","fast_bld1",
                  "smoke_bld1","bmi_bld1","totMETs_paq1","lt_fiber","lt_whgrn","lt_ahei","lt_ahei_na","lt_ahei_nowgr","lt_readmt", "lt_nut","lt_veg", "lt_fru","lt_trypt",
                  "steroid48h_bld1","steroidmth_bld1","antibio48h_bld1","antibiomth_bld1","aspirin48h_bld1","aspirinmth_bld1",
                  "acetamp48h_bld1","acetampmth_bld1","ibuprof48h_bld1","ibuprofmth_bld1","h2block48h_bld1","h2blockmth_bld1",
                  "antibio_12m_fec12","chm_12m_fec12","ims_12m_fec12","colsc_2m_fec12","ctscan_2m_fec12","cdiarr_2m_fec12",
                  "adiarr_2m_fec12","probio_2m_fec12","acid_2m_fec12","bile_2m_fec12","yog_2m_fec12","meatpref_fec12",
                  "stooltype_fec1","stooltype_fec1.1","stooltype_fec1.2","stooltype_fec1.3","stooltype_fec1.4","stooltype_fec1.5","stooltype_fec1.6",
                  "adj_choline1","adj_tmao1","adj_carnitine1","tc_plasma1","hdlc_plasma1","tg_plasma1","adj_crp1","tc.hdl1","hba1cp1",
                  "choline_std_bld1","tmao_std_bld1","carnitine_std_bld1",
                  "tc_std_bld1","hdlc_std_bld1","tg_std_bld1","tc.hdl_std_bld1","crp_std_bld1","hba1c_std_bld1")]

Basic2 = basic[,c("id","datem_bld1","dated_bld1","datem_fec2","dated_fec2","time_fec2","fast_bld1",
                  "smoke_bld1","bmi_bld1","totMETs_paq1","lt_fiber","lt_whgrn","lt_ahei","lt_ahei_na","lt_ahei_nowgr","lt_readmt","lt_nut","lt_veg", "lt_fru","lt_trypt",
                  "steroid48h_bld1","steroidmth_bld1","antibio48h_bld1","antibiomth_bld1","aspirin48h_bld1","aspirinmth_bld1",
                  "acetamp48h_bld1","acetampmth_bld1","ibuprof48h_bld1","ibuprofmth_bld1","h2block48h_bld1","h2blockmth_bld1",
                  "antibio_12m_fec12","chm_12m_fec12","ims_12m_fec12","colsc_2m_fec12","ctscan_2m_fec12","cdiarr_2m_fec12",
                  "adiarr_2m_fec12","probio_2m_fec12","acid_2m_fec12","bile_2m_fec12","yog_2m_fec12","meatpref_fec12",
                  "stooltype_fec2","stooltype_fec2.1","stooltype_fec2.2","stooltype_fec2.3","stooltype_fec2.4","stooltype_fec2.5","stooltype_fec2.6",
                  "adj_choline1","adj_tmao1","adj_carnitine1","tc_plasma1","hdlc_plasma1","tg_plasma1","adj_crp1","tc.hdl1","hba1cp1",
                  "choline_std_bld1","tmao_std_bld1","carnitine_std_bld1",
                  "tc_std_bld1","hdlc_std_bld1","tg_std_bld1","tc.hdl_std_bld1","crp_std_bld1","hba1c_std_bld1")]


Basic3 = basic[,c("id","datem_bld2","dated_bld2","datem_fec3","dated_fec3","time_fec3","fast_bld2",
                  "smoke_bld2","bmi_bld2","totMETs_paq2","lt_fiber","lt_whgrn","lt_ahei","lt_ahei_na","lt_ahei_nowgr","lt_readmt","lt_nut","lt_veg", "lt_fru","lt_trypt",
                  "steroid48h_bld2","steroidmth_bld2","antibio48h_bld2","antibiomth_bld2","aspirin48h_bld2","aspirinmth_bld2",
                  "acetamp48h_bld2","acetampmth_bld2","ibuprof48h_bld2","ibuprofmth_bld2","h2block48h_bld2","h2blockmth_bld2",
                  "antibio_12m_fec34","chm_12m_fec34","ims_12m_fec34","colsc_2m_fec34","ctscan_2m_fec34","cdiarr_2m_fec34",
                  "adiarr_2m_fec34","probio_2m_fec34","acid_2m_fec34","bile_2m_fec34","yog_2m_fec34","meatpref_fec34",
                  "stooltype_fec3","stooltype_fec3.1","stooltype_fec3.2","stooltype_fec3.3","stooltype_fec3.4","stooltype_fec3.5","stooltype_fec3.6",
                  "adj_choline2","adj_tmao2","adj_carnitine2","tc_plasma2","hdlc_plasma2","tg_plasma2","adj_crp2","tc.hdl2","hba1cp2",
                  "choline_std_bld2","tmao_std_bld2","carnitine_std_bld2",
                  "tc_std_bld2","hdlc_std_bld2","tg_std_bld2","tc.hdl_std_bld2","crp_std_bld2","hba1c_std_bld2")]

Basic4 = basic[,c("id","datem_bld2","dated_bld2","datem_fec4","dated_fec4","time_fec4","fast_bld2",
                  "smoke_bld2","bmi_bld2","totMETs_paq2","lt_fiber","lt_whgrn","lt_ahei","lt_ahei_na","lt_ahei_nowgr","lt_readmt","lt_nut","lt_veg", "lt_fru","lt_trypt",
                  "steroid48h_bld2","steroidmth_bld2","antibio48h_bld2","antibiomth_bld2","aspirin48h_bld2","aspirinmth_bld2",
                  "acetamp48h_bld2","acetampmth_bld2","ibuprof48h_bld2","ibuprofmth_bld2","h2block48h_bld2","h2blockmth_bld2",
                  "antibio_12m_fec34","chm_12m_fec34","ims_12m_fec34","colsc_2m_fec34","ctscan_2m_fec34","cdiarr_2m_fec34",
                  "adiarr_2m_fec34","probio_2m_fec34","acid_2m_fec34","bile_2m_fec34","yog_2m_fec34","meatpref_fec34",
                  "stooltype_fec4","stooltype_fec4.1","stooltype_fec4.2","stooltype_fec4.3","stooltype_fec4.4","stooltype_fec4.5","stooltype_fec4.6",
                  "adj_choline2","adj_tmao2","adj_carnitine2","tc_plasma2","hdlc_plasma2","tg_plasma2","adj_crp2","tc.hdl2","hba1cp2",
                  "choline_std_bld2","tmao_std_bld2","carnitine_std_bld2",
                  "tc_std_bld2","hdlc_std_bld2","tg_std_bld2","tc.hdl_std_bld2","crp_std_bld2","hba1c_std_bld2")]


newnames = c("id","datem_bld","dated_bld","datem_fec","dated_fec","time_fec","fast_bld",
             "smoke_bld","bmi_bld","totMETs_paq","lt_fiber","lt_whgrn","lt_ahei","lt_ahei_na","lt_ahei_nowgr","lt_readmt","lt_nut","lt_veg", "lt_fru","lt_trypt",
             "steroid48h_bld","steroidmth_bld","antibio48h_bld","antibiomth_bld","aspirin48h_bld","aspirinmth_bld",
             "acetamp48h_bld","acetampmth_bld","ibuprof48h_bld","ibuprofmth_bld","h2block48h_bld","h2blockmth_bld",
             "antibio_12m_fec","chm_12m_fec","ims_12m_fec","colsc_2m_fec","ctscan_2m_fec","cdiarr_2m_fec",
             "adiarr_2m_fec","probio_2m_fec","acid_2m_fec","bile_2m_fec","yog_2m_fec","meatpref_fec",
             "stooltype_fec","stooltype_fec.1","stooltype_fec.2","stooltype_fec.3","stooltype_fec.4","stooltype_fec.5","stooltype_fec.6",
             "adj_choline","adj_tmao","adj_carnitine","tc_plasma","hdlc_plasma","tg_plasma","adj_crp","tc.hdl","hba1cp",
             "choline_std_bld","tmao_std_bld","carnitine_std_bld",
             "tc_std_bld","hdlc_std_bld","tg_std_bld","tc.hdl_std_bld","crp_std_bld","hba1c_std_bld")

names(Basic1) = newnames
names(Basic2) = newnames
names(Basic3) = newnames
names(Basic4) = newnames

dim(Basic1) # 415 67
dim(Basic2) # 415 67
dim(Basic3) # 415 67
dim(Basic4) # 415 67


# ---------------------------------------------------------------------------;
#
#                        Processing the FFQ data
#
# ---------------------------------------------------------------------------;

ffq = read.csv("/udd/n2gji/data/FFQD.csv",header=T)
#head(ffq)
dim(ffq) # 727 616

names(ffq)

#------------
# basic adjustment
# total calories

summary(ffq$calor_fs_ffq1)
#   Min. 1st Qu.  Median    Mean 3rd Qu.    Max.    NA's 
#    864    1764    2128    2163    2479    3961      15
summary(ffq$calor_fs_ffq2)
#   Min. 1st Qu.  Median    Mean 3rd Qu.    Max.    NA's 
#    859    1711    2048    2116    2470    4162      75 

#------------
# basic adjustment - alcohol

summary(ffq$alco_fs_ffq1)
#   Min. 1st Qu.  Median    Mean 3rd Qu.    Max.    NA's 
#   0.00    3.37   13.94   17.73   23.93  118.41      15
summary(ffq$alco_fs_ffq2)
#   Min. 1st Qu.  Median    Mean 3rd Qu.    Max.    NA's 
#   0.00    3.57   13.68   16.65   22.23  112.99      75


# ---------------------------------------------------------------------------;
#                        Get FFQ file to use
# ---------------------------------------------------------------------------;

FFQ1 = ffq[,c("id","calor_fs_ffq1","alco_fs_ffq1","wholegrain_ffq1","fruit_ffq1","vege_ffq1","nut_ffq1","legume_ffq1","ssb_ffq1","dairy_ffq1","egg_ffq1","fish_ffq1",
              "meat_ffq1","livers_ffq1","meatliv_ffq1","redmt_ffq1","redliv_ffq1","redmtliv_ffq1","ptrymt_ffq1","ptryliv_ffq1","ptrymtliv_ffq1",
              "whgrn_fs_ffq1","whsub_fs_ffq1","crude_fs_ffq1","carbo_fs_ffq1","aofib_fs_ffq1","engl_fs_ffq1",
              "prot_fs_ffq1","aprot_fs_ffq1","vprot_fs_ffq1","tfat_fs_ffq1","afat_fs_ffq1","vfat_fs_ffq1",
              "betaine_fs_ffq1","frcho_fs_ffq1","gpcho_fs_ffq1","pcho_fs_ffq1","ptdcho_fs_ffq1","sphingo_fs_ffq1","choline_fs_ffq1","betchol_fs_ffq1",
              "betaine_fo_ffq1","frcho_fo_ffq1","ptdcho_fo_ffq1","choline_fo_ffq1","betchol_fo_ffq1",
              "b6_fs_ffq1","fol98_fs_ffq1","b12_fs_ffq1","b6_fo_ffq1","fol98_fo_ffq1","b12_fo_ffq1")]

FFQ2 = FFQ1

FFQ3 = ffq[,c("id","calor_fs_ffq2","alco_fs_ffq2","wholegrain_ffq2","fruit_ffq2","vege_ffq2","nut_ffq2","legume_ffq2","ssb_ffq2","dairy_ffq2","egg_ffq2","fish_ffq2",
              "meat_ffq2","livers_ffq2","meatliv_ffq2","redmt_ffq2","redliv_ffq2","redmtliv_ffq2","ptrymt_ffq2","ptryliv_ffq2","ptrymtliv_ffq2",
              "whgrn_fs_ffq2","whsub_fs_ffq2","crude_fs_ffq2","carbo_fs_ffq2","aofib_fs_ffq2","engl_fs_ffq2",
              "prot_fs_ffq2","aprot_fs_ffq2","vprot_fs_ffq2","tfat_fs_ffq2","afat_fs_ffq2","vfat_fs_ffq2",
              "betaine_fs_ffq2","frcho_fs_ffq2","gpcho_fs_ffq2","pcho_fs_ffq2","ptdcho_fs_ffq2","sphingo_fs_ffq2","choline_fs_ffq2","betchol_fs_ffq2",
              "betaine_fo_ffq2","frcho_fo_ffq2","ptdcho_fo_ffq2","choline_fo_ffq2","betchol_fo_ffq2",
              "b6_fs_ffq2","fol98_fs_ffq2","b12_fs_ffq2","b6_fo_ffq2","fol98_fo_ffq2","b12_fo_ffq2")]

FFQ4 = FFQ3

namesffq = c("id","calor_fs_ffq","alco_fs_ffq","wholegrain_ffq","fruit_ffq","vege_ffq","nut_ffq","legume_ffq","ssb_ffq","dairy_ffq","egg_ffq","fish_ffq",
             "meat_ffq","livers_ffq","meatliv_ffq","redmt_ffq","redliv_ffq","redmtliv_ffq","ptrymt_ffq","ptryliv_ffq","ptrymtliv_ffq",
             "whgrn_fs_ffq","whsub_fs_ffq","crude_fs_ffq","carbo_fs_ffq","aofib_fs_ffq","engl_fs_ffq",
             "prot_fs_ffq","aprot_fs_ffq","vprot_fs_ffq","tfat_fs_ffq","afat_fs_ffq","vfat_fs_ffq",
             "betaine_fs_ffq","frcho_fs_ffq","gpcho_fs_ffq","pcho_fs_ffq","ptdcho_fs_ffq","sphingo_fs_ffq","choline_fs_ffq","betchol_fs_ffq",
             "betaine_fo_ffq","frcho_fo_ffq","ptdcho_fo_ffq","choline_fo_ffq","betchol_fo_ffq",
             "b6_fs_ffq","fol98_fs_ffq","b12_fs_ffq","b6_fo_ffq","fol98_fo_ffq","b12_fo_ffq")


identical( substr(names(FFQ1)[-1],1,nchar(names(FFQ1)[-1])-1), 
           substr(names(FFQ3)[-1],1,nchar(names(FFQ3)[-1])-1) )

identical( substr(names(FFQ1)[-1],1,nchar(names(FFQ1)[-1])-1), 
           namesffq[-1] )

names(FFQ1) = namesffq
names(FFQ2) = namesffq
names(FFQ3) = namesffq
names(FFQ4) = namesffq

dim(FFQ1) # 727 52
dim(FFQ2) # 727 52
dim(FFQ3) # 727 52
dim(FFQ4) # 727 52



# ---------------------------------------------------------------------------;
#
#                        Processing the 7DDR data
#
# ---------------------------------------------------------------------------;

# for this analysis, we will use the weekly everage 

ddr = read.csv("/udd/n2gji/data/DDRD.csv",header=T)
ddr_hy = read.csv("/udd/hpyah/project/microb_metabo/Data/DDRD.csv",header=T)
times = c("ddr_date_w1d1","ddr_date_w1d2","ddr_date_w1d3","ddr_date_w1d4","ddr_date_w1d5","ddr_date_w1d6","ddr_date_w1d7",
          "ddr_date_w2d1","ddr_date_w2d2","ddr_date_w2d3","ddr_date_w2d4","ddr_date_w2d5","ddr_date_w2d6","ddr_date_w2d7")
ddr_hy_subset <- ddr_hy[, c("id", times)]

# 根据 id 合并
ddr <- merge(ddr, ddr_hy_subset, by = "id", all.x = TRUE)

# 检查结果
head(ddr)
#head(ddr)
dim(ddr) #  692 1769

names(ddr)

# ---------------------------------------------------------------------------;
# merge with fecal sample info to calculate corresponding records
# ---------------------------------------------------------------------------;

# facal time

facaltime = basic[,c("id","dated_fec1","dated_fec2","dated_fec3","dated_fec4")]

# use id

lname = load(file="/udd/hpyah/project/microb_metabo/Data/ReadIn_taxon.RData")
facal1 = merge(facaltime[,c("id","dated_fec1")],data.frame(id=Stands_s1[,"id"]),by="id")
facal2 = merge(facaltime[,c("id","dated_fec2")],data.frame(id=Stands_s2[,"id"]),by="id")
facal3 = merge(facaltime[,c("id","dated_fec3")],data.frame(id=Stands_s3[,"id"]),by="id")
facal4 = merge(facaltime[,c("id","dated_fec4")],data.frame(id=Stands_s4[,"id"]),by="id")

#head(facal1)
#head(facal2)
#head(facal3)
#head(facal4)

dim(facal1) # 210 2
dim(facal2) # 205 2
dim(facal3) # 256 2
dim(facal4) # 254 2


#--------------------
#--------------------
#--------------------
# process time 1
#--------------------
#--------------------
#--------------------

DDR1 = merge(facal1,ddr,by="id",all.x=T)
dim(DDR1)#210  1770
names(ddr)
# time var - 
# use all ddr time points to search for days related to facal collections

times = c("ddr_date_w1d1","ddr_date_w1d2","ddr_date_w1d3","ddr_date_w1d4","ddr_date_w1d5","ddr_date_w1d6","ddr_date_w1d7",
          "ddr_date_w2d1","ddr_date_w2d2","ddr_date_w2d3","ddr_date_w2d4","ddr_date_w2d5","ddr_date_w2d6","ddr_date_w2d7")

for(i in 1:length(times)) {
  DDR1$new = dateproc(times[i],DDR1,"days")
  names(DDR1)[which(names(DDR1)=="new")] = paste("ddr_dated_",substr(times[i],10,13),sep='')
}


# ----------------------------------- NOTE -----------------------------------
# calcualte average for the 5 days before and at the day of sample collection
# ----------------------------------------------------------------------------

fix = c("w1d1","w1d2","w1d3","w1d4","w1d5","w1d6","w1d7",
        "w2d1","w2d2","w2d3","w2d4","w2d5","w2d6","w2d7")

times = paste("ddr_dated_",fix,sep='')
calor = paste("calor_fo_dr_",fix,sep='')
fat = paste("fat_fo_dr_",fix,sep='')
carbo = paste("carbo_fo_dr_",fix,sep='')
prot = paste("prot_fo_dr_",fix,sep='')

aprot = paste("aprot_fo_dr_",fix,sep='')
vprot = paste("vprot_fo_dr_",fix,sep='')
alco = paste("alco_fo_dr_",fix,sep='')
satfat = paste("satfat_fo_dr_",fix,sep='')
monfat = paste("monfat_fo_dr_",fix,sep='')
poly = paste("poly_fo_dr_",fix,sep='')
aofib = paste("aofib_fo_dr_",fix,sep='')
wsdf = paste("wsdf_fo_dr_",fix,sep='')
ifib = paste("ifib_fo_dr_",fix,sep='')

choline_fs = paste("choline_fs_dr_",fix,sep='')
choline_fo = paste("choline_fo_dr_",fix,sep='')
choline_so = paste("choline_so_dr_",fix,sep='')
betaine_fs = paste("betaine_fs_dr_",fix,sep='')
betaine_fo = paste("betaine_fo_dr_",fix,sep='')
betaine_so = paste("betaine_so_dr_",fix,sep='')

trypt_fs = paste("trypto_fs_dr_",fix,sep='')
trypt_fo = paste("trypto_fo_dr_",fix,sep='')
trypt_so = paste("trypto_so_dr_",fix,sep='')

wholegrain = paste("wholegrain_",fix,sep='')
wholegrain_main = paste("wholegrain_main_",fix,sep='')
wholegrain_snack = paste("wholegrain_snack_",fix,sep='')
wholegrain_cake = paste("wholegrain_cake_",fix,sep='')
wgr = paste("wgr_",fix,sep='')
cer = paste("cer_",fix,sep='')
oat = paste("oat_",fix,sep='')
ryebr = paste("ryebr_",fix,sep='')
dkbr = paste("dkbr_",fix,sep='')
brice = paste("brice_",fix,sep='')
ffpop = paste("ffpop_",fix,sep='')
popc = paste("popc_",fix,sep='')
pop = paste("pop_",fix,sep='')

artfood=paste("artfood_",fix,sep='')

redmeat = paste("redmeat_",fix,sep='')
liver = paste("liver_",fix,sep='')
poultry = paste("poultry_",fix,sep='')
fish = paste("fish_",fix,sep='')
egg = paste("egg_",fix,sep='')
wdairy = paste("wdairy_",fix,sep='')
hfdairy = paste("hfdairy_",fix,sep='')
fruits = paste("fruits_",fix,sep='')
veg = paste("veg_",fix,sep='')
leg = paste("leg_",fix,sep='')
nuts = paste("nuts_",fix,sep='')
pnut07= paste("pnut07_wf_dr_",fix,sep='')
pbut07= paste("pbut07_wf_dr_",fix,sep='')
wnut07= paste("wnut07_wf_dr_",fix,sep='')
pwnut07= paste("pwnut07_wf_dr_",fix,sep='')
onut07= paste("onut07_wf_dr_",fix,sep='')
veg07= paste("veg07_wf_dr_",fix,sep='')
leg07= paste("leg07_wf_dr_",fix,sep='')
nut07= paste("nut07_wf_dr_",fix,sep='')
fruit07= paste("fruit07_wf_dr_",fix,sep='')
rmeat07= paste("rmeat07_wf_dr_",fix,sep='')
chicken07= paste("chicken07_wf_dr_",fix,sep='')
fish07= paste("fish07_wf_dr_",fix,sep='')
seafood07= paste("seafood07_wf_dr_",fix,sep='')
egg07= paste("egg07_wf_dr_",fix,sep='')
dairy07= paste("dairy07_wf_dr_",fix,sep='')
######################## For WG analysis, doesn't have to do the exact 7DDR-Fecal match ########################
######################## Weekly average is sufficient and make most use of data ########################
# ds=c(seq(1,7,1))
# dss=1
# for(i in 2:7){
#   dss = paste(dss,",",ds[i],sep='')
# }
# 
# DDR1$use_d = NA
# for(i in 1:dim(DDR1)[1]) {
#   DDR1$use_d[i] = dss
# }

DDR1$use_d = NA

var_use=c("calor_fo_dr_","fat_fo_dr_","carbo_fo_dr_","prot_fo_dr_",
          "aprot_fo_dr_","vprot_fo_dr_","alco_fo_dr_","satfat_fo_dr_","monfat_fo_dr_",
          "poly_fo_dr_","aofib_fo_dr_","wsdf_fo_dr_","ifib_fo_dr_","choline_fs_dr_",
          "choline_fo_dr_","choline_so_dr_","betaine_fs_dr_","betaine_fo_dr_","betaine_so_dr_", "trypto_fs_dr_","trypto_fo_dr_","trypto_so_dr_",
          "wgr_","cer_","oat_","ryebr_","dkbr_","brice_","ffpop_","popc_","pop_",
          "pnut07_wf_dr_","pbut07_wf_dr_","wnut07_wf_dr_","pwnut07_wf_dr_","onut07_wf_dr_",
          "veg07_wf_dr_","leg07_wf_dr_","nut07_wf_dr_","fruit07_wf_dr_","rmeat07_wf_dr_",
          "chicken07_wf_dr_","fish07_wf_dr_","seafood07_wf_dr_","egg07_wf_dr_","dairy07_wf_dr_")
var_out=c("calor_fo_","fat_fo_","carbo_fo_","prot_fo_",
          "aprot_fo_","vprot_fo_","alco_fo_","satfat_fo_","monfat_fo_",
          "poly_fo_","aofib_fo_","wsdf_fo_","ifib_fo_","choline_fs_",
          "choline_fo_","choline_so_","betaine_fs_","betaine_fo_","betaine_so_","trypto_fs_","trypto_fo_","trypto_so_",
          "wgr_","cer_","oat_","ryebr_","dkbr_","brice_","ffpop_","popc_","pop_",
          "pnut07_","pbut07_","wnut07_","pwnut07_","onut07_",
          "veg07_","leg07_","nut07_","fruit07_","rmeat07_",
          "chicken07_","fish07_","seafood07_","egg07_","dairy07_")
names(DDR1)
d_use=c(
  "1d",
  "2d",
  "3d",
  "4d",
  "5d",
  "6d",
  "7d",
  "all")

for(i in 1:dim(DDR1)[1]) {
  
  t_rec = DDR1[i,times]
  t_fac = DDR1[i,"dated_fec1"]
  
  # days to use 
  d1 = which(t_fac-t_rec<=0 & t_fac-t_rec>=0)
  d2 = which(t_fac-t_rec<=1 & t_fac-t_rec>=0)
  d3 = which(t_fac-t_rec<=2 & t_fac-t_rec>=0)
  d4 = which(t_fac-t_rec<=3 & t_fac-t_rec>=0)
  d5 = which(t_fac-t_rec<=4 & t_fac-t_rec>=0)
  d6 = which(t_fac-t_rec<=5 & t_fac-t_rec>=0)
  d7 = which(t_fac-t_rec<=6 & t_fac-t_rec>=0)
  
  all_d<-seq(1,7,1)
  
  d=list(
    d1,
    d2,
    d3,
    d4,
    d5,
    d6,
    d7,
    all_d)
  
  for(j in 1:length(d_use)) {
    for(k in 1:length(var_use)){
      v_temp=paste(var_out[k],"avg",d_use[j],sep='')
      DDR1[i,v_temp] = mean(as.numeric(DDR1[i,paste(var_use[k],fix,sep='')[unlist(d[j])]]),na.rm=TRUE)
      if (is.na(DDR1[i,v_temp])){
        DDR1[i,v_temp] = mean(as.numeric(DDR1[i,paste(var_use[k],fix,sep='')[1:7]]),na.rm=TRUE)
      }
    }
  }
}

# DDR1_nomatch=data.frame(
#   rbind(
#     c(length(which(is.na(DDR1$calor_fo_avg1d))),
#       length(which(is.na(DDR1$calor_fo_avg2d))),
#       length(which(is.na(DDR1$calor_fo_avg3d))),
#       length(which(is.na(DDR1$calor_fo_avg4d))),
#       length(which(is.na(DDR1$calor_fo_avg5d))),
#       length(which(is.na(DDR1$calor_fo_avg6d))),
#       length(which(is.na(DDR1$calor_fo_avg7d)))),
#     c(length(which(is.na(DDR1$calor_fo_avg1d))),
#       length(which(is.na(DDR1$calor_fo_avg2d))),
#       length(which(is.na(DDR1$calor_fo_avg3d))),
#       length(which(is.na(DDR1$calor_fo_avg4d))),
#       length(which(is.na(DDR1$calor_fo_avg5d))),
#       length(which(is.na(DDR1$calor_fo_avg6d))),
#       length(which(is.na(DDR1$calor_fo_avg7d))))/dim(DDR1)[1]))
# colnames(DDR1_nomatch)=c("1d","2d","3d","4d","5d","6d","7d")
# rownames(DDR1_nomatch)=c("No match days","No match proportions")
# 
# DDR1_nomatch

x=NA
for(i in 1:length(var_out)){
  x=c(x,paste(paste(var_out,"avg",sep="")[i],d_use,sep=""))
}
x=x[-1]

x
# get the vars we wanted to use

DDR1 = DDR1[,c("id",times,calor,fat,carbo,prot,aprot,vprot,alco,satfat,monfat,poly,aofib,wsdf,ifib,
               choline_fs,choline_fo,choline_so,betaine_fs,betaine_fo,betaine_so,wgr,cer,oat,ryebr,dkbr,brice,ffpop,popc,
               "use_d",x,
               "aofib_fo_dr_wtavg","aofib_so_dr_wtavg","aofib_fs_dr_wtavg",
               "wsdf_fo_dr_wtavg","wsdf_so_dr_wtavg","wsdf_fs_dr_wtavg", 
               "ifib_fo_dr_wtavg","ifib_so_dr_wtavg","ifib_fs_dr_wtavg", 
               "trypto_fo_dr_wtavg", "trypto_so_dr_wtavg", "trypto_fs_dr_wtavg",
               "pect_fo_dr_wtavg",
               "wgrain07_wf_dr_wt",
               "pop07_wf_dr_wt",
               "veg07_wf_dr_wt",
               "leg07_wf_dr_wt",
               "nut07_wf_dr_wt",
               "pnut07_wf_dr_wt",
               "pbut07_wf_dr_wt",
               "wnut07_wf_dr_wt",
               "pwnut07_wf_dr_wt",
               "onut07_wf_dr_wt",
               "fruit07_wf_dr_wt",
               "rmeat07_wf_dr_wt",
               "pmeat07_wf_dr_wt",	
               "ameat07_wf_dr_wt",
               "chicken07_wf_dr_wt",
               "fish07_wf_dr_wt",
               "seafood07_wf_dr_wt",
               "milk07_wf_dr_wt",
               "yog07_wf_dr_wt",
               "egg07_wf_dr_wt",
               "calor_fo_dr_wtavg",
               "calor_so_dr_wtavg",
               "calor_fs_dr_wtavg")]

dim(DDR1) # 210 550



#--------------------
#--------------------
#--------------------
# process time 2
#--------------------
#--------------------
#--------------------


DDR2 = merge(facal2,ddr,by="id",all.x=T)
dim(DDR2) # 205 1170

# time var - 
# use all ddr time points to search for days related to facal collections

times = c("ddr_date_w1d1","ddr_date_w1d2","ddr_date_w1d3","ddr_date_w1d4","ddr_date_w1d5","ddr_date_w1d6","ddr_date_w1d7",
          "ddr_date_w2d1","ddr_date_w2d2","ddr_date_w2d3","ddr_date_w2d4","ddr_date_w2d5","ddr_date_w2d6","ddr_date_w2d7")

for(i in 1:length(times)) {
  DDR2$new = dateproc(times[i],DDR2,"days")
  names(DDR2)[which(names(DDR2)=="new")] = paste("ddr_dated_",substr(times[i],10,13),sep='')
}


# ----------------------------------- NOTE -----------------------------------
# calcualte average for the 5 days before and at the day of sample collection
# for those without matched dates just use NA
# ----------------------------------------------------------------------------

fix = c("w1d1","w1d2","w1d3","w1d4","w1d5","w1d6","w1d7",
        "w2d1","w2d2","w2d3","w2d4","w2d5","w2d6","w2d7")
times = paste("ddr_dated_",fix,sep='')

# ds=c(seq(1,7,1))
# dss=1
# for(i in 2:7){
#   dss = paste(dss,",",ds[i],sep='')
# }
# 
# DDR2$use_d = NA
# for(i in 1:dim(DDR2)[1]) {
#   DDR2$use_d[i] = dss
# }

DDR2$use_d = NA

for(i in 1:dim(DDR2)[1]) {
  
  t_rec = DDR2[i,times]
  t_fac = DDR2[i,"dated_fec2"]
  
  # days to use 
  d1 = which(t_fac-t_rec<=0 & t_fac-t_rec>=0)
  d2 = which(t_fac-t_rec<=1 & t_fac-t_rec>=0)
  d3 = which(t_fac-t_rec<=2 & t_fac-t_rec>=0)
  d4 = which(t_fac-t_rec<=3 & t_fac-t_rec>=0)
  d5 = which(t_fac-t_rec<=4 & t_fac-t_rec>=0)
  d6 = which(t_fac-t_rec<=5 & t_fac-t_rec>=0)
  d7 = which(t_fac-t_rec<=6 & t_fac-t_rec>=0)
  
  all_d<-seq(1,7,1)
  
  d=list(
    d1,
    d2,
    d3,
    d4,
    d5,
    d6,
    d7,
    all_d)
  
  for(j in 1:length(d_use)) {
    for(k in 1:length(var_use)){
      v_temp=paste(var_out[k],"avg",d_use[j],sep='')
      DDR2[i,v_temp] = mean(as.numeric(DDR2[i,paste(var_use[k],fix,sep='')[unlist(d[j])]]),na.rm=TRUE)
      if (is.na(DDR2[i,v_temp])){
        DDR2[i,v_temp] = mean(as.numeric(DDR2[i,paste(var_use[k],fix,sep='')[1:7]]),na.rm=TRUE)
      }
    }
  }
}

# DDR2_nomatch=data.frame(
#   rbind(
#     c(length(which(is.na(DDR2$calor_fo_avg1d))),
#       length(which(is.na(DDR2$calor_fo_avg2d))),
#       length(which(is.na(DDR2$calor_fo_avg3d))),
#       length(which(is.na(DDR2$calor_fo_avg4d))),
#       length(which(is.na(DDR2$calor_fo_avg5d))),
#       length(which(is.na(DDR2$calor_fo_avg6d))),
#       length(which(is.na(DDR2$calor_fo_avg7d)))),
#     c(length(which(is.na(DDR2$calor_fo_avg1d))),
#       length(which(is.na(DDR2$calor_fo_avg2d))),
#       length(which(is.na(DDR2$calor_fo_avg3d))),
#       length(which(is.na(DDR2$calor_fo_avg4d))),
#       length(which(is.na(DDR2$calor_fo_avg5d))),
#       length(which(is.na(DDR2$calor_fo_avg6d))),
#       length(which(is.na(DDR2$calor_fo_avg7d))))/dim(DDR2)[1]))
# colnames(DDR2_nomatch)=c("1d","2d","3d","4d","5d","6d","7d")
# rownames(DDR2_nomatch)=c("No match days","No match proportions")
# 
# DDR2_nomatch

# get the vars we wanted to use

DDR2 = DDR2[,c("id",times,calor,fat,carbo,prot,aprot,vprot,alco,satfat,monfat,poly,aofib,wsdf,ifib,
               choline_fs,choline_fo,choline_so,betaine_fs,betaine_fo,betaine_so,wgr,cer,oat,ryebr,dkbr,brice,ffpop,popc,
               "use_d",x,
               "aofib_fo_dr_wtavg","aofib_so_dr_wtavg","aofib_fs_dr_wtavg",
               "wsdf_fo_dr_wtavg","wsdf_so_dr_wtavg","wsdf_fs_dr_wtavg", 
               "ifib_fo_dr_wtavg","ifib_so_dr_wtavg","ifib_fs_dr_wtavg", 
               "trypto_fo_dr_wtavg", "trypto_so_dr_wtavg", "trypto_fs_dr_wtavg",
               "pect_fo_dr_wtavg",
               "wgrain07_wf_dr_wt",
               "pop07_wf_dr_wt",
               "veg07_wf_dr_wt",
               "leg07_wf_dr_wt",
               "nut07_wf_dr_wt",
               "pnut07_wf_dr_wt",
               "pbut07_wf_dr_wt",
               "wnut07_wf_dr_wt",
               "pwnut07_wf_dr_wt",
               "onut07_wf_dr_wt",
               "fruit07_wf_dr_wt",
               "rmeat07_wf_dr_wt",
               "pmeat07_wf_dr_wt",	
               "ameat07_wf_dr_wt",
               "chicken07_wf_dr_wt",
               "fish07_wf_dr_wt",
               "seafood07_wf_dr_wt",
               "milk07_wf_dr_wt",
               "yog07_wf_dr_wt",
               "egg07_wf_dr_wt",
               "calor_fo_dr_wtavg",
               "calor_so_dr_wtavg",
               "calor_fs_dr_wtavg")]

dim(DDR2) # 205 550



#--------------------
#--------------------
#--------------------
# process time 3
#--------------------
#--------------------
#--------------------

DDR3 = merge(facal3,ddr,by="id",all.x=T)
dim(DDR3) # 205 1770

# time var - 
# use all ddr time points to search for days related to facal collections

times = c("ddr_date_w1d1","ddr_date_w1d2","ddr_date_w1d3","ddr_date_w1d4","ddr_date_w1d5","ddr_date_w1d6","ddr_date_w1d7",
          "ddr_date_w2d1","ddr_date_w2d2","ddr_date_w2d3","ddr_date_w2d4","ddr_date_w2d5","ddr_date_w2d6","ddr_date_w2d7")

for(i in 1:length(times)) {
  DDR3$new = dateproc(times[i],DDR3,"days")
  names(DDR3)[which(names(DDR3)=="new")] = paste("ddr_dated_",substr(times[i],10,13),sep='')
}


# ----------------------------------- NOTE -----------------------------------
# calcualte average for the 5 days before and at the day of sample collection
# for those without matched dates just use NA in this step and process later in the second step
# ----------------------------------------------------------------------------


fix = c("w1d1","w1d2","w1d3","w1d4","w1d5","w1d6","w1d7",
        "w2d1","w2d2","w2d3","w2d4","w2d5","w2d6","w2d7")
times = paste("ddr_dated_",fix,sep='')

# ds=c(seq(1,7,1))
# dss=1
# for(i in 2:7){
#   dss = paste(dss,",",ds[i],sep='')
# }
# 
# DDR2$use_d = NA
# for(i in 1:dim(DDR2)[1]) {
#   DDR2$use_d[i] = dss
# }

DDR3$use_d = NA

for(i in 1:dim(DDR3)[1]) {
  
  t_rec = DDR3[i,times]
  t_fac = DDR3[i,"dated_fec3"]
  
  # days to use 
  d1 = which(t_fac-t_rec<=0 & t_fac-t_rec>=0)
  d2 = which(t_fac-t_rec<=1 & t_fac-t_rec>=0)
  d3 = which(t_fac-t_rec<=2 & t_fac-t_rec>=0)
  d4 = which(t_fac-t_rec<=3 & t_fac-t_rec>=0)
  d5 = which(t_fac-t_rec<=4 & t_fac-t_rec>=0)
  d6 = which(t_fac-t_rec<=5 & t_fac-t_rec>=0)
  d7 = which(t_fac-t_rec<=6 & t_fac-t_rec>=0)
  
  all_d<-seq(8,14,1)
  
  d=list(
    d1,
    d2,
    d3,
    d4,
    d5,
    d6,
    d7,
    all_d)
  
  for(j in 1:length(d_use)) {
    for(k in 1:length(var_use)){
      v_temp=paste(var_out[k],"avg",d_use[j],sep='')
      DDR3[i,v_temp] = mean(as.numeric(DDR3[i,paste(var_use[k],fix,sep='')[unlist(d[j])]]),na.rm=TRUE)
      if (is.na(DDR3[i,v_temp])){
        DDR3[i,v_temp] = mean(as.numeric(DDR3[i,paste(var_use[k],fix,sep='')[8:14]]),na.rm=TRUE)
      }
    }
  }
}

# DDR3_nomatch=data.frame(
#   rbind(
#     c(length(which(is.na(DDR3$calor_fo_avg1d))),
#       length(which(is.na(DDR3$calor_fo_avg2d))),
#       length(which(is.na(DDR3$calor_fo_avg3d))),
#       length(which(is.na(DDR3$calor_fo_avg4d))),
#       length(which(is.na(DDR3$calor_fo_avg5d))),
#       length(which(is.na(DDR3$calor_fo_avg6d))),
#       length(which(is.na(DDR3$calor_fo_avg7d)))),
#     c(length(which(is.na(DDR3$calor_fo_avg1d))),
#       length(which(is.na(DDR3$calor_fo_avg2d))),
#       length(which(is.na(DDR3$calor_fo_avg3d))),
#       length(which(is.na(DDR3$calor_fo_avg4d))),
#       length(which(is.na(DDR3$calor_fo_avg5d))),
#       length(which(is.na(DDR3$calor_fo_avg6d))),
#       length(which(is.na(DDR3$calor_fo_avg7d))))/dim(DDR3)[1]))
# colnames(DDR3_nomatch)=c("1d","2d","3d","4d","5d","6d","7d")
# rownames(DDR3_nomatch)=c("No match days","No match proportions")

#DDR3_nomatch


# get the vars we wanted to use

DDR3 = DDR3[,c("id",times,calor,fat,carbo,prot,aprot,vprot,alco,satfat,monfat,poly,aofib,wsdf,ifib,
               choline_fs,choline_fo,choline_so,betaine_fs,betaine_fo,betaine_so,wgr,cer,oat,ryebr,dkbr,brice,ffpop,popc,
               "use_d",x,
               "aofib_fo_dr_wtavg","aofib_so_dr_wtavg","aofib_fs_dr_wtavg",
               "wsdf_fo_dr_wtavg","wsdf_so_dr_wtavg","wsdf_fs_dr_wtavg", 
               "ifib_fo_dr_wtavg","ifib_so_dr_wtavg","ifib_fs_dr_wtavg", 
               "trypto_fo_dr_wtavg", "trypto_so_dr_wtavg", "trypto_fs_dr_wtavg",
               "pect_fo_dr_wtavg",
               "wgrain07_wf_dr_wt",
               "pop07_wf_dr_wt",
               "veg07_wf_dr_wt",
               "leg07_wf_dr_wt",
               "nut07_wf_dr_wt",
               "pnut07_wf_dr_wt",
               "pbut07_wf_dr_wt",
               "wnut07_wf_dr_wt",
               "pwnut07_wf_dr_wt",
               "onut07_wf_dr_wt",
               "fruit07_wf_dr_wt",
               "rmeat07_wf_dr_wt",
               "pmeat07_wf_dr_wt",	
               "ameat07_wf_dr_wt",
               "chicken07_wf_dr_wt",
               "fish07_wf_dr_wt",
               "seafood07_wf_dr_wt",
               "milk07_wf_dr_wt",
               "yog07_wf_dr_wt",
               "egg07_wf_dr_wt",
               "calor_fo_dr_wtavg",
               "calor_so_dr_wtavg",
               "calor_fs_dr_wtavg")]

dim(DDR3) # 256 550



#--------------------
#--------------------
#--------------------
# process time 4
#--------------------
#--------------------
#--------------------


DDR4 = merge(facal4,ddr,by="id",all.x=T)
dim(DDR4) # 205 447

# time var - 
# use all ddr time points to search for days related to facal collections

times = c("ddr_date_w1d1","ddr_date_w1d2","ddr_date_w1d3","ddr_date_w1d4","ddr_date_w1d5","ddr_date_w1d6","ddr_date_w1d7",
          "ddr_date_w2d1","ddr_date_w2d2","ddr_date_w2d3","ddr_date_w2d4","ddr_date_w2d5","ddr_date_w2d6","ddr_date_w2d7")

for(i in 1:length(times)) {
  DDR4$new = dateproc(times[i],DDR4,"days")
  names(DDR4)[which(names(DDR4)=="new")] = paste("ddr_dated_",substr(times[i],10,13),sep='')
}


# ----------------------------------- NOTE -----------------------------------
# calcualte average for the 5 days before and at the day of sample collection
# for those without matched dates just use NA in this step and process later in the second step
# ----------------------------------------------------------------------------

fix = c("w1d1","w1d2","w1d3","w1d4","w1d5","w1d6","w1d7",
        "w2d1","w2d2","w2d3","w2d4","w2d5","w2d6","w2d7")
times = paste("ddr_dated_",fix,sep='')

# ds=c(seq(1,7,1))
# dss=1
# for(i in 2:7){
#   dss = paste(dss,",",ds[i],sep='')
# }
# 
# DDR2$use_d = NA
# for(i in 1:dim(DDR2)[1]) {
#   DDR2$use_d[i] = dss
# }

DDR4$use_d = NA

for(i in 1:dim(DDR4)[1]) {
  
  t_rec = DDR4[i,times]
  t_fac = DDR4[i,"dated_fec4"]
  
  # days to use 
  d1 = which(t_fac-t_rec<=0 & t_fac-t_rec>=0)
  d2 = which(t_fac-t_rec<=1 & t_fac-t_rec>=0)
  d3 = which(t_fac-t_rec<=2 & t_fac-t_rec>=0)
  d4 = which(t_fac-t_rec<=3 & t_fac-t_rec>=0)
  d5 = which(t_fac-t_rec<=4 & t_fac-t_rec>=0)
  d6 = which(t_fac-t_rec<=5 & t_fac-t_rec>=0)
  d7 = which(t_fac-t_rec<=6 & t_fac-t_rec>=0)
  
  all_d<-seq(8,14,1)
  
  d=list(
    d1,
    d2,
    d3,
    d4,
    d5,
    d6,
    d7,
    all_d)
  
  for(j in 1:length(d_use)) {
    for(k in 1:length(var_use)){
      v_temp=paste(var_out[k],"avg",d_use[j],sep='')
      DDR4[i,v_temp] = mean(as.numeric(DDR4[i,paste(var_use[k],fix,sep='')[unlist(d[j])]]),na.rm=TRUE)
      if (is.na(DDR4[i,v_temp])){
        DDR4[i,v_temp] = mean(as.numeric(DDR4[i,paste(var_use[k],fix,sep='')[8:14]]),na.rm=TRUE)
      }
    }
  }
}

# DDR4_nomatch=data.frame(
#   rbind(
#     c(length(which(is.na(DDR4$calor_fo_avg1d))),
#       length(which(is.na(DDR4$calor_fo_avg2d))),
#       length(which(is.na(DDR4$calor_fo_avg3d))),
#       length(which(is.na(DDR4$calor_fo_avg4d))),
#       length(which(is.na(DDR4$calor_fo_avg5d))),
#       length(which(is.na(DDR4$calor_fo_avg6d))),
#       length(which(is.na(DDR4$calor_fo_avg7d)))),
#     c(length(which(is.na(DDR4$calor_fo_avg1d))),
#       length(which(is.na(DDR4$calor_fo_avg2d))),
#       length(which(is.na(DDR4$calor_fo_avg3d))),
#       length(which(is.na(DDR4$calor_fo_avg4d))),
#       length(which(is.na(DDR4$calor_fo_avg5d))),
#       length(which(is.na(DDR4$calor_fo_avg6d))),
#       length(which(is.na(DDR4$calor_fo_avg7d))))/dim(DDR4)[1]))
# colnames(DDR4_nomatch)=c("1d","2d","3d","4d","5d","6d","7d")
# rownames(DDR4_nomatch)=c("No match days","No match proportions")
# 
# DDR4_nomatch

# get the vars we wanted to use

DDR4 = DDR4[,c("id",times,calor,fat,carbo,prot,aprot,vprot,alco,satfat,monfat,poly,aofib,wsdf,ifib,
               choline_fs,choline_fo,choline_so,betaine_fs,betaine_fo,betaine_so,wgr,cer,oat,ryebr,dkbr,brice,ffpop,popc,
               "use_d",x,
               "aofib_fo_dr_wtavg","aofib_so_dr_wtavg","aofib_fs_dr_wtavg",
               "wsdf_fo_dr_wtavg","wsdf_so_dr_wtavg","wsdf_fs_dr_wtavg", 
               "ifib_fo_dr_wtavg","ifib_so_dr_wtavg","ifib_fs_dr_wtavg", 
               "trypto_fo_dr_wtavg", "trypto_so_dr_wtavg", "trypto_fs_dr_wtavg",
               "pect_fo_dr_wtavg",
               "wgrain07_wf_dr_wt",
               "pop07_wf_dr_wt",
               "veg07_wf_dr_wt",
               "leg07_wf_dr_wt",
               "nut07_wf_dr_wt",
               "pnut07_wf_dr_wt",
               "pbut07_wf_dr_wt",
               "wnut07_wf_dr_wt",
               "pwnut07_wf_dr_wt",
               "onut07_wf_dr_wt",
               "fruit07_wf_dr_wt",
               "rmeat07_wf_dr_wt",
               "pmeat07_wf_dr_wt",	
               "ameat07_wf_dr_wt",
               "chicken07_wf_dr_wt",
               "fish07_wf_dr_wt",
               "seafood07_wf_dr_wt",
               "milk07_wf_dr_wt",
               "yog07_wf_dr_wt",
               "egg07_wf_dr_wt",
               "calor_fo_dr_wtavg",
               "calor_so_dr_wtavg",
               "calor_fs_dr_wtavg")]

dim(DDR4) # 254 550




# ---------------------------------------------------------------------------;
#
#                   Load microbiome information for merging
#
# ---------------------------------------------------------------------------;
# for taxon -- use standardized data for analysis
#-----------------------------------------------------------------------------

lname = load(file="/udd/hpxil/3.source/hpfs_microbiome_biobakery4/ReadIn_taxon.RData")
print(lname)
#"annot_taxon" "taxon_total" "Stands_s1"   "Stands_s2"   "Stands_s3"   "Stands_s4" 
# stands - sum is 100 - standize to 1 by /100

Stands_s1[,-1] = Stands_s1[,-1]/100
Stands_s2[,-1] = Stands_s2[,-1]/100
Stands_s3[,-1] = Stands_s3[,-1]/100
Stands_s4[,-1] = Stands_s4[,-1]/100
dim(Stands_s1)
dim(Stands_s2)
dim(Stands_s3)
dim(Stands_s4)

Stands_s1[1:3,1:20]
Stands_s2[1:3,1:20]
Stands_s3[1:3,1:20]
Stands_s4[1:3,1:20]

Taxon_std_s1 = Stands_s1
Taxon_std_s2 = Stands_s2
Taxon_std_s3 = Stands_s3
Taxon_std_s4 = Stands_s4

Annot_Taxon =  annot_taxon
#head(Annot_Taxon)
dim(Annot_Taxon) # 6965 12
# ---------------------------------------------------------------------------;
#
#                         Merging for DNA pathway analysis
# 
# ---------------------------------------------------------------------------;
lname = load(file="/udd/hpxil/3.source/hpfs_microbiome_biobakery4/ReadIn_stratified_pathway.RData")
print(lname)
#[1] "annot_stra_pathway" "spath_raw_total"    "spath_raw_s1"       "spath_raw_s2"       "spath_raw_s3"      
#[6] "spath_raw_s4"       "spath_std_total"    "spath_std_s1"       "spath_std_s2"       "spath_std_s3"      
#[11] "spath_std_s4"   

lname_un = load(file="/udd/hpxil/3.source/hpfs_microbiome_biobakery4/ReadIn_unstratified_pathway.RData")
print(lname_un)
#[1] "annot_unstra_pathway" "unpath_raw_total"     "unpath_raw_s1"        "unpath_raw_s2"        "unpath_raw_s3"       
#[6] "unpath_raw_s4"        "unpath_std_total"     "unpath_std_s1"        "unpath_std_s2"        "unpath_std_s3"       
#[11] "unpath_std_s4"       
# for raw data

dna_pth_raw_s1 = spath_raw_s1
dna_pth_raw_s2 = spath_raw_s2
dna_pth_raw_s3 = spath_raw_s3
dna_pth_raw_s4 = spath_raw_s4

dna_pth_std_s1 = spath_std_s1
dna_pth_std_s2 = spath_std_s2
dna_pth_std_s3 = spath_std_s3
dna_pth_std_s4 = spath_std_s4

dna_upth_raw_s1 = unpath_raw_s1
dna_upth_raw_s2 = unpath_raw_s2
dna_upth_raw_s3 = unpath_raw_s3
dna_upth_raw_s4 = unpath_raw_s4

dna_upth_std_s1 = unpath_std_s1
dna_upth_std_s2 = unpath_std_s2
dna_upth_std_s3 = unpath_std_s3
dna_upth_std_s4 = unpath_std_s4



# delete the "_s1,2,3 or 4" at the back at every bug so that we can rbind data from 4 times into long format

getsubchr = function(a,destart,deend) {	substr(a,destart+1,nchar(a)-deend) }
names(dna_pth_raw_s1)[-1] = getsubchr(names(dna_pth_raw_s1)[-1],0,3)
names(dna_pth_raw_s2)[-1] = getsubchr(names(dna_pth_raw_s2)[-1],0,3)
names(dna_pth_raw_s3)[-1] = getsubchr(names(dna_pth_raw_s3)[-1],0,3)
names(dna_pth_raw_s4)[-1] = getsubchr(names(dna_pth_raw_s4)[-1],0,3)

names(dna_pth_std_s1)[-1] = getsubchr(names(dna_pth_std_s1)[-1],0,3)
names(dna_pth_std_s2)[-1] = getsubchr(names(dna_pth_std_s2)[-1],0,3)
names(dna_pth_std_s3)[-1] = getsubchr(names(dna_pth_std_s3)[-1],0,3)
names(dna_pth_std_s4)[-1] = getsubchr(names(dna_pth_std_s4)[-1],0,3)

names(dna_upth_raw_s1)[-1] = getsubchr(names(dna_upth_raw_s1)[-1],0,3)
names(dna_upth_raw_s2)[-1] = getsubchr(names(dna_upth_raw_s2)[-1],0,3)
names(dna_upth_raw_s3)[-1] = getsubchr(names(dna_upth_raw_s3)[-1],0,3)
names(dna_upth_raw_s4)[-1] = getsubchr(names(dna_upth_raw_s4)[-1],0,3)

names(dna_upth_std_s1)[-1] = getsubchr(names(dna_upth_std_s1)[-1],0,3)
names(dna_upth_std_s2)[-1] = getsubchr(names(dna_upth_std_s2)[-1],0,3)
names(dna_upth_std_s3)[-1] = getsubchr(names(dna_upth_std_s3)[-1],0,3)
names(dna_upth_std_s4)[-1] = getsubchr(names(dna_upth_std_s4)[-1],0,3)




#-----------------------------------------------------------------------------
# for DNA enzyme -- use both raw and std,  
lname = load(file="/udd/hpxil/3.source/hpfs_microbiome_biobakery4/ReadIn_stratified_enzyme.RData")
lname_un = load(file="/udd/hpxil/3.source/hpfs_microbiome_biobakery4/ReadIn_unstratified_enzyme.RData")
# 查看RData文件中包含的对象名称
print(lname)
#[1] "annot_stra_enzyme" "senzy_raw_total"   "senzy_raw_s1"      "senzy_raw_s2"      "senzy_raw_s3"      "senzy_raw_s4"     
#[7] "senzy_std_total"   "senzy_std_s1"      "senzy_std_s2"      "senzy_std_s3"      "senzy_std_s4"     
print(lname_un)
#[1] "annot_unstra_enzyme" "unenzy_raw_total"    "unenzy_raw_s1"       "unenzy_raw_s2"       "unenzy_raw_s3"      
#[6] "unenzy_raw_s4"       "unenzy_std_total"    "unenzy_std_s1"       "unenzy_std_s2"       "unenzy_std_s3"      
#[11] "unenzy_std_s4"
#-----------------------------------------------------------------------------

dim(senzy_raw_s1) # 154819
dim(senzy_raw_s2) 
dim(senzy_raw_s3) 
dim(senzy_raw_s4) 

dim(senzy_std_s1) 
dim(senzy_std_s2) 
dim(senzy_std_s3) 
dim(senzy_std_s4) 

senzy_raw_s1[1:3,1:20]
senzy_raw_s2[1:3,1:20]
senzy_raw_s3[1:3,1:20]
senzy_raw_s4[1:3,1:20]

senzy_std_s1[1:3,1:20]
senzy_std_s2[1:3,1:20]
senzy_std_s3[1:3,1:20]
senzy_std_s4[1:3,1:20]

#head(Annot_DNA)
dim(annot_stra_enzyme) # 154818 2

dna_enz_raw_s1 = senzy_raw_s1
dna_enz_raw_s2 = senzy_raw_s2
dna_enz_raw_s3 = senzy_raw_s3
dna_enz_raw_s4 = senzy_raw_s4

dna_enz_std_s1 = senzy_std_s1
dna_enz_std_s2 = senzy_std_s2
dna_enz_std_s3 = senzy_std_s3
dna_enz_std_s4 = senzy_std_s4

identical(names(senzy_std_s1),names(senzy_std_s2))


dim(unenzy_raw_s1) # 2424
dim(unenzy_raw_s2) 
dim(unenzy_raw_s3) 
dim(unenzy_raw_s4) 

dim(unenzy_std_s1) 
dim(unenzy_std_s2) 
dim(unenzy_std_s3) 
dim(unenzy_std_s4) 

unenzy_raw_s1[1:3,1:20]
unenzy_raw_s2[1:3,1:20]
unenzy_raw_s3[1:3,1:20]
unenzy_raw_s4[1:3,1:20]

unenzy_std_s1[1:3,1:20]
unenzy_std_s2[1:3,1:20]
unenzy_std_s3[1:3,1:20]
unenzy_std_s4[1:3,1:20]

#head(Annot_DNA)
dim(annot_unstra_enzyme) # 2434

dna_unenz_raw_s1 = unenzy_raw_s1
dna_unenz_raw_s2 = unenzy_raw_s2
dna_unenz_raw_s3 = unenzy_raw_s3
dna_unenz_raw_s4 = unenzy_raw_s4

dna_unenz_std_s1 = unenzy_std_s1
dna_unenz_std_s2 = unenzy_std_s2
dna_unenz_std_s3 = unenzy_std_s3
dna_unenz_std_s4 = unenzy_std_s4

identical(names(unenzy_std_s1),names(unenzy_std_s2))

# ---------------------------------------------------------------------------;
#
#                          Merging for analysis
#
# ---------------------------------------------------------------------------;
#                         Basic phenotype files
# ---------------------------------------------------------------------------;

# load(file="CheckPointData_foranalysis.RData")


#-------------
# Time point 1

dim(Basic1) # 415 69
dim(FFQ1) # 727 52
dim(DDR1) # 210 550

S1 = merge(data.frame(id = Basic1$id, SampleID = paste("S",Basic1$id,sep=''), IDTime = paste("S",Basic1$id,".1",sep=''), AnalySet=1), Basic1,by.x="id",by.y="id")
S1 = merge(S1,FFQ1,by.x="id",by.y="id",all=T)
S1 = merge(S1,DDR1,by.x="id",by.y="id",all.y=T)

#head(S1)
dim(S1) # 210 672
#-------------
# Time point 2

dim(Basic2) #  415  69
dim(FFQ2) # 727  52
dim(DDR2) #205 550

S2 = merge(data.frame(id = Basic2$id, SampleID = paste("S",Basic2$id,sep=''), IDTime = paste("S",Basic2$id,".2",sep=''), AnalySet=2), Basic2,by.x="id",by.y="id")
S2 = merge(S2,FFQ2,by.x="id",by.y="id",all=T)
S2 = merge(S2,DDR2,by.x="id",by.y="id",all.y=T)

#head(S2)
dim(S2) # 205 672

#-------------
# Time point 3

dim(Basic3) # 415 69
dim(FFQ3) # 727 52
dim(DDR3) # 256 550

S3 = merge(data.frame(id = Basic3$id, SampleID = paste("S",Basic3$id,sep=''), IDTime = paste("S",Basic3$id,".3",sep=''), AnalySet=3), Basic3,by.x="id",by.y="id")
S3 = merge(S3,FFQ3,by.x="id",by.y="id",all=T)
S3 = merge(S3,DDR3,by.x="id",by.y="id",all.y=T)

#head(S3)
dim(S3) # 256 672

#-------------
# Time point 4

dim(Basic4) # 415  69
dim(FFQ4) # 727  52
dim(DDR4) # 254 550

S4 = merge(data.frame(id = Basic4$id, SampleID = paste("S",Basic4$id,sep=''), IDTime = paste("S",Basic4$id,".4",sep=''), AnalySet=4), Basic4,by.x="id",by.y="id")
S4 = merge(S4,FFQ4,by.x="id",by.y="id",all=T)
S4 = merge(S4,DDR4,by.x="id",by.y="id",all.y=T)

#head(S4)
dim(S4) # 254 672



# -------------------------------------------
# merging data

AllPhenos = rbind(rbind(rbind(S1,S2),S3),S4)
#head(AllPhenos)
dim(AllPhenos) # 925 672


# -------------------------------------------
# check on missing

# check vars to be adjusted

AdjVars = c("bmi_bld","time_fec","fast_bld","totMETs_paq","alco_fs_ffq","calor_fs_ffq","lt_fiber",
            "lt_whgrn","lt_ahei","lt_ahei_na","lt_ahei_nowgr","lt_readmt","fish_ffq", "nut_ffq",
            "probio_2m_fec","antibio_12m_fec","colsc_2m_fec","ims_12m_fec","acid_2m_fec","chm_12m_fec","adiarr_2m_fec",
            "stooltype_fec.1","stooltype_fec.2","stooltype_fec.3","stooltype_fec.4","stooltype_fec.5","stooltype_fec.6")

# check missing

check = data.frame(var=AdjVars,missing=NA)
for(i in 1:length(AdjVars)) {
  check$missing[i]=length(AllPhenos[is.na(AllPhenos[,AdjVars[i]]),AdjVars[i]])
}
print(check[which(check$missing!=0),])


#              var missing
#2       time_fec      23
#5    alco_fs_ffq      30
#6   calor_fs_ffq      30
#13      fish_ffq       4
#14       nut_ffq       4
#15 probio_2m_fec      27
#17  colsc_2m_fec      25
#19   acid_2m_fec      26
#21 adiarr_2m_fec      23

AllPhenos$time_fec[is.na(AllPhenos$time_fec)] = 1  # assuming morning, which is the most choosen time

AllPhenos$calor_fs_ffq[is.na(AllPhenos$calor_fs_ffq)] = median(AllPhenos$calor_fs_ffq,na.rm=TRUE) 
AllPhenos$alco_fs_ffq[is.na(AllPhenos$alco_fs_ffq)] = median(AllPhenos$alco_fs_ffq,na.rm=TRUE) 
AllPhenos$fish_ffq[is.na(AllPhenos$fish_ffq)] = median(AllPhenos$fish_ffq,na.rm=TRUE)
AllPhenos$nut_ffq[is.na(AllPhenos$nut_ffq)] = median(AllPhenos$nut_ffq,na.rm=TRUE) 
AllPhenos$nut07_avg5d[is.na(AllPhenos$nut07_avg5d)] = median(AllPhenos$nut07_avg5d,na.rm=TRUE) 

AllPhenos$fish_avg5d[is.na(AllPhenos$fish_avg5d)] = median(AllPhenos$fish_avg5d,na.rm=TRUE) 

AllPhenos$probio_2m_fec[is.na(AllPhenos$probio_2m_fec)] = 0
AllPhenos$colsc_2m_fec[is.na(AllPhenos$colsc_2m_fec)] = 0
AllPhenos$acid_2m_fec[is.na(AllPhenos$acid_2m_fec)] = 0
AllPhenos$adiarr_2m_fec[is.na(AllPhenos$adiarr_2m_fec)] = 0


# -------------------------------------------
# Make major dietary variables to deciles for certain analysis
# _dm means median of deciles

varlist = c("calor_fs_ffq","wholegrain_ffq","fruit_ffq","vege_ffq","nut_ffq","legume_ffq","ssb_ffq","dairy_ffq","egg_ffq","fish_ffq",
            "b6_fs_ffq","fol98_fs_ffq","b12_fs_ffq","b6_fo_ffq","fol98_fo_ffq","b12_fo_ffq","whgrn_fs_ffq","whsub_fs_ffq","crude_fs_ffq",
            "carbo_fs_ffq","aofib_fs_ffq","engl_fs_ffq","prot_fs_ffq","aprot_fs_ffq","vprot_fs_ffq","tfat_fs_ffq","afat_fs_ffq","vfat_fs_ffq",
            "meat_ffq","livers_ffq","meatliv_ffq","redmt_ffq","redliv_ffq","redmtliv_ffq","ptrymt_ffq","ptryliv_ffq","ptrymtliv_ffq",
            "betaine_fs_ffq","frcho_fs_ffq","gpcho_fs_ffq","pcho_fs_ffq","ptdcho_fs_ffq","sphingo_fs_ffq","choline_fs_ffq","betchol_fs_ffq",
            "betaine_fo_ffq","frcho_fo_ffq","ptdcho_fo_ffq","choline_fo_ffq","betchol_fo_ffq", # ffq
            x, # 7 ddr
            "bmi_bld","totMETs_paq","lt_fiber","lt_whgrn","lt_ahei","lt_ahei_na","lt_ahei_nowgr","lt_readmt","lt_nut","lt_veg", "lt_fru","lt_trypt") # others


for(i in 1:length(varlist)) {
  
  vari = varlist[i]
  
  AllPhenos$new = makequantiles(vari,AllPhenos,10,"median")
  names(AllPhenos)[which(names(AllPhenos)=="new")] = paste(vari,"_dm",sep='')
  
  
}

#names(AllPhenos)
dim(AllPhenos) # 925854
# save all phenos
# save(AllPhenos,file="AllPhenos.RData")
AllPhenos$pbut07_avg1d
summary(AllPhenos[,c("pnut07_avg1d","pbut07_avg1d","onut07_avg1d","wnut07_avg1d","pwnut07_avg1d","nut07_avg1d")])
summary(AllPhenos[,c("pnut07_avg2d","pbut07_avg2d","onut07_avg2d","wnut07_avg2d","pwnut07_avg2d","nut07_avg2d")])
summary(AllPhenos[,c("pnut07_avg3d","pbut07_avg3d","onut07_avg3d","wnut07_avg3d","pwnut07_avg3d","nut07_avg3d")])
summary(AllPhenos[,c("pnut07_avg4d","pbut07_avg4d","onut07_avg4d","wnut07_avg4d","pwnut07_avg4d","nut07_avg4d")])
summary(AllPhenos[,c("pnut07_avg5d","pbut07_avg5d","onut07_avg5d","wnut07_avg5d","pwnut07_avg5d","nut07_avg5d")])
dim(AllPhenos)
# ---------------------------------------------------------------------------;
#
#                         Merging for Taxon analysis
#
# ---------------------------------------------------------------------------;

# bind all taxon data

Taxon_std_s1$AnalySet = 1
Taxon_std_s2$AnalySet = 2
Taxon_std_s3$AnalySet = 3
Taxon_std_s4$AnalySet = 4

identical(names(Taxon_std_s1),names(Taxon_std_s2))
identical(names(Taxon_std_s1),names(Taxon_std_s3))
identical(names(Taxon_std_s1),names(Taxon_std_s4))

Taxon_std = rbind(rbind(rbind(Taxon_std_s1,Taxon_std_s2),Taxon_std_s3),Taxon_std_s4)
#head(Taxon_std)[,1:5]
#head(Taxon_std)[,(dim(Taxon_std)[2]-5):dim(Taxon_std)[2]]
dim(Taxon_std) # 925 6970


# merge with phenos

Taxon_use = merge(AllPhenos,Taxon_std,by=c("id","AnalySet"))
#names(Taxon_use)
#head(Taxon_use)[,1:5]
#head(Taxon_use)[,(dim(Taxon_use)[2]-5):dim(Taxon_use)[2]]
dim(Taxon_use) #  925 7822

names(Taxon_use)
# ---------------------------------------------------------------------------;
#
#                         Merging for DNA pathway analysis
# 
# ---------------------------------------------------------------------------;

dna_pth_raw_s1$AnalySet = 1
dna_pth_raw_s2$AnalySet = 2
dna_pth_raw_s3$AnalySet = 3
dna_pth_raw_s4$AnalySet = 4

identical(names(dna_pth_raw_s1),names(dna_pth_raw_s2))
identical(names(dna_pth_raw_s1),names(dna_pth_raw_s3))
identical(names(dna_pth_raw_s1),names(dna_pth_raw_s4))
dna_pth_raw = rbind(rbind(rbind(dna_pth_raw_s1,dna_pth_raw_s2),dna_pth_raw_s3),dna_pth_raw_s4)
#head(dna_pth_raw)[,1:5]
#head(dna_pth_raw)[,(dim(dna_pth_raw)[2]-5):dim(dna_pth_raw)[2]]
dim(dna_pth_raw) # 925 20334
table(AllPhenos$AnalySet)
table(dna_pth_raw$AnalySet)
dna_pth_raw_use = merge(AllPhenos,dna_pth_raw,by=c("id","AnalySet"))
#names(dna_pth_raw_use)
#head(dna_pth_raw_use)[,1:5]
#head(dna_pth_raw_use)[,(dim(dna_pth_raw_use)[2]-5):dim(dna_pth_raw_use)[2]]
dim(dna_pth_raw_use) #  925 21186

# for std data

dna_pth_std_s1$AnalySet = 1
dna_pth_std_s2$AnalySet = 2
dna_pth_std_s3$AnalySet = 3
dna_pth_std_s4$AnalySet = 4

identical(names(dna_pth_std_s1),names(dna_pth_std_s2))
identical(names(dna_pth_std_s1),names(dna_pth_std_s3))
identical(names(dna_pth_std_s1),names(dna_pth_std_s4))

dna_pth_std = rbind(rbind(rbind(dna_pth_std_s1,dna_pth_std_s2),dna_pth_std_s3),dna_pth_std_s4)
#head(dna_pth_std)[,1:5]
#head(dna_pth_std)[,(dim(dna_pth_std)[2]-5):dim(dna_pth_std)[2]]
dim(dna_pth_std) # 925 9147

dna_pth_std_use = merge(AllPhenos,dna_pth_std,by=c("id","AnalySet"))
#names(dna_pth_std_use)
#head(dna_pth_std_use)[,1:5]
#head(dna_pth_std_use)[,(dim(dna_pth_std_use)[2]-5):dim(dna_pth_std_use)[2]]
dim(dna_pth_std_use) #  925 21186

# ----------------------Unstratified Pathway-----------------------------------------------------;
dna_upth_raw_s1$AnalySet = 1
dna_upth_raw_s2$AnalySet = 2
dna_upth_raw_s3$AnalySet = 3
dna_upth_raw_s4$AnalySet = 4

identical(names(dna_upth_raw_s1),names(dna_upth_raw_s2))
identical(names(dna_upth_raw_s1),names(dna_upth_raw_s3))
identical(names(dna_upth_raw_s1),names(dna_upth_raw_s4))
dna_upth_raw = rbind(rbind(rbind(dna_upth_raw_s1,dna_upth_raw_s2),dna_upth_raw_s3),dna_upth_raw_s4)
#head(dna_upth_raw)[,1:5]
#head(dna_upth_raw)[,(dim(dna_upth_raw)[2]-5):dim(dna_upth_raw)[2]]
dim(dna_upth_raw) # 925 509
table(AllPhenos$AnalySet)
table(dna_upth_raw$AnalySet)
dna_upth_raw_use = merge(AllPhenos,dna_upth_raw,by=c("id","AnalySet"))
#names(dna_upth_raw_use)
#head(dna_upth_raw_use)[,1:5]
#head(dna_upth_raw_use)[,(dim(dna_upth_raw_use)[2]-5):dim(dna_upth_raw_use)[2]]
dim(dna_upth_raw_use) #  925 1361

# for std data

dna_upth_std_s1$AnalySet = 1
dna_upth_std_s2$AnalySet = 2
dna_upth_std_s3$AnalySet = 3
dna_upth_std_s4$AnalySet = 4

identical(names(dna_upth_std_s1),names(dna_upth_std_s2))
identical(names(dna_upth_std_s1),names(dna_upth_std_s3))
identical(names(dna_upth_std_s1),names(dna_upth_std_s4))

dna_upth_std = rbind(rbind(rbind(dna_upth_std_s1,dna_upth_std_s2),dna_upth_std_s3),dna_upth_std_s4)
#head(dna_upth_std)[,1:5]
#head(dna_upth_std)[,(dim(dna_upth_std)[2]-5):dim(dna_upth_std)[2]]
dim(dna_upth_std) # 925 9147

dna_upth_std_use = merge(AllPhenos,dna_upth_std,by=c("id","AnalySet"))
#names(dna_upth_std_use)
#head(dna_upth_std_use)[,1:5]
#head(dna_upth_std_use)[,(dim(dna_upth_std_use)[2]-5):dim(dna_upth_std_use)[2]]
dim(dna_upth_std_use) #  925 1857


# ---------------------------------------------------------------------------;
#
#                       Merging for DNA Enzyme analysis
#sftp://127.99.136.12/udd/hpxil/3.source/hpfs_microbiome_biobakery4/ReadIn_stratified_enzyme.RData
# ---------------------------------------------------------------------------;

# for raw data
dna_enz_raw_s1$AnalySet = 1
dna_enz_raw_s2$AnalySet = 2
dna_enz_raw_s3$AnalySet = 3
dna_enz_raw_s4$AnalySet = 4

identical(names(dna_enz_raw_s1),names(dna_enz_raw_s2))
identical(names(dna_enz_raw_s1),names(dna_enz_raw_s3))
identical(names(dna_enz_raw_s1),names(dna_enz_raw_s4))

dna_enz_raw = rbind(rbind(rbind(dna_enz_raw_s1,dna_enz_raw_s2),dna_enz_raw_s3),dna_enz_raw_s4)
#head(dna_enz_raw)[,1:5]
#head(dna_enz_raw)[,(dim(dna_enz_raw)[2]-5):dim(dna_enz_raw)[2]]
dim(dna_enz_raw) #  925 154820
dim(AllPhenos)
dna_enz_raw_use = merge(AllPhenos,dna_enz_raw,by=c("id","AnalySet"))
#names(dna_enz_raw_use)
#head(dna_enz_raw_use)[,1:5]
#head(dna_enz_raw_use)[,(dim(dna_enz_raw_use)[2]-5):dim(dna_enz_raw_use)[2]]
dim(dna_enz_raw_use) #  new=925 155672


# for std data
#[1] "annot_stra_enzyme" "senzy_raw_total"   "senzy_raw_s1"      "senzy_raw_s2"      "senzy_raw_s3"      "senzy_raw_s4"     
#[7] "senzy_std_total"   "senzy_std_s1"      "senzy_std_s2"      "senzy_std_s3"      "senzy_std_s4"     

dna_enz_std_s1$AnalySet = 1
dna_enz_std_s2$AnalySet = 2
dna_enz_std_s3$AnalySet = 3
dna_enz_std_s4$AnalySet = 4

identical(names(dna_enz_std_s1),names(dna_enz_std_s2))
identical(names(dna_enz_std_s1),names(dna_enz_std_s3))
identical(names(dna_enz_std_s1),names(dna_enz_std_s4))

dna_enz_std = rbind(rbind(rbind(dna_enz_std_s1,dna_enz_std_s2),dna_enz_std_s3),dna_enz_std_s4)
#head(dna_enz_std)[,1:5]
#head(dna_enz_std)[,(dim(dna_enz_std)[2]-5):dim(dna_enz_std)[2]]
dim(dna_enz_std) # 925 94861  new= 925 154820

dna_enz_std_use = merge(AllPhenos,dna_enz_std,by=c("id","AnalySet"))
#names(dna_enz_std_use)
#head(dna_enz_std_use)[,1:5]
#head(dna_enz_std_use)[,(dim(dna_enz_std_use)[2]-5):dim(dna_enz_std_use)[2]]
dim(dna_enz_std_use) #  925 155672

# save DNA and tax data in case there's problem in next steps
# save(Taxon_use, dna_pth_raw_use, dna_pth_std_use, dna_enz_raw_use, dna_enz_std_use, Annot_Taxon, Annot_DNApath, Annot_DNA,
# 	file="Taxon_DNA_temp.RData")


# for raw data
dna_unenz_raw_s1$AnalySet = 1
dna_unenz_raw_s2$AnalySet = 2
dna_unenz_raw_s3$AnalySet = 3
dna_unenz_raw_s4$AnalySet = 4

identical(names(dna_unenz_raw_s1),names(dna_unenz_raw_s2))
identical(names(dna_unenz_raw_s1),names(dna_unenz_raw_s3))
identical(names(dna_unenz_raw_s1),names(dna_unenz_raw_s4))

dna_unenz_raw = rbind(rbind(rbind(dna_unenz_raw_s1,dna_unenz_raw_s2),dna_unenz_raw_s3),dna_unenz_raw_s4)
#head(dna_enz_raw)[,1:5]
#head(dna_enz_raw)[,(dim(dna_enz_raw)[2]-5):dim(dna_enz_raw)[2]]
dim(dna_unenz_raw) # 925 2425
dim(AllPhenos)
dna_unenz_raw_use = merge(AllPhenos,dna_unenz_raw,by=c("id","AnalySet"))
#names(dna_enz_raw_use)
#head(dna_enz_raw_use)[,1:5]
#head(dna_enz_raw_use)[,(dim(dna_enz_raw_use)[2]-5):dim(dna_enz_raw_use)[2]]
dim(dna_unenz_raw_use) # 925 3277


# for std data
#[1] "annot_stra_enzyme" "unenzy_raw_total"   "unenzy_raw_s1"      "unenzy_raw_s2"      "unenzy_raw_s3"      "unenzy_raw_s4"     
#[7] "unenzy_std_total"   "unenzy_std_s1"      "unenzy_std_s2"      "unenzy_std_s3"      "unenzy_std_s4"     

dna_unenz_std_s1$AnalySet = 1
dna_unenz_std_s2$AnalySet = 2
dna_unenz_std_s3$AnalySet = 3
dna_unenz_std_s4$AnalySet = 4

identical(names(dna_unenz_std_s1),names(dna_unenz_std_s2))
identical(names(dna_unenz_std_s1),names(dna_unenz_std_s3))
identical(names(dna_unenz_std_s1),names(dna_unenz_std_s4))

dna_unenz_std = rbind(rbind(rbind(dna_unenz_std_s1,dna_unenz_std_s2),dna_unenz_std_s3),dna_unenz_std_s4)
#head(dna_enz_std)[,1:5]
#head(dna_enz_std)[,(dim(dna_enz_std)[2]-5):dim(dna_enz_std)[2]]
dim(dna_unenz_std) # 925 2425

dna_unenz_std_use = merge(AllPhenos,dna_unenz_std,by=c("id","AnalySet"))
#names(dna_enz_std_use)
#head(dna_enz_std_use)[,1:5]
#head(dna_enz_std_use)[,(dim(dna_enz_std_use)[2]-5):dim(dna_enz_std_use)[2]]
dim(dna_unenz_std_use) #  925 3277

# taxon
#head(Taxon_use)[1:2,1:20]
dim(Taxon_use) # 925 8318
names(Taxon_use)
names(annot_taxon)

identical(names(Taxon_use)[1351:8315],annot_taxon$taxon_new)  
#head(Annot_Taxon)
dim(annot_taxon) # 6965   12

# DNA
identical(rownames(dna_pth_raw_use),rownames(dna_pth_std_use))
identical(names(dna_pth_raw_use),names(dna_pth_std_use))

#head(dna_pth_raw_use)[1:3,1:50]
dim(dna_pth_raw_use) # 925 21682

#head(dna_pth_std_use)[1:3,1:50]
dim(dna_pth_std_use) #  925 21682
names(dna_pth_std_use)
identical(names(dna_pth_std_use)[855:21186],annot_stra_pathway$microName)
head(Annot_DNApath)
dim(Annot_DNApath) # 9145 9


identical(rownames(dna_enz_raw_use),rownames(dna_enz_std_use))
identical(names(dna_enz_raw_use),names(dna_enz_std_use))

#head(dna_enz_raw_use)[1:3,1:50]
dim(dna_enz_raw_use) # 925  155672

#head(dna_enz_std_use)[1:3,1:50]
dim(dna_enz_std_use) # 925  155672
names(dna_enz_std_use)
identical(names(dna_enz_std_use)[855:155672],annot_stra_enzyme$enzyme_id)

# ---------------------------------------------------------------------------;
#
#                   Check and save all the files
#
# ---------------------------------------------------------------------------;


# taxon
#head(Taxon_use)[1:2,1:20]
dim(Taxon_use) # 925 8674

#head(Annot_Taxon)
dim(Annot_Taxon) # 6965 12

#-----------------------------------------------------------------
# save the data for final use
#-----------------------------------------------------------------
save( Taxon_use,Annot_Taxon,
      dna_pth_raw_use,dna_pth_std_use,dna_upth_raw_use,dna_upth_std_use,annot_stra_pathway,annot_unstra_pathway,
      dna_enz_raw_use,dna_enz_std_use,dna_unenz_raw_use,dna_unenz_std_use,annot_stra_enzyme,annot_unstra_enzyme,
      lvs_hpfs_pm,lvs_hpfs_fea,
      file="/udd/n2gji/micro/nut_final/data/Final_Data_2nd.RData")

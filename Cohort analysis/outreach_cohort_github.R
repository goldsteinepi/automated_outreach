##################
# Automated outreach cohort analysis
# Citation: Kilburn K, Goldstein ND, Biddle T, Glaycomb C, Bilinski A, Kirchner HL, Chabris CF, Goren A. Evaluating Automated Post-Discharge Outreach Among Surgical Patients Using Retrospective EHR Data and Complementary Causal Methods. Manuscript in preparation.
# 2/21/26 -- Neal Goldstein
##################


### FUNCTIONS ###

library(readxl) #read xlsx file
library(cobalt) #iptw diagnostics
library(psych) #describe, describeBy
library(gmodels) #CrossTable
library(forestplot) #forest plots for regression results
library(sandwich) #robust SEs
library(marginaleffects) #average marginal effects aka risk difference
library(lme4) #mixed effects modeling, note: need to install from source: install.packages("lme4", type = "source")

#robust standard errors for RR estimation of a binary outcome using Poisson (Zou, AJE, 2004; PMID: 15033648)
#modified Poisson with cluster-robust standard errors equivalent to random intercept (Yelland, Salter, Ryan, AJE, 2011; PMID: 21841157)
robust_ci = function(model, group=NA) {
  
  #determine whether clustered errors are needed
  if (length(group)>1) {
    
    #robust covariance matrix with clustering
    cov_matrix = vcovCL(model, cluster=group, type="HC0")
    
  } else {
    
    #robust covariance matrix w/out clustering
    cov_matrix = vcovHC(model, type="HC0")
    
  }
  
  #robust SE
  robust_se = sqrt(diag(cov_matrix))
  
  #coefficients
  betas = coef(model)
  
  #95% CIs on log scale
  lower = betas - 1.96 * robust_se
  upper = betas + 1.96 * robust_se
  
  #return CI
  cbind(lower, upper)
}

#define null-coalescing operator (used by lme4 but not installed due to lme4/Matrix conflict with present R version)
`%||%` <- function(a, b) if (!is.null(a)) a else b


### READ DATA ###

setwd("/Automated Outreach/Surgical/Data Pulls/")

#Epic/Clarity discharge data
discharges = read.csv("All discharges/1117GPEKK1205250504260735-Post Surgery Outreach-v1.08.csv", as.is=T, na.strings="")
str(discharges)

#recast variable data types as needed
#discharges$PAT_ENC_CSN_ID = as.numeric(gsub(",","",discharges$PAT_ENC_CSN_ID))
discharges$ADMSN_DT = as.POSIXct(discharges$ADMSN_DT, format="%Y-%m-%d", tz="UTC")
discharges$DISCH_DT = as.POSIXct(discharges$DISCH_DT, format="%Y-%m-%d", tz="UTC")
discharges$Death_Dt = as.POSIXct(discharges$Death_Dt, format="%Y-%m-%d", tz="UTC")
discharges$BMI = as.numeric(discharges$BMI) #empty strings will cause a warning
discharges$Frst_Unplanned_Readmit_Dt  = as.POSIXct(discharges$Frst_Unplanned_Readmit_Dt, format="%Y-%m-%d", tz="UTC")
discharges$Frst_30_Day_ED_Adm_Dt  = as.POSIXct(discharges$Frst_30_Day_ED_Adm_Dt, format="%Y-%m-%d", tz="UTC")

#PDC dashboard export (for automated outreach identification); open CSV as saved from Tablaeu and save as XLSX to fix import error
pdc_dashboard <- read_excel("Dashboards/PDC Dashboard full 20-Feb-2026.xlsx")
str(pdc_dashboard)

#OR dashboard export (for surgery identification); open CSV as saved from Tablaeu and save as XLSX to fix import error
#or_dashboard <- read_excel("Dashboards/OR Dashboard 2025 Surgeries 5-Jan-2026.xlsx")
#str(or_dashboard)

#Geisinger facilities
Geisinger_EDs = c("17756","18510","18711","17740","17866","17822","17044","17822","18702","17815") #https://locations.geisinger.org/?k=&z=&lt=Emergency+department
Geisinger_Hospitals = c("18510","17740","17822","17822","17044","18702","17866","18711","17815","17756") #https://locations.geisinger.org/?k=&z=&lt=Hospital

#Geisinger facility regions
#https://www.geisinger.org/-/media/OneGeisinger/pdfs/ghs/about-geisinger/pdfs/Geisinger-System-Map---August-2020.pdf

#Surgical department mappings data
surgical_depts <- read_excel("../R/Department matchings_2.19.26.xlsx")
str(surgical_depts)


# ### COMPARISON to PDC DASHBOARD ###
# 
# #outreach comparison
# sum(discharges$Pat_Outreach_Flg, na.rm=T) #12719
# nrow(pdc_dashboard) #76726
# 
# #proportion PDC dashboard CSNs found in discharge CSNs (Pat_Outreach_Flg)
# sum(pdc_dashboard$`Pat Enc Csn Id` %in% discharges$PAT_ENC_CSN_ID[discharges$Pat_Outreach_Flg==1]) / nrow(pdc_dashboard) * 100
# #7%
# 
# #sample a few CSNs not found for Nate
# csns_not_found = pdc_dashboard$`Pat Enc Csn Id`[which(!(pdc_dashboard$`Pat Enc Csn Id` %in% discharges$PAT_ENC_CSN_ID[discharges$Pat_Outreach_Flg==1]))]
# sample(csns_not_found, 10)
# 
# #proportion PDC dashboard patient IDs found in discharge patient IDs (Pat_Outreach_Flg)
# sum(pdc_dashboard$`Pat Id` %in% discharges$PAT_ID[discharges$Pat_Outreach_Flg==1]) / nrow(pdc_dashboard) * 100
# #8%
# 
# #sample a few Patient IDs not found for Nate
# patids_not_found = pdc_dashboard$`Pat Id`[which(!(pdc_dashboard$`Pat Id` %in% discharges$PAT_ID[discharges$Pat_Outreach_Flg==1]))]
# sample(patids_not_found, 10)
# 
# #proportion PDC dashboard CSNs found in discharge CSNs (all)
# sum(pdc_dashboard$`Pat Enc Csn Id` %in% discharges$PAT_ENC_CSN_ID) / nrow(pdc_dashboard) * 100
# #56%
# 
# #sample a few CSNs not found for Nate
# csns_not_found = pdc_dashboard$`Pat Enc Csn Id`[which(!(pdc_dashboard$`Pat Enc Csn Id` %in% discharges$PAT_ENC_CSN_ID))]
# sample(csns_not_found, 10)
# 
# #proportion PDC dashboard patient IDs found in discharge patient IDs (all)
# sum(pdc_dashboard$`Pat Id` %in% discharges$PAT_ID) / nrow(pdc_dashboard) * 100
# #79%
# 
# #sample a few Patient IDs not found for Nate
# patids_not_found = pdc_dashboard$`Pat Id`[which(!(pdc_dashboard$`Pat Id` %in% discharges$PAT_ID))]
# sample(patids_not_found, 10)
# 
# 
# ### COMPARISON to OR DASHBOARD ###
# 
# #surgeries comparison (need to limit to same dates since OR dashboard is just 2025)
# sum(discharges$OR_Dash_Case[discharges$ADMSN_DT>=as.POSIXct("2025-01-01") & discharges$ADMSN_DT<=as.POSIXct("2025-12-31")], na.rm=T) #105220
# nrow(or_dashboard) #80657
# 
# #proportion OR dashboard CSNs found in discharge CSNs (OR_Dash_Case)
# sum(or_dashboard$`Pat Enc Csn Id` %in% discharges$PAT_ENC_CSN_ID[discharges$OR_Dash_Case==1]) / nrow(or_dashboard) * 100
# #<1%
# 
# #sample a few CSNs not found for Nate
# csns_not_found = or_dashboard$`Pat Enc Csn Id`[which(!(or_dashboard$`Pat Enc Csn Id` %in% discharges$PAT_ENC_CSN_ID[discharges$OR_Dash_Case==1]))]
# sample(csns_not_found, 10)
# 
# #proportion OR dashboard patient IDs are found in discharge patient IDs (OR_Dash_Case)
# sum(or_dashboard$`Pat Id` %in% discharges$PAT_ID[discharges$OR_Dash_Case==1]) / nrow(or_dashboard) * 100
# #79%
# 
# #sample a few Patient IDs not found for Nate
# patids_not_found = or_dashboard$`Pat Id`[which(!(or_dashboard$`Pat Id` %in% discharges$PAT_ID[discharges$OR_Dash_Case==1]))]
# sample(patids_not_found, 10)
# 
# #proportion OR dashboard CSNs found in discharge CSNs (all)
# sum(or_dashboard$`Pat Enc Csn Id` %in% discharges$PAT_ENC_CSN_ID) / nrow(or_dashboard) * 100
# #2%
# 
# #sample a few CSNs not found for Nate
# csns_not_found = or_dashboard$`Pat Enc Csn Id`[which(!(or_dashboard$`Pat Enc Csn Id` %in% discharges$PAT_ENC_CSN_ID))]
# sample(csns_not_found, 10)
# 
# #proportion OR dashboard patient IDs are found in discharge patient IDs (all)
# sum(or_dashboard$`Pat Id` %in% discharges$PAT_ID) / nrow(or_dashboard) * 100
# #84%
# 
# #sample a few Patient IDs not found for Nate
# patids_not_found = or_dashboard$`Pat Id`[which(!(or_dashboard$`Pat Id` %in% discharges$PAT_ID))]
# sample(patids_not_found, 10)


### LINK DASHBOARDS to DISCHARGES ###

#create unique keys in the discharge and PDC datasets by patient ID and discharge date
discharges$patient_id_discharge_date = paste(discharges$PAT_ID, discharges$DISCH_DT)
pdc_dashboard$patient_id_discharge_date = paste(pdc_dashboard$`Pat Id`, pdc_dashboard$`Disch Dt`)

#identify the matches based on type of outreach
discharges$pdc_dashboard_flag_toc = ifelse(discharges$patient_id_discharge_date %in% pdc_dashboard$patient_id_discharge_date[grep("PCP TOC", pdc_dashboard$`Outreach Name`)], 1, 0)
discharges$pdc_dashboard_flag_surg = ifelse(discharges$patient_id_discharge_date %in% pdc_dashboard$patient_id_discharge_date[grep("Surg", pdc_dashboard$`Outreach Name`)], 1, 0)
discharges$pdc_dashboard_flag_copd = ifelse(discharges$patient_id_discharge_date %in% pdc_dashboard$patient_id_discharge_date[grep("COPD", pdc_dashboard$`Outreach Name`)], 1, 0)
discharges$pdc_dashboard_flag_diabetes = ifelse(discharges$patient_id_discharge_date %in% pdc_dashboard$patient_id_discharge_date[grep("Diabetes", pdc_dashboard$`Outreach Name`)], 1, 0)
pdc_dashboard$discharges_flag = ifelse(pdc_dashboard$patient_id_discharge_date %in% discharges$patient_id_discharge_date, 1, 0)

table(discharges$Pat_Outreach_Flg)
table(discharges$pdc_dashboard_flag_toc)
table(discharges$pdc_dashboard_flag_surg)
table(discharges$pdc_dashboard_flag_copd)
table(discharges$pdc_dashboard_flag_diabetes)

# #for the OR dashboard, check if the surgery date falls within the inpatient hospitalization
# merged_data = merge(discharges, or_dashboard, by.x="PAT_ID", by.y="Pat Id")
# 
# #keep only valid date range matches
# valid_matches = merged_data[(merged_data$`Surg Date` >= merged_data$ADMSN_DT) & (merged_data$`Surg Date` <= merged_data$DISCH_DT), ]
# 
# #explore distribution of surgeries w/ and w/out PDC outreach
# surgeries_pdc = as.matrix(table(valid_matches$Service[valid_matches$pdc_dashboard_flag==1]))
# surgeries_all = as.matrix(table(valid_matches$Service[valid_matches$Service %in% row.names(surgeries_pdc)]))
# surgeries_prop = (surgeries_pdc/surgeries_all)
# 
# #only retain surgeries w/ >= 10% outreach (a proxy for higher risk of readmission)
# #quantile(surgeries_prop, probs=c(0.25,0.75))
# #valid_matches = valid_matches[valid_matches$Service %in% row.names(surgeries_prop)[which(surgeries_prop>0.1)], ]
# valid_matches = valid_matches[valid_matches$Service %in% row.names(surgeries_prop)[which(surgeries_prop>0.03 & surgeries_prop<0.24)], ]
# 
# #create unique keys for matching
# key1_valid = unique(paste(valid_matches$PAT_ID, valid_matches$ADMSN_DT, valid_matches$DISCH_DT))
# key1_all = paste(discharges$PAT_ID, discharges$ADMSN_DT, discharges$DISCH_DT)
# key2_valid = unique(paste(valid_matches$PAT_ID, valid_matches$`Surg Date`))
# key2_all = paste(or_dashboard$`Pat Id`, or_dashboard$`Surg Date`)
# 
# #identify the matches
# discharges$or_dashboard_flag = as.integer(key1_all %in% key1_valid)
# or_dashboard$discharges_flag = as.integer(key2_all %in% key2_valid)
# 
# prop.table(table(or_dashboard$discharges_flag))
# #91% match rate from or_dashboard to discharges
# 
# #cleanup
# rm(merged_data, valid_matches, key1_all, key1_valid, key2_all, key2_valid, surgeries_pdc, surgeries_all, surgeries_prop)


### CREATE COHORT using AUTOMATED OUTREACH ELIGIBILITY CRITERIA ###

#adult >=18
discharges = discharges[discharges$AGE_AT_DISCHARGE >= 18, ]

#language preference English, Spanish, Unknown (missing)
table(discharges$LANG_NAME, useNA="always")
discharges = discharges[discharges$LANG_NAME == "English" | discharges$LANG_NAME == "Spanish", ]

#discharge date of January 14, 2025 or later (systemwide go-live)
discharges = discharges[discharges$DISCH_DT >= as.POSIXct("2025-01-14", tz="UTC"), ]

#evaluation end date of December 31, 2025
discharges = discharges[discharges$DISCH_DT <= as.POSIXct("2025-12-31", tz="UTC"), ]

#surgical department classification
discharges$surgery_dept = surgical_depts$`Outreach name`[match(discharges$PRIM_WORKDAY_DEPT_NAME, surgical_depts$PRIM_WORKDAY_DEPT_NAME)]

#identify potential immortal time bias: automated outreach occurs next business day (or two) -- cannot have mortality or outcome before exposure
#hist((as.numeric(discharges$Frst_Unplanned_Readmit_Dt - discharges$DISCH_DT)/(24*60*60)), breaks="fd")
discharges$outreach_effect_1weekday = as.POSIXct(ifelse(weekdays(discharges$DISCH_DT)=="Friday", (discharges$DISCH_DT + 3*86400), 
                                                        ifelse(weekdays(discharges$DISCH_DT)=="Saturday", (discharges$DISCH_DT + 2*86400), 
                                                               (discharges$DISCH_DT + 1*86400))), tz="UTC") #use seconds to shift day
discharges$outreach_effect_2weekday = as.POSIXct(ifelse(weekdays(discharges$DISCH_DT)=="Thursday", (discharges$DISCH_DT + 4*86400),
                                                        ifelse(weekdays(discharges$DISCH_DT)=="Friday", (discharges$DISCH_DT + 4*86400), 
                                                               ifelse(weekdays(discharges$DISCH_DT)=="Saturday", (discharges$DISCH_DT + 3*86400), 
                                                                      (discharges$DISCH_DT + 2*86400)))), tz="UTC") #use seconds to shift day

#bariatric surgeries have exposure in week 2 (add additional 7 day offset)
discharges$outreach_effect_1weekday[!is.na(discharges$surgery_dept) & discharges$surgery_dept=="Bariatric Surgery"] = discharges$outreach_effect_1weekday[!is.na(discharges$surgery_dept) & discharges$surgery_dept=="Bariatric Surgery"] + 7*86400 #use seconds to shift day
discharges$outreach_effect_2weekday[!is.na(discharges$surgery_dept) & discharges$surgery_dept=="Bariatric Surgery"] = discharges$outreach_effect_2weekday[!is.na(discharges$surgery_dept) & discharges$surgery_dept=="Bariatric Surgery"] + 7*86400 #use seconds to shift day

#create immortal time bias indicators
discharges$death_before_outreach_1day = ifelse((discharges$Death_Dt <= discharges$outreach_effect_1weekday) & !(is.na(discharges$Death_Dt)), 1, 0) 
discharges$readmit_before_outreach_1day = ifelse((discharges$Frst_Unplanned_Readmit_Dt <= discharges$outreach_effect_1weekday) & !(is.na(discharges$Frst_Unplanned_Readmit_Dt)), 1, 0) #use seconds to shift day
discharges$ed_before_outreach_1day = ifelse((discharges$Frst_30_Day_ED_Adm_Dt <= discharges$outreach_effect_1weekday) & !(is.na(discharges$Frst_30_Day_ED_Adm_Dt)), 1, 0) #use seconds to shift day
discharges$death_before_outreach_2day = ifelse((discharges$Death_Dt <= discharges$outreach_effect_2weekday) & !(is.na(discharges$Death_Dt)), 1, 0) 
discharges$readmit_before_outreach_2day = ifelse((discharges$Frst_Unplanned_Readmit_Dt <= discharges$outreach_effect_2weekday) & !(is.na(discharges$Frst_Unplanned_Readmit_Dt)), 1, 0) #use seconds to shift day
discharges$ed_before_outreach_2day = ifelse((discharges$Frst_30_Day_ED_Adm_Dt <= discharges$outreach_effect_2weekday) & !(is.na(discharges$Frst_30_Day_ED_Adm_Dt)), 1, 0) #use seconds to shift day

table(discharges$death_before_outreach_1day, useNA="always")
table(discharges$readmit_before_outreach_1day, useNA="always")
table(discharges$ed_before_outreach_1day, useNA="always")
table(discharges$death_before_outreach_2day, useNA="always")
table(discharges$readmit_before_outreach_2day, useNA="always")
table(discharges$ed_before_outreach_2day, useNA="always")

#other criteria (language, encounter type, discharge disposition, and ADT patient class) were already applied during Nate's data export
#Language Preference: English, Unable to Obtain; updated to add Spanish
#Encounter Type: Hospital Encounter w/ Admission and Discharge Times populated
#Discharge Disposition: exclude Expired; updated to only Home - Self Care, Home with Services, Home with Services - Planned Readmission, Home - Planned Readmission, AMA
#ADT Pat Class: Inpatient, Outpatient, Emergency, Observation, SORU (23 hr Surgical Overnight)

#cohort definition one: limit to surgeries with targeted outreach
discharges$cohort_surgeries_targeted = ifelse(discharges$PRIM_WORKDAY_DEPT_NAME %in% surgical_depts$PRIM_WORKDAY_DEPT_NAME[surgical_depts$Treatment==1], 1, 0)
table(discharges$cohort_surgeries_targeted, useNA="always")

#cohort definition two: limit to non-targeted surgeries
#discharges$cohort_surgeries_notargeted = ifelse(discharges$PRIM_WORKDAY_DEPT_NAME %in% surgical_depts$PRIM_WORKDAY_DEPT_NAME, 1, 0)
#table(discharges$cohort_surgeries_notargeted, useNA="always")

#cohort definition three: medical discharges (PCP TOC, Diabetes PDC, or COPD PDC)
#discharges$cohort_medical = ifelse(cohort_analysis$ADT_PAT_CLASS_NAME=="Inpatient" & )

#create cohort with targeted surgeries and immortal time bias restrictions
#note: being overly restrictive by dropping both outcomes, but keeps it simpler to present results as a single cohort
cohort_analysis = discharges[discharges$cohort_surgeries_targeted==1 & discharges$death_before_outreach_1day==0 & discharges$readmit_before_outreach_1day==0 & discharges$ed_before_outreach_1day==0, ]
#cohort_analysis = discharges[discharges$cohort_surgeries_targeted==1 & discharges$death_before_outreach_2day==0 & discharges$readmit_before_outreach_2day==0 & discharges$ed_before_outreach_2day==0, ]

#define intervention condition and comparison condition (exclude other outreach efforts from comparison group)
cohort_analysis$Pat_Outreach_Flg = cohort_analysis$pdc_dashboard_flag_surg
cohort_analysis = cohort_analysis[cohort_analysis$pdc_dashboard_flag_copd==0 & cohort_analysis$pdc_dashboard_flag_diabetes==0 & cohort_analysis$pdc_dashboard_flag_toc==0, ]


### RECODES and DATA OPERATIONALIZATION ###

#sex
table(cohort_analysis$SEX_NAME, useNA="always")
cohort_analysis$SEX_NAME = ifelse(cohort_analysis$SEX_NAME=="Unknown", NA, cohort_analysis$SEX_NAME)

#payors, GHP vs not, Medicare vs not, Medicaid vs not
table(cohort_analysis$Prim_Payor_name, useNA="always")

cohort_analysis$Prim_Payor_GHP = 0
cohort_analysis$Prim_Payor_Medicare = 0
cohort_analysis$Prim_Payor_Medicaid = 0

cohort_analysis$Prim_Payor_GHP[grep("GHP", cohort_analysis$Prim_Payor_name)] = 1
cohort_analysis$Prim_Payor_Medicare[grep("MEDICARE", cohort_analysis$Prim_Payor_name, ignore.case=T)] = 1
cohort_analysis$Prim_Payor_Medicaid[grep("MEDICAID", cohort_analysis$Prim_Payor_name, ignore.case=T)] = 1
cohort_analysis$Prim_Payor = ifelse(cohort_analysis$Prim_Payor_GHP==1, 0, ifelse(cohort_analysis$Prim_Payor_Medicare==1, 1, ifelse(cohort_analysis$Prim_Payor_Medicaid==1, 2, 3)))
table(cohort_analysis$Prim_Payor)

#remove invalid ages
#hist(cohort_analysis$AGE_AT_DISCHARGE, breaks="fd")
#cohort_analysis$AGE_AT_DISCHARGE[cohort_analysis$AGE_AT_DISCHARGE>=110] = NA

#transform LOS
hist(cohort_analysis$LOS, breaks="fd")
cohort_analysis$LOS_log = log(cohort_analysis$LOS)
cohort_analysis$LOS_log[which(is.infinite(cohort_analysis$LOS_log))] = NA

#missing assumed zero for CCI calculation
table(cohort_analysis$CCI, useNA="always")
cohort_analysis$CCI[is.na(cohort_analysis$CCI)] = 0

#binary primary outcomes
cohort_analysis$Unplanned_Readmit_30_Day = as.integer(!is.na(cohort_analysis$Unplanned_Readmit_30_Day_Cnt))
cohort_analysis$ED_30_Day = as.integer(!is.na(cohort_analysis$ED_30_Day_Cnt))

#count outcomes
table(cohort_analysis$Unplanned_Readmit_30_Day_Cnt, useNA="always")
cohort_analysis$Unplanned_Readmit_30_Day_Cnt = ifelse(is.na(cohort_analysis$Unplanned_Readmit_30_Day_Cnt), 0, cohort_analysis$Unplanned_Readmit_30_Day_Cnt)
table(cohort_analysis$ED_30_Day_Cnt, useNA="always")
cohort_analysis$ED_30_Day_Cnt = ifelse(is.na(cohort_analysis$ED_30_Day_Cnt), 0, cohort_analysis$ED_30_Day_Cnt)

#race, 0=white, 1=black, 2=other
table(cohort_analysis$RACE_NAME, useNA="always")
cohort_analysis$RACE_NAME = ifelse(cohort_analysis$RACE_NAME=="Patient Declined to Provide" | cohort_analysis$RACE_NAME=="Unable to Obtain", NA, cohort_analysis$RACE_NAME)
cohort_analysis$RACE_NAME_3cat = ifelse(cohort_analysis$RACE_NAME=="White", 0, ifelse(cohort_analysis$RACE_NAME=="Black or African American", 1, 2))

#ethnicity, 0=non-hispanic, 1=hispanic
table(cohort_analysis$ETHNICITY_NAME, useNA="always")
cohort_analysis$ETHNICITY_NAME = ifelse(cohort_analysis$ETHNICITY_NAME=="Patient Declined to Provide" | cohort_analysis$ETHNICITY_NAME=="Unable to Obtain", NA, cohort_analysis$ETHNICITY_NAME)
cohort_analysis$ETHNICITY_NAME_2cat = ifelse(cohort_analysis$ETHNICITY_NAME=="Hispanic or Latino" | cohort_analysis$ETHNICITY_NAME=="Hispanic,Asian/Pacific", 1, 0)

#Geisinger facility in ZIP
cohort_analysis$Geisinger_Hospital = as.integer(cohort_analysis$PRIM_ZIP %in% Geisinger_Hospitals)
cohort_analysis$Geisinger_ED = as.integer(cohort_analysis$PRIM_ZIP %in% Geisinger_EDs)
table(cohort_analysis$Geisinger_Hospital, useNA="always")
table(cohort_analysis$Geisinger_ED, useNA="always")

#patient type
table(cohort_analysis$ADT_PAT_CLASS_NAME, useNA="always")
cohort_analysis$Pt_type = ifelse(cohort_analysis$ADT_PAT_CLASS_NAME=="Inpatient", 0, ifelse(cohort_analysis$ADT_PAT_CLASS_NAME=="Outpatient", 1, 2))
cohort_analysis$Pt_type_inpatient = ifelse(cohort_analysis$ADT_PAT_CLASS_NAME=="Inpatient", 1, 0)
cohort_analysis$Pt_type_outpatient = ifelse(cohort_analysis$ADT_PAT_CLASS_NAME=="Outpatient", 1, 0)
cohort_analysis$Pt_type_emergency = ifelse(cohort_analysis$ADT_PAT_CLASS_NAME=="Emergency", 1, 0)
table(cohort_analysis$Pt_type_inpatient, useNA="always")
table(cohort_analysis$Pt_type_outpatient, useNA="always")
table(cohort_analysis$Pt_type_emergency, useNA="always")

#month of visit (trend)
cohort_analysis$DISCH_DT_month = as.integer(format(cohort_analysis$DISCH_DT, "%m"))

#automated outreach
#table(cohort_analysis$Pat_Outreach_Flg, useNA="always")
#cohort_analysis$Pat_Outreach_Flg = ifelse(is.na(cohort_analysis$Pat_Outreach_Flg), 0, 1)

#readmission risk scores, due to large amount of missingness we are not using for now
#table(is.na(cohort_analysis$Admn_Readm_Risk_Score)) #last available RRS prior to the admission timestamp
#table(is.na(cohort_analysis$Pre_Dsch_Readm_Risk_Score)) #last available RRS that was calculated during the duration of the patient admission prior to discharge
#table(is.na(cohort_analysis$Post_Dsch_Readm_Risk_Score)) #first available RRS following the patient discharge timestamp
#cohort_analysis$Admn_Readm_Risk_Score_missing = ifelse(is.na(cohort_analysis$Admn_Readm_Risk_Score), 1, 0)
#CrossTable(cohort_analysis$Pat_Outreach_Flg, cohort_analysis$Admn_Readm_Risk_Score_missing, prop.r=F, prop.t=F, prop.chisq=F, chisq=T)
#CrossTable(cohort_analysis$Unplanned_Readmit_30_Day, cohort_analysis$Admn_Readm_Risk_Score_missing, prop.r=F, prop.t=F, prop.chisq=F, chisq=T)



# ### DESCRIPTIVES for GEISINGER SYMPOSIUM POSTER ###
# 
# #time-series by year-month
# cohort_analysis$year_month = format(cohort_analysis$DISCH_DT, "%Y-%m")
# 
# #aggregate outreach
# outreach_ts = aggregate(
#   list(count = cohort_analysis$or_dashboard_flag[cohort_analysis$pdc_dashboard_flag==1]),
#   by = list(year_month = cohort_analysis$year_month[cohort_analysis$pdc_dashboard_flag==1],
#             surgery = cohort_analysis$or_dashboard_flag[cohort_analysis$pdc_dashboard_flag==1]),
#   FUN = length
# )
# 
# #outreach plot
# plot(x=1:length(outreach_ts$year_month[outreach_ts$surgery==0]), y=outreach_ts$count[outreach_ts$surgery==0], type="l", xlab="", ylab="Outreaches (n)", xaxt="n", font.lab=2, font=2, ylim=c(1000,2000))
# lines(x=1:length(outreach_ts$year_month[outreach_ts$surgery==1]), y=outreach_ts$count[outreach_ts$surgery==1], lty=2)
# axis(1, at=1:length(outreach_ts$year_month[outreach_ts$surgery==0]), labels=paste(outreach_ts$year_month[outreach_ts$surgery==0]), las=2, font=2)
# legend("topleft", lty=c(1,2), legend=c("Medical","Surgical"), text.font=2, cex=0.8)
# 
# #descriptives for outreach
# sum(cohort_analysis$pdc_dashboard_flag, na.rm=T)
# prop.table(table(cohort_analysis$or_dashboard_flag[cohort_analysis$pdc_dashboard_flag==1]))
# 
# #aggregate readmissions
# readmits_ts = aggregate(
#   list(count = cohort_analysis$or_dashboard_flag[cohort_analysis$Unplanned_Readmit_30_Day_Cnt>0]),
#   by = list(year_month = cohort_analysis$year_month[cohort_analysis$Unplanned_Readmit_30_Day_Cnt>0],
#             outreach = cohort_analysis$pdc_dashboard_flag[cohort_analysis$Unplanned_Readmit_30_Day_Cnt>0]),
#   FUN = length
# )
# 
# #aggregate ED visits
# edvisits_ts = aggregate(
#   list(count = cohort_analysis$or_dashboard_flag[cohort_analysis$ED_30_Day_Cnt>0]),
#   by = list(year_month = cohort_analysis$year_month[cohort_analysis$ED_30_Day_Cnt>0],
#             outreach = cohort_analysis$pdc_dashboard_flag[cohort_analysis$ED_30_Day_Cnt>0]),
#   FUN = length
# )
# 
# #denominator for proportion calculations
# monthly_totals = aggregate(
#   x = list(count = cohort_analysis$pdc_dashboard_flag),
#   by = list(
#     year_month = cohort_analysis$year_month,
#     outreach = cohort_analysis$pdc_dashboard_flag
#   ),
#   FUN = length
# )
# 
# #remove invalid date
# #monthly_totals = monthly_totals[-15, ]
# 
# #readmission plot by outreach group
# 
# #calculate proportion
# readmits_ts$discharges = monthly_totals$count
# readmits_ts$count_prop = readmits_ts$count / readmits_ts$discharges
# 
# #create matrix of counts by outreach (rows) and year_month (columns)
# bar_matrix = with(readmits_ts, tapply(count_prop, list(outreach, year_month), sum))
# 
# barplot(bar_matrix, beside = T, col = c("steelblue", "tomato"), ylab = "Proportion", las = 2, ylim=c(0,0.25))
# legend("topleft", legend = c("Automated", "Standard of care"), fill = c("steelblue", "tomato"), bty = "n", horiz=T)
# 
# #descriptives for readmissions
# sum(readmits_ts$count, na.rm=T)
# 
# #ED visits plot by outreach group
# 
# #calculate proportion
# edvisits_ts$discharges = monthly_totals$count
# edvisits_ts$count_prop = edvisits_ts$count / edvisits_ts$discharges
# 
# #create matrix of counts by outreach (rows) and year_month (columns)
# bar_matrix = with(edvisits_ts, tapply(count_prop, list(outreach, year_month), sum))
# 
# barplot(bar_matrix, beside = T, col = c("steelblue", "tomato"), ylab = "Proportion", las = 2, ylim=c(0,0.25))
# legend("topleft", legend = c("Automated", "Standard of care"), fill = c("steelblue", "tomato"), bty = "n", horiz=T)
# 
# #descriptives for ED visitgs
# sum(edvisits_ts$count, na.rm=T)


### DESCRIPTIVES ###

#number of patients
length(unique(cohort_analysis$PAT_ID))

CrossTable(cohort_analysis$Pat_Outreach_Flg)
describe(cohort_analysis$AGE_AT_DISCHARGE)
CrossTable(cohort_analysis$SEX_NAME)
CrossTable(cohort_analysis$LANG_NAME)
CrossTable(cohort_analysis$RACE_NAME_3cat)
CrossTable(cohort_analysis$ETHNICITY_NAME_2cat)
table(cohort_analysis$Prim_Payor_name, useNA="always") #note large number of missing
CrossTable(cohort_analysis$Prim_Payor)
describe(cohort_analysis$LOS)
describe(cohort_analysis$CCI)
CrossTable(cohort_analysis$Geisinger_ED)
CrossTable(cohort_analysis$Geisinger_Hospital) #since Hospital and ED are equal, just need to model one
CrossTable(cohort_analysis$Pt_type)
describe(cohort_analysis$DISCH_DT_month)
CrossTable(cohort_analysis$surgery_dept)
#describe(cohort_analysis$Admn_Readm_Risk_Score)
#describe(cohort_analysis$Pre_Dsch_Readm_Risk_Score)
#describe(cohort_analysis$Post_Dsch_Readm_Risk_Score)
#describe(cohort_analysis$Unplanned_Readmit_30_Day_Cnt)
CrossTable(cohort_analysis$Unplanned_Readmit_30_Day)
#describe(cohort_analysis$ED_30_Day_Cnt)
CrossTable(cohort_analysis$ED_30_Day)

#trends over time
plot(aggregate(Unplanned_Readmit_30_Day_Cnt ~ DISCH_DT_month, FUN=sum, data=cohort_analysis), type="l", ylim=c(0,500), xlab="Month (in 2025)", ylab="Unplanned readmits (30 day)", lwd=2)
plot(aggregate(ED_30_Day_Cnt ~ DISCH_DT_month, FUN=sum, data=cohort_analysis), type="l", ylim=c(0,500), xlab="Month (in 2025)", ylab="ED visits (30 day)", lwd=2)
plot(aggregate(Pat_Outreach_Flg ~ DISCH_DT_month, FUN=sum, data=cohort_analysis), type="l", xlab="Month (in 2025)", ylab="Automated outreach", lwd=2, ylim=c(0,2000))

#averaged outcomes
sum(cohort_analysis$Unplanned_Readmit_30_Day, na.rm=T) / nrow(cohort_analysis) * 100
sum(cohort_analysis$ED_30_Day, na.rm=T) / nrow(cohort_analysis) * 100
sum(cohort_analysis$Unplanned_Readmit_30_Day[cohort_analysis$Pat_Outreach_Flg==1], na.rm=T) / nrow(cohort_analysis[cohort_analysis$Pat_Outreach_Flg==1,]) * 100
sum(cohort_analysis$Unplanned_Readmit_30_Day[cohort_analysis$Pat_Outreach_Flg==0], na.rm=T) / nrow(cohort_analysis[cohort_analysis$Pat_Outreach_Flg==0,]) * 100
sum(cohort_analysis$ED_30_Day[cohort_analysis$Pat_Outreach_Flg==1], na.rm=T) / nrow(cohort_analysis[cohort_analysis$Pat_Outreach_Flg==1,]) * 100
sum(cohort_analysis$ED_30_Day[cohort_analysis$Pat_Outreach_Flg==0], na.rm=T) / nrow(cohort_analysis[cohort_analysis$Pat_Outreach_Flg==0,]) * 100

#by outreach
describeBy(cohort_analysis$AGE_AT_DISCHARGE, cohort_analysis$Pat_Outreach_Flg)
CrossTable(cohort_analysis$SEX_NAME, cohort_analysis$Pat_Outreach_Flg, prop.r=F, prop.t=F, prop.chisq=F, chisq=T)
CrossTable(cohort_analysis$LANG_NAME, cohort_analysis$Pat_Outreach_Flg, prop.r=F, prop.t=F, prop.chisq=F, chisq=T)
CrossTable(cohort_analysis$RACE_NAME_3cat, cohort_analysis$Pat_Outreach_Flg, prop.r=F, prop.t=F, prop.chisq=F, chisq=T)
CrossTable(cohort_analysis$ETHNICITY_NAME_2cat, cohort_analysis$Pat_Outreach_Flg, prop.r=F, prop.t=F, prop.chisq=F, chisq=T)
CrossTable(cohort_analysis$Prim_Payor, cohort_analysis$Pat_Outreach_Flg, prop.r=F, prop.t=F, prop.chisq=F, chisq=T)
describeBy(cohort_analysis$LOS, cohort_analysis$Pat_Outreach_Flg)
describeBy(cohort_analysis$CCI, cohort_analysis$Pat_Outreach_Flg)   #wilcox.test(cohort_analysis$CCI ~ cohort_analysis$Pat_Outreach_Flg)
CrossTable(cohort_analysis$Geisinger_Hospital, cohort_analysis$Pat_Outreach_Flg, prop.r=F, prop.t=F, prop.chisq=F, chisq=T)
CrossTable(cohort_analysis$Pt_type, cohort_analysis$Pat_Outreach_Flg, prop.r=F, prop.t=F, prop.chisq=F, chisq=T)
describeBy(cohort_analysis$DISCH_DT_month, cohort_analysis$Pat_Outreach_Flg)
CrossTable(cohort_analysis$surgery_dept, cohort_analysis$Pat_Outreach_Flg, prop.r=F, prop.t=F, prop.chisq=F, chisq=T)
#describeBy(cohort_analysis$Admn_Readm_Risk_Score, cohort_analysis$Pat_Outreach_Flg)
#describeBy(cohort_analysis$Unplanned_Readmit_30_Day_Cnt, cohort_analysis$Pat_Outreach_Flg)
CrossTable(cohort_analysis$Unplanned_Readmit_30_Day, cohort_analysis$Pat_Outreach_Flg, prop.r=F, prop.t=F, prop.chisq=F, chisq=T)
#describeBy(cohort_analysis$ED_30_Day_Cnt, cohort_analysis$Pat_Outreach_Flg)
CrossTable(cohort_analysis$ED_30_Day, cohort_analysis$Pat_Outreach_Flg, prop.r=F, prop.t=F, prop.chisq=F, chisq=T)

#by readmission, binary
CrossTable(cohort_analysis$Pat_Outreach_Flg, cohort_analysis$Unplanned_Readmit_30_Day, prop.r=F, prop.t=F, prop.chisq=F, chisq=T)
describeBy(cohort_analysis$AGE_AT_DISCHARGE, cohort_analysis$Unplanned_Readmit_30_Day)
CrossTable(cohort_analysis$SEX_NAME, cohort_analysis$Unplanned_Readmit_30_Day, prop.r=F, prop.t=F, prop.chisq=F, chisq=T)
CrossTable(cohort_analysis$LANG_NAME, cohort_analysis$Unplanned_Readmit_30_Day, prop.r=F, prop.t=F, prop.chisq=F, chisq=T)
CrossTable(cohort_analysis$RACE_NAME_3cat, cohort_analysis$Unplanned_Readmit_30_Day, prop.r=F, prop.t=F, prop.chisq=F, chisq=T)
CrossTable(cohort_analysis$ETHNICITY_NAME_2cat, cohort_analysis$Unplanned_Readmit_30_Day, prop.r=F, prop.t=F, prop.chisq=F, chisq=T)
CrossTable(cohort_analysis$Prim_Payor, cohort_analysis$Unplanned_Readmit_30_Day, prop.r=F, prop.t=F, prop.chisq=F, chisq=T)
describeBy(cohort_analysis$LOS, cohort_analysis$Unplanned_Readmit_30_Day)
describeBy(cohort_analysis$CCI, cohort_analysis$Unplanned_Readmit_30_Day)   #wilcox.test(cohort_analysis$CCI ~ cohort_analysis$Unplanned_Readmit_30_Day)
CrossTable(cohort_analysis$Geisinger_Hospital, cohort_analysis$Unplanned_Readmit_30_Day, prop.r=F, prop.t=F, prop.chisq=F, chisq=T)
CrossTable(cohort_analysis$Pt_type, cohort_analysis$Unplanned_Readmit_30_Day, prop.r=F, prop.t=F, prop.chisq=F, chisq=T)
describeBy(cohort_analysis$DISCH_DT_month, cohort_analysis$Unplanned_Readmit_30_Day)
CrossTable(cohort_analysis$surgery_dept, cohort_analysis$Unplanned_Readmit_30_Day, prop.r=F, prop.t=F, prop.chisq=F, chisq=T)
#describeBy(cohort_analysis$Admn_Readm_Risk_Score, cohort_analysis$Unplanned_Readmit_30_Day)

#by ED visit, binary
CrossTable(cohort_analysis$Pat_Outreach_Flg, cohort_analysis$ED_30_Day, prop.r=F, prop.t=F, prop.chisq=F, chisq=T)
describeBy(cohort_analysis$AGE_AT_DISCHARGE, cohort_analysis$ED_30_Day)
CrossTable(cohort_analysis$SEX_NAME, cohort_analysis$ED_30_Day, prop.r=F, prop.t=F, prop.chisq=F, chisq=T)
CrossTable(cohort_analysis$LANG_NAME, cohort_analysis$ED_30_Day, prop.r=F, prop.t=F, prop.chisq=F, chisq=T)
CrossTable(cohort_analysis$RACE_NAME_3cat, cohort_analysis$ED_30_Day, prop.r=F, prop.t=F, prop.chisq=F, chisq=T)
CrossTable(cohort_analysis$ETHNICITY_NAME_2cat, cohort_analysis$ED_30_Day, prop.r=F, prop.t=F, prop.chisq=F, chisq=T)
CrossTable(cohort_analysis$Prim_Payor, cohort_analysis$ED_30_Day, prop.r=F, prop.t=F, prop.chisq=F, chisq=T)
describeBy(cohort_analysis$LOS, cohort_analysis$ED_30_Day)
describeBy(cohort_analysis$CCI, cohort_analysis$ED_30_Day)   #wilcox.test(cohort_analysis$CCI ~ cohort_analysis$ED_30_Day)
CrossTable(cohort_analysis$Geisinger_Hospital, cohort_analysis$ED_30_Day, prop.r=F, prop.t=F, prop.chisq=F, chisq=T)
CrossTable(cohort_analysis$Pt_type, cohort_analysis$ED_30_Day, prop.r=F, prop.t=F, prop.chisq=F, chisq=T)
describeBy(cohort_analysis$DISCH_DT_month, cohort_analysis$ED_30_Day)
CrossTable(cohort_analysis$surgery_dept, cohort_analysis$ED_30_Day, prop.r=F, prop.t=F, prop.chisq=F, chisq=T)
#describeBy(cohort_analysis$Admn_Readm_Risk_Score, cohort_analysis$ED_30_Day)


### IPTW MODEL for AUTOMATED OUTREACH ###

#predict automated outreach for confounding effect of who received it
model_ps_iptw = glm(Pat_Outreach_Flg ~ CCI + as.factor(RACE_NAME_3cat) + ETHNICITY_NAME_2cat + LANG_NAME + as.factor(Prim_Payor) + AGE_AT_DISCHARGE + LOS_log + DISCH_DT_month, family=binomial(link="logit"), data=cohort_analysis)

summary(model_ps_iptw)
exp(coef(model_ps_iptw))
exp(confint(model_ps_iptw))

#predicted probabilities of automated outreach
probs = predict(model_ps_iptw, type="response")
cohort_analysis$PS_IPTW = NA
cohort_analysis[names(probs), "PS_IPTW"] = probs
hist(cohort_analysis$PS_IPTW, breaks="fd")
boxplot(cohort_analysis$PS_IPTW)

#inverse probability of treatment weights
cohort_analysis$IPTW = ifelse(cohort_analysis$Pat_Outreach_Flg==1, 1/cohort_analysis$PS_IPTW, 1/(1-cohort_analysis$PS_IPTW))
hist(cohort_analysis$IPTW, breaks="fd")
boxplot(cohort_analysis$IPTW)

#stabilize the weights
cohort_analysis$IPTW_stabil = ifelse(cohort_analysis$Pat_Outreach_Flg==1, mean(cohort_analysis$Pat_Outreach_Flg, na.rm=T)/cohort_analysis$PS_IPTW, (1-mean(cohort_analysis$Pat_Outreach_Flg, na.rm=T))/(1-cohort_analysis$PS_IPTW))
hist(cohort_analysis$IPTW_stabil, breaks="fd")
boxplot(cohort_analysis$IPTW_stabil)

#trim the extreme values
IPTW_quantile = quantile(cohort_analysis$IPTW_stabil, c(0.01, 0.99), na.rm=T)
cohort_analysis$IPTW_stabil_trimmed = ifelse(cohort_analysis$IPTW_stabil < IPTW_quantile[1], IPTW_quantile[1], cohort_analysis$IPTW_stabil)
cohort_analysis$IPTW_stabil_trimmed = ifelse(cohort_analysis$IPTW_stabil > IPTW_quantile[2], IPTW_quantile[2], cohort_analysis$IPTW_stabil)
hist(cohort_analysis$IPTW_stabil_trimmed, breaks="fd")
boxplot(cohort_analysis$IPTW_stabil_trimmed)

#IPTW diagnostics
cohort_analysis_complete = cohort_analysis[!is.na(cohort_analysis$PS_IPTW), ]
plot_labels = c("Race, White","Race, Black","Race, Other","Hispanic","Language, Spanish","Payor, Geisinger","Payor,Medicare","Payor, Medicaid","Payor, Other","Age","Length of stay","Discharge month")
names(plot_labels) = c("as.factor(RACE_NAME_3cat)_0","as.factor(RACE_NAME_3cat)_1","as.factor(RACE_NAME_3cat)_2","ETHNICITY_NAME_2cat","LANG_NAME_Spanish","as.factor(Prim_Payor)_0","as.factor(Prim_Payor)_1","as.factor(Prim_Payor)_2","as.factor(Prim_Payor)_3","AGE_AT_DISCHARGE","LOS_log","DISCH_DT_month")
bal.tab(Pat_Outreach_Flg ~ as.factor(RACE_NAME_3cat) + ETHNICITY_NAME_2cat + LANG_NAME + as.factor(Prim_Payor) + AGE_AT_DISCHARGE + LOS_log + DISCH_DT_month, data=cohort_analysis_complete, weights=cohort_analysis_complete$IPTW_stabil, method="weighting")
love.plot(Pat_Outreach_Flg ~ as.factor(RACE_NAME_3cat) + ETHNICITY_NAME_2cat + LANG_NAME + as.factor(Prim_Payor) + AGE_AT_DISCHARGE + LOS_log + DISCH_DT_month, data=cohort_analysis_complete, weights=cohort_analysis_complete$IPTW_stabil, threshold=0.1, var.names=plot_labels)

rm(probs, model_ps_iptw, cohort_analysis_complete)

#if not using IPSW
IPSW_stabil = 1


# ### IPSW MODEL for AUTOMATED OUTREACH
# 
# #predict automated outreach for selection bias inclusion/exclusion of surgeries
# model_ps_ipsw = glm(Pat_Outreach_Flg ~ as.factor(surgery_dept), family=binomial(link="logit"), data=cohort_analysis)
# 
# summary(model_ps_ipsw)
# exp(coef(model_ps_ipsw))
# exp(confint(model_ps_ipsw))
# 
# #predicted probabilities of automated outreach
# probs = predict(model_ps_ipsw, type="response")
# cohort_analysis$PS_IPSW = NA
# cohort_analysis[names(probs), "PS_IPSW"] = probs
# hist(cohort_analysis$PS_IPSW, breaks="fd")
# boxplot(cohort_analysis$PS_IPSW)
# 
# #inverse probability of treatment weights
# cohort_analysis$IPSW = ifelse(cohort_analysis$Pat_Outreach_Flg==1, 1/cohort_analysis$PS_IPSW, 1/(1-cohort_analysis$PS_IPSW))
# hist(cohort_analysis$IPSW, breaks="fd")
# boxplot(cohort_analysis$IPSW)
# 
# #stabilize the weights
# cohort_analysis$IPSW_stabil = ifelse(cohort_analysis$Pat_Outreach_Flg==1, mean(cohort_analysis$Pat_Outreach_Flg, na.rm=T)/cohort_analysis$PS_IPSW, (1-mean(cohort_analysis$Pat_Outreach_Flg, na.rm=T))/(1-cohort_analysis$PS_IPSW))
# hist(cohort_analysis$IPSW_stabil, breaks="fd")
# boxplot(cohort_analysis$IPSW_stabil)
# 
# #trim the extreme values
# #IPSW_quantile = quantile(cohort_analysis$IPSW_stabil, c(0.02, 0.98), na.rm=T)
# #cohort_analysis$IPSW_stabil_trimmed = ifelse(cohort_analysis$IPSW_stabil < IPSW_quantile[1], IPSW_quantile[1], cohort_analysis$IPSW_stabil)
# #cohort_analysis$IPSW_stabil_trimmed = ifelse(cohort_analysis$IPSW_stabil > IPSW_quantile[2], IPSW_quantile[2], cohort_analysis$IPSW_stabil)
# #hist(cohort_analysis$IPSW_stabil_trimmed, breaks="fd")
# #boxplot(cohort_analysis$IPSW_stabil_trimmed)
# 
# #IPSW diagnostics
# cohort_analysis_complete = cohort_analysis[!is.na(cohort_analysis$PS_IPSW), ]
# bal.tab(Pat_Outreach_Flg ~ as.factor(surgery_dept), data=cohort_analysis_complete, weights=cohort_analysis_complete$IPSW_stabil, method="weighting")
# love.plot(Pat_Outreach_Flg ~ as.factor(surgery_dept), data=cohort_analysis_complete, weights=cohort_analysis_complete$IPSW_stabil, threshold=0.1)
# 
# rm(probs, model_ps_ipsw, cohort_analysis_complete)


### PRIMARY OUTCOME MODELS ###

#empty model to check for clustering by patient
model_lme4_readmit = glmer(Unplanned_Readmit_30_Day ~ (1 | PAT_ID), family=poisson(), data=cohort_analysis)
summary(model_lme4_readmit) #clustering effect, proceed with clustered robust SEs

#adjustment model, binary outcome, readmission
model_binary_adjusted_readmit = glm(Unplanned_Readmit_30_Day ~ Pat_Outreach_Flg + Pt_type_inpatient + CCI + SEX_NAME + as.factor(RACE_NAME_3cat) + ETHNICITY_NAME_2cat + LANG_NAME + as.factor(Prim_Payor) + AGE_AT_DISCHARGE + LOS_log + Geisinger_Hospital + DISCH_DT_month, family=poisson(), data=cohort_analysis)

summary(model_binary_adjusted_readmit) #confirm no overdispersion (residual deviance >> df)
exp(coef(model_binary_adjusted_readmit))
exp(robust_ci(model_binary_adjusted_readmit, group=cohort_analysis$PAT_ID)) #for no clustering, omit group argument

#ipw model, binary outcome, readmission
model_binary_ipw_readmit = glm(Unplanned_Readmit_30_Day ~ Pat_Outreach_Flg + Pt_type_inpatient + CCI + SEX_NAME + as.factor(RACE_NAME_3cat) + ETHNICITY_NAME_2cat + LANG_NAME + as.factor(Prim_Payor) + AGE_AT_DISCHARGE + LOS_log + Geisinger_Hospital + DISCH_DT_month, family=poisson(), data=cohort_analysis, weights=IPTW_stabil_trimmed*IPSW_stabil)

summary(model_binary_ipw_readmit) #confirm no overdispersion (residual deviance >> df)
exp(coef(model_binary_ipw_readmit))
exp(robust_ci(model_binary_ipw_readmit, group=cohort_analysis$PAT_ID)) #for no clustering, omit group argument
avg_comparisons(model_binary_ipw_readmit, variables="Pat_Outreach_Flg") #average marginal effect/risk difference

#adjustment model, binary outcome, ED visits
model_binary_adjusted_ed = glm(ED_30_Day ~ Pat_Outreach_Flg + Pt_type_inpatient + CCI + SEX_NAME + as.factor(RACE_NAME_3cat) + ETHNICITY_NAME_2cat + LANG_NAME + as.factor(Prim_Payor) + AGE_AT_DISCHARGE + LOS_log + Geisinger_ED + DISCH_DT_month, family=poisson(), data=cohort_analysis)

summary(model_binary_adjusted_ed) #confirm no overdispersion (residual deviance >> df)
exp(coef(model_binary_adjusted_ed))
exp(robust_ci(model_binary_adjusted_ed, group=cohort_analysis$PAT_ID)) #for no clustering, omit group argument

#ipw model, binary outcome, ED visits
model_binary_ipw_ed = glm(ED_30_Day ~ Pat_Outreach_Flg + Pt_type_inpatient + CCI + SEX_NAME + as.factor(RACE_NAME_3cat) + ETHNICITY_NAME_2cat + LANG_NAME + as.factor(Prim_Payor) + AGE_AT_DISCHARGE + LOS_log + Geisinger_ED + DISCH_DT_month, family=poisson(), data=cohort_analysis, weights=IPTW_stabil_trimmed*IPSW_stabil)

summary(model_binary_ipw_ed) #confirm no overdispersion (residual deviance >> df)
exp(coef(model_binary_ipw_ed))
exp(robust_ci(model_binary_ipw_ed, group=cohort_analysis$PAT_ID)) #for no clustering, omit group argument
avg_comparisons(model_binary_ipw_ed, variables="Pat_Outreach_Flg") #average marginal effect/risk difference

#adjustment model, count outcome, readmission
model_count_adjusted_readmit = glm(Unplanned_Readmit_30_Day_Cnt ~ Pat_Outreach_Flg + Pt_type_inpatient + CCI + SEX_NAME + as.factor(RACE_NAME_3cat) + ETHNICITY_NAME_2cat + LANG_NAME + as.factor(Prim_Payor) + AGE_AT_DISCHARGE + LOS_log + Geisinger_Hospital + DISCH_DT_month, family=poisson(), data=cohort_analysis)

summary(model_count_adjusted_readmit) #confirm no overdispersion (residual deviance >> df)
exp(coef(model_count_adjusted_readmit))
exp(robust_ci(model_count_adjusted_readmit, group=cohort_analysis$PAT_ID)) #for no clustering, use standard confint function

#ipw model, count outcome, readmission
model_count_ipw_readmit = glm(Unplanned_Readmit_30_Day_Cnt ~ Pat_Outreach_Flg + Pt_type_inpatient + CCI + SEX_NAME + as.factor(RACE_NAME_3cat) + ETHNICITY_NAME_2cat + LANG_NAME + as.factor(Prim_Payor) + AGE_AT_DISCHARGE + LOS_log + Geisinger_Hospital + DISCH_DT_month, family=poisson(), data=cohort_analysis, weights=IPTW_stabil_trimmed*IPSW_stabil)

summary(model_count_ipw_readmit) #confirm no overdispersion (residual deviance >> df)
exp(coef(model_count_ipw_readmit))
exp(robust_ci(model_count_ipw_readmit, group=cohort_analysis$PAT_ID)) #for no clustering, use standard confint function

#adjustment model, count outcome, ED visits
model_count_adjusted_ed = glm(ED_30_Day_Cnt ~ Pat_Outreach_Flg + Pt_type_inpatient + CCI + SEX_NAME + as.factor(RACE_NAME_3cat) + ETHNICITY_NAME_2cat + LANG_NAME + as.factor(Prim_Payor) + AGE_AT_DISCHARGE + LOS_log + Geisinger_ED + DISCH_DT_month, family=poisson(), data=cohort_analysis)

summary(model_count_adjusted_ed) #confirm no overdispersion (residual deviance >> df)
exp(coef(model_count_adjusted_ed))
exp(robust_ci(model_count_adjusted_ed, group=cohort_analysis$PAT_ID)) #for no clustering, use standard confint function

#ipw model, count outcome, ED visits
model_count_ipw_ed = glm(ED_30_Day_Cnt ~ Pat_Outreach_Flg + Pt_type_inpatient + CCI + SEX_NAME + as.factor(RACE_NAME_3cat) + ETHNICITY_NAME_2cat + LANG_NAME + as.factor(Prim_Payor) + AGE_AT_DISCHARGE + LOS_log + Geisinger_ED + DISCH_DT_month, family=poisson(), data=cohort_analysis, weights=IPTW_stabil_trimmed*IPSW_stabil)

summary(model_count_ipw_ed) #confirm no overdispersion (residual deviance >> df)
exp(coef(model_count_ipw_ed))
exp(robust_ci(model_count_ipw_ed, group=cohort_analysis$PAT_ID)) #for no clustering, use standard confint function


### FOREST PLOTS for PRIMARY OUTCOME MODELS ###

#binary outcome models
model_binary_ipw_readmit = glm(Unplanned_Readmit_30_Day ~ Pat_Outreach_Flg + Pt_type_inpatient + CCI + SEX_NAME + as.factor(RACE_NAME_3cat) + ETHNICITY_NAME_2cat + LANG_NAME + as.factor(Prim_Payor) + AGE_AT_DISCHARGE + LOS_log + Geisinger_Hospital + DISCH_DT_month, family=poisson(), data=cohort_analysis, weights=IPTW_stabil_trimmed*IPSW_stabil)
model_binary_ipw_ed = glm(ED_30_Day ~ Pat_Outreach_Flg + Pt_type_inpatient + CCI + SEX_NAME + as.factor(RACE_NAME_3cat) + ETHNICITY_NAME_2cat + LANG_NAME + as.factor(Prim_Payor) + AGE_AT_DISCHARGE + LOS_log + Geisinger_ED + DISCH_DT_month, family=poisson(), data=cohort_analysis, weights=IPTW_stabil_trimmed*IPSW_stabil)

#extract estimates, dropping intercept
plot_mean_readmit = exp(coef(model_binary_ipw_readmit))[-1]
plot_ci_readmit = exp(robust_ci(model_binary_ipw_readmit, group=cohort_analysis$PAT_ID))[-1, ]

#extract estimates, dropping intercept
plot_mean_ed = exp(coef(model_binary_ipw_ed))[-1]
plot_ci_ed = exp(robust_ci(model_binary_ipw_ed, group=cohort_analysis$PAT_ID))[-1, ]

forestplot(
  labeltext = c("Automated Outreach","Inpatient (ref: Out)","CCI","Male (ref: Female)","Black (ref: White)","Other (ref: White)","Hispanic","Spanish (ref: English)","Medicare (ref: GHP)","Medicaid (ref: GHP)","Other (ref: GHP)","Age","LOS (log)","Geisinger Facility","Discharge month"),
  mean = cbind(plot_mean_readmit, plot_mean_ed),
  lower = cbind(plot_ci_readmit[, 1], plot_ci_ed[, 1]),
  upper = cbind(plot_ci_readmit[, 2], plot_ci_ed[, 2]),
  zero = 1,
  xlab = "Relative Risk",
  fn.ci_norm = c(fpDrawCircleCI,fpDrawDiamondCI),
  boxsize = 0.2,
  legend = c("Readmissions","ED visits"),
  txt_gp = fpTxtGp(ticks=gpar(cex=0.75),xlab=gpar(cex=1),label=gpar(cex=0.7)),
  xticks = seq(0.6, 1.5, 0.1),
  clip = c(0.5,1.6)
)

#count outcome models
model_count_ipw_readmit = glm(Unplanned_Readmit_30_Day_Cnt ~ Pat_Outreach_Flg + Pt_type_inpatient + CCI + SEX_NAME + as.factor(RACE_NAME_3cat) + ETHNICITY_NAME_2cat + LANG_NAME + as.factor(Prim_Payor) + AGE_AT_DISCHARGE + LOS_log + Geisinger_Hospital + DISCH_DT_month, family=poisson(), data=cohort_analysis, weights=IPTW_stabil_trimmed*IPSW_stabil)
model_count_ipw_ed = glm(ED_30_Day_Cnt ~ Pat_Outreach_Flg + Pt_type_inpatient + CCI + SEX_NAME + as.factor(RACE_NAME_3cat) + ETHNICITY_NAME_2cat + LANG_NAME + as.factor(Prim_Payor) + AGE_AT_DISCHARGE + LOS_log + Geisinger_ED + DISCH_DT_month, family=poisson(), data=cohort_analysis, weights=IPTW_stabil_trimmed*IPSW_stabil)

#extract estimates, dropping intercept
plot_mean_readmit = exp(coef(model_count_ipw_readmit))[-1]
plot_ci_readmit = exp(robust_ci(model_count_ipw_readmit, group=cohort_analysis$PAT_ID))[-1, ]

#extract estimates, dropping intercept
plot_mean_ed = exp(coef(model_count_ipw_ed))[-1]
plot_ci_ed = exp(robust_ci(model_count_ipw_ed, group=cohort_analysis$PAT_ID))[-1, ]

forestplot(
  #labeltext = c("Automated Outreach","Surgery","CCI","Male (ref: Female)","Black (ref: White)","Other (ref: White)","Hispanic","Spanish (ref: English)","Medicare (ref: GHP)","Medicaid (ref: GHP)","Other (ref: GHP)","Age","LOS (log)","Geisinger Facility","Inpatient (ref: Out)","Discharge month","Outreach w/Surgery"),
  labeltext = c("Automated Outreach","Inpatient (ref: Out)","CCI","Male (ref: Female)","Black (ref: White)","Other (ref: White)","Hispanic","Spanish (ref: English)","Medicare (ref: GHP)","Medicaid (ref: GHP)","Other (ref: GHP)","Age","LOS (log)","Geisinger Facility","Discharge month","Breast surgery (ref: Bariatric)","Colorectal surgery (ref: Bariatric)","General surgery (ref: Bariatric)","Oncology surgery (ref: Bariatric)","Trauma surgery (ref: Bariatric)","Urologic surgery (ref: Bariatric)","Vascular surgery (ref: Bariatric)","Outreach w/Inpatient"),
  mean = cbind(plot_mean_readmit, plot_mean_ed),
  lower = cbind(plot_ci_readmit[, 1], plot_ci_ed[, 1]),
  upper = cbind(plot_ci_readmit[, 2], plot_ci_ed[, 2]),
  zero = 1,
  xlab = "Relative Risk",
  fn.ci_norm=c(fpDrawCircleCI,fpDrawDiamondCI),
  boxsize=0.2,
  legend=c("Readmissions","ED visits"),
  txt_gp = fpTxtGp(ticks=gpar(cex=0.75),xlab=gpar(cex=1),label=gpar(cex=0.7)),
  xticks = seq(0.6, 1.5, 0.1),
  clip = c(0.5,1.6)
)


### SURGERY SPECIFIC OUTCOME MODELS ###

#iterate through each surgery type, note for breast and oncology, need to collapse due to small n
cohort_surgery_results = data.frame("Cohort"=unique(cohort_analysis$surgery_dept), "n"=NA, "RR_readmit"=NA, "RR_readmit_lo"=NA, "RR_readmit_hi"=NA, "RD_readmit"=NA, "RD_readmit_lo"=NA, "RD_readmit_hi"=NA, "RR_ed"=NA, "RR_ed_lo"=NA, "RR_ed_hi"=NA, "RD_ed"=NA, "RD_ed_lo"=NA, "RD_ed_hi"=NA, stringsAsFactors=F)
cohort_surgery_results$Cohort[7] = "Breast or Oncology"
cohort_surgery_results = cohort_surgery_results[-8, ]

for (i in 1:nrow(cohort_surgery_results)) {
  
  #subset to this surgery
  if (cohort_surgery_results$Cohort[i] == "Breast or Oncology") {
    cohort_surgery = cohort_analysis[cohort_analysis$surgery_dept=="Breast Surgery" | cohort_analysis$surgery_dept=="Surgical Oncology", ]
  } else {
    cohort_surgery = cohort_analysis[cohort_analysis$surgery_dept==cohort_surgery_results$Cohort[i], ]
  }
  
  cohort_surgery_results$n[i] = nrow(cohort_surgery)
  
  #construct iptw
  model_ps_iptw = glm(Pat_Outreach_Flg ~ CCI + as.factor(RACE_NAME_3cat) + ETHNICITY_NAME_2cat + LANG_NAME + as.factor(Prim_Payor) + AGE_AT_DISCHARGE + LOS_log + DISCH_DT_month, family=binomial(link="logit"), data=cohort_surgery)
  
  #summary(model_ps_iptw)
  #exp(coef(model_ps_iptw))
  #exp(confint(model_ps_iptw))
  
  #predicted probabilities of automated outreach
  probs = predict(model_ps_iptw, type="response")
  cohort_surgery$PS_IPTW = NA
  cohort_surgery[names(probs), "PS_IPTW"] = probs
  #hist(cohort_surgery$PS_IPTW, breaks="fd")
  #boxplot(cohort_surgery$PS_IPTW)
  
  #inverse probability of treatment weights
  cohort_surgery$IPTW = ifelse(cohort_surgery$Pat_Outreach_Flg==1, 1/cohort_surgery$PS_IPTW, 1/(1-cohort_surgery$PS_IPTW))
  #hist(cohort_surgery$IPTW, breaks="fd")
  #boxplot(cohort_surgery$IPTW)
  
  #stabilize the weights
  cohort_surgery$IPTW_stabil = ifelse(cohort_surgery$Pat_Outreach_Flg==1, mean(cohort_surgery$Pat_Outreach_Flg, na.rm=T)/cohort_surgery$PS_IPTW, (1-mean(cohort_surgery$Pat_Outreach_Flg, na.rm=T))/(1-cohort_surgery$PS_IPTW))
  #hist(cohort_surgery$IPTW_stabil, breaks="fd")
  #boxplot(cohort_surgery$IPTW_stabil)
  
  #trim the extreme values
  IPTW_quantile = quantile(cohort_surgery$IPTW_stabil, c(0.01, 0.99), na.rm=T)
  cohort_surgery$IPTW_stabil_trimmed = ifelse(cohort_surgery$IPTW_stabil < IPTW_quantile[1], IPTW_quantile[1], cohort_surgery$IPTW_stabil)
  cohort_surgery$IPTW_stabil_trimmed = ifelse(cohort_surgery$IPTW_stabil > IPTW_quantile[2], IPTW_quantile[2], cohort_surgery$IPTW_stabil)
  #hist(cohort_surgery$IPTW_stabil_trimmed, breaks="fd")
  #boxplot(cohort_surgery$IPTW_stabil_trimmed)
  
  #IPTW diagnostics
  cohort_surgery_complete = cohort_surgery[!is.na(cohort_surgery$PS_IPTW), ]
  #bal.tab(Pat_Outreach_Flg ~ as.factor(RACE_NAME_3cat) + ETHNICITY_NAME_2cat + LANG_NAME + as.factor(Prim_Payor) + AGE_AT_DISCHARGE + LOS_log + DISCH_DT_month, data=cohort_surgery_complete, weights=cohort_surgery_complete$IPTW_stabil_trimmed, method="weighting")
  lp = love.plot(Pat_Outreach_Flg ~ as.factor(RACE_NAME_3cat) + ETHNICITY_NAME_2cat + LANG_NAME + as.factor(Prim_Payor) + AGE_AT_DISCHARGE + LOS_log + DISCH_DT_month, data=cohort_surgery_complete, weights=cohort_surgery_complete$IPTW_stabil_trimmed, threshold=0.1)
  print(lp) #workdaround needed when plotting is embedded in a for loop
  
  #outcome models
  model_binary_ipw_readmit = glm(Unplanned_Readmit_30_Day ~ Pat_Outreach_Flg + Pt_type_inpatient + CCI + SEX_NAME + as.factor(RACE_NAME_3cat) + ETHNICITY_NAME_2cat + LANG_NAME + as.factor(Prim_Payor) + AGE_AT_DISCHARGE + LOS_log + Geisinger_Hospital + DISCH_DT_month, family=poisson(), data=cohort_surgery, weights=IPTW_stabil_trimmed)
  model_binary_ipw_ed = glm(ED_30_Day ~ Pat_Outreach_Flg + Pt_type_inpatient + CCI + SEX_NAME + as.factor(RACE_NAME_3cat) + ETHNICITY_NAME_2cat + LANG_NAME + as.factor(Prim_Payor) + AGE_AT_DISCHARGE + LOS_log + Geisinger_ED + DISCH_DT_month, family=poisson(), data=cohort_surgery, weights=IPTW_stabil_trimmed)
  
  #extract estimates
  cohort_surgery_results$RR_readmit[i] = exp(coef(model_binary_ipw_readmit))["Pat_Outreach_Flg"]
  cohort_surgery_results$RR_readmit_lo[i] = exp(robust_ci(model_binary_ipw_readmit, group=cohort_surgery$PAT_ID))["Pat_Outreach_Flg", "lower"]
  cohort_surgery_results$RR_readmit_hi[i] = exp(robust_ci(model_binary_ipw_readmit, group=cohort_surgery$PAT_ID))["Pat_Outreach_Flg", "upper"]
  
  cohort_surgery_results$RR_ed[i] = exp(coef(model_binary_ipw_ed))["Pat_Outreach_Flg"]
  cohort_surgery_results$RR_ed_lo[i] = exp(robust_ci(model_binary_ipw_ed, group=cohort_surgery$PAT_ID))["Pat_Outreach_Flg", "lower"]
  cohort_surgery_results$RR_ed_hi[i] = exp(robust_ci(model_binary_ipw_ed, group=cohort_surgery$PAT_ID))["Pat_Outreach_Flg", "upper"]
  
  #average marginal effect/risk difference
  cohort_surgery_results$RD_readmit[i] = avg_comparisons(model_binary_ipw_readmit, variables="Pat_Outreach_Flg")$estimate
  cohort_surgery_results$RD_readmit_lo[i] = avg_comparisons(model_binary_ipw_readmit, variables="Pat_Outreach_Flg")$conf.low
  cohort_surgery_results$RD_readmit_hi[i] = avg_comparisons(model_binary_ipw_readmit, variables="Pat_Outreach_Flg")$conf.high
  
  cohort_surgery_results$RD_ed[i] = avg_comparisons(model_binary_ipw_ed, variables="Pat_Outreach_Flg")$estimate
  cohort_surgery_results$RD_ed_lo[i] = avg_comparisons(model_binary_ipw_ed, variables="Pat_Outreach_Flg")$conf.low
  cohort_surgery_results$RD_ed_hi[i] = avg_comparisons(model_binary_ipw_ed, variables="Pat_Outreach_Flg")$conf.high
}
rm(i, probs, model_ps_iptw, cohort_surgery_complete, lp)

#plot results
forestplot(
  labeltext = paste(cohort_surgery_results$Cohort),
  mean = cohort_surgery_results[c("RR_readmit","RR_ed")],
  lower = cohort_surgery_results[c("RR_readmit_lo","RR_ed_lo")],
  upper = cohort_surgery_results[c("RR_readmit_hi","RR_ed_hi")],
  zero = 1,
  xlab = "Relative Risk",
  fn.ci_norm = c(fpDrawCircleCI,fpDrawDiamondCI),
  boxsize = 0.2,
  legend = c("Readmissions","ED visits"),
  txt_gp = fpTxtGp(ticks=gpar(cex=0.75),xlab=gpar(cex=1),label=gpar(cex=0.7)),
  xticks = seq(0.3, 2.0, 0.3),
  clip = c(0.2,2.1)
)


### INPATIENT/OUTPATIENT SPECIFIC OUTCOME MODELS ###

#separate inpatient/outpatient
cohort_inpatient_results = data.frame("Cohort"=c("Inpatient","Outpatient"), "n"=NA, "RR_readmit"=NA, "RR_readmit_lo"=NA, "RR_readmit_hi"=NA, "RD_readmit"=NA, "RD_readmit_lo"=NA, "RD_readmit_hi"=NA, "RR_ed"=NA, "RR_ed_lo"=NA, "RR_ed_hi"=NA, "RD_ed"=NA, "RD_ed_lo"=NA, "RD_ed_hi"=NA, stringsAsFactors=F)

for (i in 1:nrow(cohort_inpatient_results)) {
  
  #subset to this surgery
  if (cohort_inpatient_results$Cohort[i] == "Inpatient") {
    cohort_inpatient = cohort_analysis[cohort_analysis$Pt_type_inpatient==1, ]
  } else {
    cohort_inpatient = cohort_analysis[cohort_analysis$Pt_type_outpatient==1, ]
  }
  
  cohort_inpatient_results$n[i] = nrow(cohort_inpatient)
  
  #construct iptw
  model_ps_iptw = glm(Pat_Outreach_Flg ~ CCI + as.factor(RACE_NAME_3cat) + ETHNICITY_NAME_2cat + LANG_NAME + as.factor(Prim_Payor) + AGE_AT_DISCHARGE + LOS_log + DISCH_DT_month, family=binomial(link="logit"), data=cohort_inpatient)
  
  #summary(model_ps_iptw)
  #exp(coef(model_ps_iptw))
  #exp(confint(model_ps_iptw))
  
  #predicted probabilities of automated outreach
  probs = predict(model_ps_iptw, type="response")
  cohort_inpatient$PS_IPTW = NA
  cohort_inpatient[names(probs), "PS_IPTW"] = probs
  #hist(cohort_inpatient$PS_IPTW, breaks="fd")
  #boxplot(cohort_inpatient$PS_IPTW)
  
  #inverse probability of treatment weights
  cohort_inpatient$IPTW = ifelse(cohort_inpatient$Pat_Outreach_Flg==1, 1/cohort_inpatient$PS_IPTW, 1/(1-cohort_inpatient$PS_IPTW))
  #hist(cohort_inpatient$IPTW, breaks="fd")
  #boxplot(cohort_inpatient$IPTW)
  
  #stabilize the weights
  cohort_inpatient$IPTW_stabil = ifelse(cohort_inpatient$Pat_Outreach_Flg==1, mean(cohort_inpatient$Pat_Outreach_Flg, na.rm=T)/cohort_inpatient$PS_IPTW, (1-mean(cohort_inpatient$Pat_Outreach_Flg, na.rm=T))/(1-cohort_inpatient$PS_IPTW))
  #hist(cohort_inpatient$IPTW_stabil, breaks="fd")
  #boxplot(cohort_inpatient$IPTW_stabil)
  
  #trim the extreme values
  IPTW_quantile = quantile(cohort_inpatient$IPTW_stabil, c(0.01, 0.99), na.rm=T)
  cohort_inpatient$IPTW_stabil_trimmed = ifelse(cohort_inpatient$IPTW_stabil < IPTW_quantile[1], IPTW_quantile[1], cohort_inpatient$IPTW_stabil)
  cohort_inpatient$IPTW_stabil_trimmed = ifelse(cohort_inpatient$IPTW_stabil > IPTW_quantile[2], IPTW_quantile[2], cohort_inpatient$IPTW_stabil)
  #hist(cohort_inpatient$IPTW_stabil_trimmed, breaks="fd")
  #boxplot(cohort_inpatient$IPTW_stabil_trimmed)
  
  #IPTW diagnostics
  cohort_inpatient_complete = cohort_inpatient[!is.na(cohort_inpatient$PS_IPTW), ]
  #bal.tab(Pat_Outreach_Flg ~ as.factor(RACE_NAME_3cat) + ETHNICITY_NAME_2cat + LANG_NAME + as.factor(Prim_Payor) + AGE_AT_DISCHARGE + LOS_log + DISCH_DT_month, data=cohort_inpatient_complete, weights=cohort_inpatient_complete$IPTW_stabil_trimmed, method="weighting")
  lp = love.plot(Pat_Outreach_Flg ~ as.factor(RACE_NAME_3cat) + ETHNICITY_NAME_2cat + LANG_NAME + as.factor(Prim_Payor) + AGE_AT_DISCHARGE + LOS_log + DISCH_DT_month, data=cohort_inpatient_complete, weights=cohort_inpatient_complete$IPTW_stabil_trimmed, threshold=0.1)
  print(lp) #workdaround needed when plotting is embedded in a for loop
  
  #outcome models
  model_binary_ipw_readmit = glm(Unplanned_Readmit_30_Day ~ Pat_Outreach_Flg + CCI + SEX_NAME + as.factor(RACE_NAME_3cat) + ETHNICITY_NAME_2cat + LANG_NAME + as.factor(Prim_Payor) + AGE_AT_DISCHARGE + LOS_log + Geisinger_Hospital + DISCH_DT_month, family=poisson(), data=cohort_inpatient, weights=IPTW_stabil_trimmed)
  model_binary_ipw_ed = glm(ED_30_Day ~ Pat_Outreach_Flg + CCI + SEX_NAME + as.factor(RACE_NAME_3cat) + ETHNICITY_NAME_2cat + LANG_NAME + as.factor(Prim_Payor) + AGE_AT_DISCHARGE + LOS_log + Geisinger_ED + DISCH_DT_month, family=poisson(), data=cohort_inpatient, weights=IPTW_stabil_trimmed)
  
  #extract estimates
  cohort_inpatient_results$RR_readmit[i] = exp(coef(model_binary_ipw_readmit))["Pat_Outreach_Flg"]
  cohort_inpatient_results$RR_readmit_lo[i] = exp(robust_ci(model_binary_ipw_readmit, group=cohort_inpatient$PAT_ID))["Pat_Outreach_Flg", "lower"]
  cohort_inpatient_results$RR_readmit_hi[i] = exp(robust_ci(model_binary_ipw_readmit, group=cohort_inpatient$PAT_ID))["Pat_Outreach_Flg", "upper"]
  
  cohort_inpatient_results$RR_ed[i] = exp(coef(model_binary_ipw_ed))["Pat_Outreach_Flg"]
  cohort_inpatient_results$RR_ed_lo[i] = exp(robust_ci(model_binary_ipw_ed, group=cohort_inpatient$PAT_ID))["Pat_Outreach_Flg", "lower"]
  cohort_inpatient_results$RR_ed_hi[i] = exp(robust_ci(model_binary_ipw_ed, group=cohort_inpatient$PAT_ID))["Pat_Outreach_Flg", "upper"]
  
  #average marginal effect/risk difference
  cohort_inpatient_results$RD_readmit[i] = avg_comparisons(model_binary_ipw_readmit, variables="Pat_Outreach_Flg")$estimate
  cohort_inpatient_results$RD_readmit_lo[i] = avg_comparisons(model_binary_ipw_readmit, variables="Pat_Outreach_Flg")$conf.low
  cohort_inpatient_results$RD_readmit_hi[i] = avg_comparisons(model_binary_ipw_readmit, variables="Pat_Outreach_Flg")$conf.high
  
  cohort_inpatient_results$RD_ed[i] = avg_comparisons(model_binary_ipw_ed, variables="Pat_Outreach_Flg")$estimate
  cohort_inpatient_results$RD_ed_lo[i] = avg_comparisons(model_binary_ipw_ed, variables="Pat_Outreach_Flg")$conf.low
  cohort_inpatient_results$RD_ed_hi[i] = avg_comparisons(model_binary_ipw_ed, variables="Pat_Outreach_Flg")$conf.high
}
rm(i, probs, model_ps_iptw, cohort_inpatient_complete, lp)

#plot results
forestplot(
  labeltext = c(paste(cohort_inpatient_results$Cohort),NA), #some sort of bug with the forestplot code and only two groups, adding a third group fixes it
  mean = rbind(cohort_inpatient_results[c("RR_readmit","RR_ed")],NA),
  lower = rbind(cohort_inpatient_results[c("RR_readmit_lo","RR_ed_lo")],NA),
  upper = rbind(cohort_inpatient_results[c("RR_readmit_hi","RR_ed_hi")],NA),
  zero = 1,
  xlab = "Relative Risk",
  fn.ci_norm = c(fpDrawCircleCI,fpDrawDiamondCI),
  boxsize = 0.1,
  legend = c("Readmissions","ED visits"),
  txt_gp = fpTxtGp(ticks=gpar(cex=0.75),xlab=gpar(cex=1),label=gpar(cex=0.7)),
  xticks = seq(0.3, 2.0, 0.3),
  clip = c(0.2,2.1)
)

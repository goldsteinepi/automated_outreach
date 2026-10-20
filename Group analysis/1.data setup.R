#Set-up analysis file 
rm(list=ls()) 

library(tidyverse) 
library(janitor)
library(readxl)
library(writexl)
library(fixest)
library(ggplot2)
library(writexl)
library(patchwork)
library(dplyr)
library(did)
library(modelsummary)
library(ggfixest)
library(skimr)

setwd("C:/Users/kkilburn/OneDrive - Geisinger System Services/GPE and BIT Shared Directory-GPE shared data (PHI) - Documents/Automated Discharge/Data/")
#NDG: setwd("~/OneDrive - Geisinger System Services/GPE and BIT Shared Directory-GPE shared data (PHI) - Documents/Automated Discharge/Data/")

# =============================================================================
# 1. Data set-up
# =============================================================================

#---------- Pull in data from Clarity, PDC dashboard, and department names
df_raw <- read_csv("Data Pulls/All discharges/1117GPEKK1205250504260735-Post Surgery Outreach-v1.08-Y12969.csv", 
              )  %>% clean_names() %>%
          mutate(
            admsn_dt = as.Date(admsn_dt, format = "%m/%d/%Y"),
            disch_dt = as.Date(disch_dt, format = "%m/%d/%Y"),
            disch_prov_id = as.numeric(disch_prov_id)
          )


df_pdc <- read_excel("Data Pulls/Dashboards/PDC Dashboard full 20-Feb-2026.xlsx") %>%
                        clean_names() %>%
          mutate(
            rr_score = ifelse(!readm_risk_score %in% c("46183", "46341"), # I believe these are null values
              readm_risk_score,NA))

df_matches <- read_excel("R/Department matchings_2.19.26.xlsx") %>% clean_names()  #original department names

df_rematches <- read_excel("R/Department matchings_4.4.26.xlsx") %>% #pooling some department names
          clean_names() %>% rename(department = prim_workday_dept_name) %>% 
          select(department, pooled_department, outreach_new)


#------- Pull in corrected 2023 data for ED visits as it was missing from earlier pull
df_23edvisits <- read_csv("Data Pulls/All discharges/1117GPEKK1205250514262020-Post Surgery Outreach-v1.09.csv")  %>% 
        clean_names() %>% 
        mutate(
          admsn_dt = as.Date(admsn_dt, format = "%m/%d/%Y"),
          disch_dt = as.Date(disch_dt, format = "%m/%d/%Y"),
          disch_prov_id = as.numeric(disch_prov_id)
        ) %>%  
        mutate(year=year(disch_dt), 
               month = month(disch_dt),
               ed23 = ed_30_day_cnt,
               eddate = frst_30_day_ed_adm_dt) %>%
        filter(year == 2023) %>%
        select("pat_id", "pat_enc_csn_id", year, disch_dt, ed23, eddate)


# =============================================================================
# 2. Data merging and recoding
# =============================================================================

#--------Merging '23 ED visits with main data frame


df_main <- df_raw %>%
  left_join(df_23edvisits, by = c("pat_id", "pat_enc_csn_id", "disch_dt")) %>%
  mutate(
    ed_30_day_cnt = if_else(!is.na(year) & year == 2023, ed23, ed_30_day_cnt), #year is missing for non-2023 so need to specify it this way to keep other years
    frst_30_day_ed_adm_dt = if_else(!is.na(year) & year == 2023, eddate, frst_30_day_ed_adm_dt)
  ) %>%
  select(-ed23, -eddate, -year)



#-------Setting up Analysis data frame with discharge/patient level data

##look at all relevant variables = sex_name, month, or_dash_case, readmit, planned_readmit, unplanned_readmit, ed_30, department

  df_analysis <- merge(df_main, df_matches, by = "prim_workday_dept_name") %>%
          mutate(year=year(disch_dt), 
                 month = month(disch_dt),
                 or_dash_case = if_else(is.na(or_dash_case), 0, or_dash_case),
                 readmit_30_day_cnt = if_else(is.na(readmit_30_day_cnt), 0, readmit_30_day_cnt ),
                 planned_readmit_30_day_cnt = if_else(is.na(planned_readmit_30_day_cnt), 0, planned_readmit_30_day_cnt),
                 unplanned_readmit_30_day_cnt = if_else(is.na(unplanned_readmit_30_day_cnt), 0, unplanned_readmit_30_day_cnt),
                 ed_30_day_cnt = if_else(is.na(ed_30_day_cnt), 0, ed_30_day_cnt ),
                 female = if_else(sex_name == "Female", 1, 0),
                 pat_outreach_flg = ifelse(disch_dt < as.Date("2024-05-09"), 0, pat_outreach_flg)) %>%
          rename(department = prim_workday_dept_name,
                 adt_class = adt_pat_class_name) %>%
         # filter(year>2023 & year <2026) %>% # limiting to 24 months 2024 to 2025 
          filter((lang_name == "English" | lang_name == "Spanish") & age_at_discharge>= 18) #filtering for criteria to receive outreach)


        # # Checking Neal's patient CSNs for department correct
        # check<- df_analysis%>% filter(pat_enc_csn_id %in% c(310168846980,
        #                                                     310178716511,310142837287,
        #                                                     310158168945, 310165255406,
        #                                                     310167768318))
        #checking out discharges that had outreach not in treatment group
          # otherstocheck<- df_analysis %>% filter(treatment == 0 & pat_outreach_flg==1)
          # crossover <- merge(otherstocheck, df_pdc, by = "pat_enc_csn_id")
              ## ONLY 16/3234 discharges that have outreach from control group show up in the PDC data 

#--------Building some group level aggregated data/ desc statistics        

  df_group <- df_analysis %>%
    group_by(department, outreach_name, treatment, month, year) %>%
    summarise(
      n = n(),
      pct_or_dash_case = round(mean(or_dash_case, na.rm = TRUE) * 100, 0),
      .groups = "drop"
    ) %>%
    group_by(department) %>%
    add_tally(wt = n, name = "total_dep") %>%  # sum of n within department
    ungroup() %>%
    mutate(total_n = sum(n),  # overall total across all departments/treatments/months
           ratio_dis = n/total_dep)

         
#--------Final individual sample data set for analysis       

  df_analysis <- merge(df_analysis, df_group)


# #--------Main group data set for analysis
#   
#   df_test <- df_analysis %>%
#     mutate(readmit_yn = ifelse(readmit_30_day_cnt>0, 1, 0),
#            edvisit_yn = ifelse(ed_30_day_cnt>0, 1, 0),
#            unplanned_yn = ifelse(unplanned_readmit_30_day_cnt>0, 1, 0))%>%
#     group_by(department, outreach_name, treatment, month, year) %>%
#     summarise(
#       n = n(),
#       reads = sum(readmit_yn, na.rm = TRUE),
#       unplanned_reads = sum(unplanned_yn, na.rm = TRUE),
#       edvisits = sum(edvisit_yn, na.rm = TRUE),
#       flg_outreach = sum(replace_na(pat_outreach_flg,0), na.rm = TRUE),
#       pct_reads = round((reads / n) * 100, 0),
#       pct_unplanned_reads = round((unplanned_reads / n) * 100, 0),
#       pct_edvisits = round((edvisits/n)* 100, 0),
#       pct_outreach = round((flg_outreach/n)*100, 0),
#       avg_cci = mean(cci, na.rm = TRUE),
#       pct_or_dash_case = round(mean(or_dash_case, na.rm = TRUE) * 100, 0),
#       avg_bmi = mean(bmi, na.rm= TRUE),
#       avg_age = mean(age_at_discharge, na.rm = TRUE),
#       female_total = sum(female, na.rm = TRUE),
#       pct_female = round((female_total / n)*100, 0),
#       los = mean(los, na.rm = TRUE),
#       pre_rrs = mean(admn_readm_risk_score, na.rm = TRUE),
#       dis_rrs = mean(pre_dsch_readm_risk_score, na.rm = TRUE),
#       post_rrs = mean(post_dsch_readm_risk_score, na.rm = TRUE),
#       .groups = "drop"
#     ) %>% 
#     group_by(department) %>%
#     add_tally(wt = n, name = "total_dep") %>%  # sum of n within department
#     ungroup() %>%
#     mutate(total_dis = sum(n),  # overall total across all departments/treatments/months
#            ratio_dis = n/total_dep) %>%
#     mutate(m24 = case_when(
#       year == 2024 ~ month,
#       year == 2025 ~ month + 12,
#       TRUE ~ NA_real_ ),  # ensure month index is integer 1..24 for main analysis time frame
#       treatment = factor(treatment, levels = c(0, 1),
#                          labels = c("Control", "Treatment")),
#       dept_id = match(department, unique(department))) %>%
#      mutate(month_index = (year - 2023) * 12 + month)   # new index variable for all months 2023-2026 we have
#   write_xlsx(df_test, "R/group analysis/df_test.xlsx")
#   
#   overall <-df_test %>% filter(month_index<25) %>% group_by(outreach_name)%>% 
#     summarise(pct_ureads = mean(pct_unplanned_reads),
#               mean_los = mean(los))
#   print(overall)
# 
# #--------------Reviewing main data elements
#   #Checking for event-study units with n<10
#  smallns <- df_test %>% filter(n < 10) %>% mutate(smalln_flag = 1)
#  table(smallns$department, smallns$n)
#  df_test <- left_join(df_test, smallns)  # Reassigned departments with many small ns to other locations (see design choices doc)
 
    #Pool departments with a lot of small ns - bringing in df_rematches to reassign some departments
       
       df_analysis <- merge(df_analysis, df_rematches, by = "department")  %>% 
         mutate(department = coalesce(pooled_department, department),
                outreach_name = coalesce(outreach_new, outreach_name)) 
       write_xlsx(df_analysis, "R/group analysis/df_analysis.xlsx")    

#--------Main group data set for analysis
       
       df_test <- df_analysis %>%
         mutate(readmit_yn = ifelse(readmit_30_day_cnt>0, 1, 0),
                edvisit_yn = ifelse(ed_30_day_cnt>0, 1, 0),
                unplanned_yn = ifelse(unplanned_readmit_30_day_cnt>0, 1, 0))%>%
         group_by(department, outreach_name, treatment, month, year) %>%
         summarise(
           n = n(),
           reads = sum(readmit_yn, na.rm = TRUE),
           unplanned_reads = sum(unplanned_yn, na.rm = TRUE),
           edvisits = sum(edvisit_yn, na.rm = TRUE),
           flg_outreach = sum(pat_outreach_flg, na.rm = TRUE),
           pct_reads = round((reads / n) * 100, 0),
           pct_unplanned_reads = round((unplanned_reads / n) * 100, 0),
           pct_edvisits = round((edvisits/n)* 100, 0),
           pct_outreach = round((flg_outreach/n)*100, 0),
           avg_cci = mean(cci, na.rm = TRUE),
           pct_or_dash_case = round(mean(or_dash_case, na.rm = TRUE) * 100, 0),
           avg_bmi = mean(bmi, na.rm= TRUE),
           avg_age = mean(age_at_discharge, na.rm = TRUE),
           female_total = sum(female, na.rm = TRUE),
           pct_female = round((female_total / n)*100, 0),
           los = mean(los, na.rm = TRUE),
           pre_rrs = mean(admn_readm_risk_score, na.rm = TRUE),
           dis_rrs = mean(pre_dsch_readm_risk_score, na.rm = TRUE),
           post_rrs = mean(post_dsch_readm_risk_score, na.rm = TRUE),
           .groups = "drop"
         ) %>% 
         group_by(department) %>%
         add_tally(wt = n, name = "total_dep") %>%  # sum of n within department
         ungroup() %>%
         mutate(total_dis = sum(n),  # overall total across all departments/treatments/months
                ratio_dis = n/total_dep) %>%
         mutate(m24 = case_when(
           year == 2024 ~ month,
           year == 2025 ~ month + 12,
           TRUE ~ NA_real_ ),  # ensure month index is integer 1..24 for main analysis time frame
           treatment = factor(treatment, levels = c(0, 1),
                              labels = c("Control", "Treatment")),
           dept_id = match(department, unique(department))) %>%
         mutate(month_index = (year - 2023) * 12 + month)   # new index variable for all months 2023-2026 we have
       write_xlsx(df_test, "R/group analysis/df_test.xlsx")    
       
       #check to ensure pooling worked
       smallns <- df_test %>% filter(n < 10) %>% mutate(smalln_flag = 1)
       table(smallns$department, smallns$n)
       
       #additionally remove months where n<=5
       df_test <- df_test %>% filter(n>5)
       table(df_test$department)
       
       
 ##-----OPTIONAL - change up the comparison group by removing low risk departments/ same day surgery ones?
        overall <-df_test %>% filter(month_index<25) %>% group_by(outreach_name)%>% 
          summarise(pct_ureads = mean(pct_unplanned_reads),
                      mean_los = mean(los))
        print(overall)
    #Option 1    
         df_test_drop2 <- df_test %>% 
           filter(!outreach_name %in% c("Ophthamology surgery", "Pain management"))
    #Option 2    
         df_test_drop4 <- df_test %>% 
           filter(!outreach_name %in% c("Ophthamology surgery", "Pain management", "Oral Surgery", "Otolaryngology surgery")) 

         
# =============================================================================
###Reviewing data
# =============================================================================
         
df_analysis %>% group_by(year) %>% skim()
table(df_analysis$department)
unique(df_analysis$department)
table(df_analysis$or_dash_case)
table(df_analysis$treatment)
table (df_analysis$pat_outreach_flg)
table (df_analysis$female)
summary(df_analysis$pre_dsch_readm_risk_score)
summary(df_analysis$post_dsch_readm_risk_score)
summary(df_analysis$admn_readm_risk_score)

# =============================================================================
#  3. Data preparation for models
# =============================================================================
#-----------choose which comparison sample (df_test, df_test_drop2, df_test_drop4)
        # df_test<-df_test_drop2 #drop 2 comparison groups
        #df_test<-df_test_drop4 #drop 4 comparison groups

#--------Setting up staggered DiD groups


    df_did <- df_test %>%
      mutate(
        # First-treatment period (g): 0 = never treated, all other months are base off starting date of Jan '23
        g = case_when(
          treatment == "Control"                                         ~ 0L,
          outreach_name %in% c("Colorectal Surgery", "Urology Surgery") ~ 17L,  #switches on May '24
          outreach_name == "Vascular Surgery"                           ~ 31L,  #switches on July '25
          treatment == "Treatment"                                       ~ 25L, #switches on Jan '25
          TRUE                                                           ~ NA_integer_
        ),
        month_index = (year - 2023) * 12 + month,  # new index variable for all months 2023-2026 we have
        # Month-of-year FE (1-12) for seasonality
        month_fe = factor(((month_index - 1) %% 12) + 1L),
        # Relative event time (NA for never-treated)
        rel_time = ifelse(g > 0, month_index - g, NA_real_)
      )
    df_did$g <- as.numeric(df_did$g)
    # Data check: no units with missing g
    stopifnot(!anyNA(df_did$g))

#--------Setting up ITS groups
    
    # Remove Colorectal, Urology, and Vascular so all treated units share m24 = 13
    df_its <- df_test %>%
      filter(!outreach_name %in% c("Colorectal Surgery", "Urology Surgery",
                                   "Vascular Surgery")) %>%  
      filter(year>2023 & year <2026) %>% # limiting to 24 months 2024 to 2025 
      mutate(
        treated    = as.integer(treatment == "Treatment"),
        post       = as.integer(m24 >= 13),
        time_since = if_else(m24 >= 13, m24 - 13L, 0L),  # 0 pre, 1,2,3... post
        time_pre   = if_else(m24 < 13, m24, 12L),  # 1,2,3... post is 12
        month_fe   = factor(((m24 - 1) %% 12) + 1L)       # month-of-year FE
      )
    
    #Sample after exclusions:
    print(table(df_its$treatment))
    #Departments retained:
    print(sort(unique(df_its$outreach_name)))

######### Reviewing matches between our departments and PDC dashboard categories
  #df_check <- merge(df_analysis, df_pdc, by = "pat_enc_csn_id") 
    
    
###############################################################################
# Subgroup analysis by Inpatient status
###############################################################################    
    
# #-------Filter for inpatient discharges per Cory's request (subgroup analysis)
      df_analysis_ips <- merge(df_analysis, df_group) %>%
        filter(adt_class == "Inpatient")
        #filter(adt_class %in%  c("Inpatient", "SORU (23 hr Surgical Overnight)"))
    
# #-------Set up two data frames - one for surgery group analysis and one for dept group analysis
    #concern is that inpatients are a much smaller group (~20% of total discharges)
    
    df_test_ips_dept <- df_analysis_ips %>%
      mutate(readmit_yn = ifelse(readmit_30_day_cnt>0, 1, 0),
             edvisit_yn = ifelse(ed_30_day_cnt>0, 1, 0),
             unplanned_yn = ifelse(unplanned_readmit_30_day_cnt>0, 1, 0))%>%
      group_by(department, outreach_name, treatment, month, year) %>%
      summarise(
        n = n(),
        reads = sum(readmit_yn, na.rm = TRUE),
        unplanned_reads = sum(unplanned_yn, na.rm = TRUE),
        edvisits = sum(edvisit_yn, na.rm = TRUE),
        flg_outreach = sum(pat_outreach_flg, na.rm = TRUE),
        pct_reads = round((reads / n) * 100, 0),
        pct_unplanned_reads = round((unplanned_reads / n) * 100, 0),
        pct_edvisits = round((edvisits/n)* 100, 0),
        pct_outreach = round((flg_outreach/n)*100, 0),
        avg_cci = mean(cci, na.rm = TRUE),
        pct_or_dash_case = round(mean(or_dash_case, na.rm = TRUE) * 100, 0),
        avg_bmi = mean(bmi, na.rm= TRUE),
        avg_age = mean(age_at_discharge, na.rm = TRUE),
        female_total = sum(female, na.rm = TRUE),
        pct_female = round((female_total / n)*100, 0),
        los = mean(los, na.rm = TRUE),
        pre_rrs = mean(admn_readm_risk_score, na.rm = TRUE),
        dis_rrs = mean(pre_dsch_readm_risk_score, na.rm = TRUE),
        post_rrs = mean(post_dsch_readm_risk_score, na.rm = TRUE),
        .groups = "drop"
      ) %>% 
      group_by(department) %>%
      add_tally(wt = n, name = "total_dep") %>%  # sum of n within department
      ungroup() %>%
      mutate(total_dis = sum(n),  # overall total across all departments/treatments/months
             ratio_dis = n/total_dep) %>%
      mutate(m24 = case_when(
        year == 2024 ~ month,
        year == 2025 ~ month + 12,
        TRUE ~ NA_real_ ),  # ensure month index is integer 1..24 for main analysis time frame
        treatment = factor(treatment, levels = c(0, 1),
                           labels = c("Control", "Treatment")),
        dept_id = match(department, unique(department))) %>%
      mutate(month_index = (year - 2023) * 12 + month) %>%   # new index variable for all months 2023-2026 we have
    filter(n>5)
    
    table(df_test_ips_dept$department)
    
    df_test_ips_out <- df_analysis_ips %>%
      mutate(readmit_yn = ifelse(readmit_30_day_cnt>0, 1, 0),
             edvisit_yn = ifelse(ed_30_day_cnt>0, 1, 0),
             unplanned_yn = ifelse(unplanned_readmit_30_day_cnt>0, 1, 0))%>%
      group_by(outreach_name, treatment, month, year) %>%
      summarise(
        n = n(),
        reads = sum(readmit_yn, na.rm = TRUE),
        unplanned_reads = sum(unplanned_yn, na.rm = TRUE),
        edvisits = sum(edvisit_yn, na.rm = TRUE),
        flg_outreach = sum(pat_outreach_flg, na.rm = TRUE),
        pct_reads = round((reads / n) * 100, 0),
        pct_unplanned_reads = round((unplanned_reads / n) * 100, 0),
        pct_edvisits = round((edvisits/n)* 100, 0),
        pct_outreach = round((flg_outreach/n)*100, 0),
        avg_cci = mean(cci, na.rm = TRUE),
        pct_or_dash_case = round(mean(or_dash_case, na.rm = TRUE) * 100, 0),
        avg_bmi = mean(bmi, na.rm= TRUE),
        avg_age = mean(age_at_discharge, na.rm = TRUE),
        female_total = sum(female, na.rm = TRUE),
        pct_female = round((female_total / n)*100, 0),
        los = mean(los, na.rm = TRUE),
        pre_rrs = mean(admn_readm_risk_score, na.rm = TRUE),
        dis_rrs = mean(pre_dsch_readm_risk_score, na.rm = TRUE),
        post_rrs = mean(post_dsch_readm_risk_score, na.rm = TRUE),
        .groups = "drop"
      ) %>% 
      group_by(outreach_name) %>%
      add_tally(wt = n, name = "total_dep") %>%  # sum of n within department
      ungroup() %>%
      mutate(total_dis = sum(n),  # overall total across all departments/treatments/months
             ratio_dis = n/total_dep) %>%
      mutate(m24 = case_when(
        year == 2024 ~ month,
        year == 2025 ~ month + 12,
        TRUE ~ NA_real_ ),  # ensure month index is integer 1..24 for main analysis time frame
        treatment = factor(treatment, levels = c(0, 1),
                           labels = c("Control", "Treatment")),
        out_id = match(outreach_name, unique(outreach_name))) %>%
      mutate(month_index = (year - 2023) * 12 + month)   # new index variable for all months 2023-2026 we have
    
    
    table(df_test_ips_out$outreach_name)
    
    
    #--------Setting up ITS groups
    
    # Remove Colorectal, Urology, and Vascular so all treated units share m24 = 13
    df_its_ips_out <- df_test_ips_out %>%
      filter(!outreach_name %in% c("Colorectal Surgery", "Urology Surgery",
                                   "Vascular Surgery")) %>%  
      filter(year>2023 & year <2026) %>% # limiting to 24 months 2024 to 2025 
      mutate(
        treated    = as.integer(treatment == "Treatment"),
        post       = as.integer(m24 >= 13),
        time_since = if_else(m24 >= 13, m24 - 13L, 0L),  # 0 pre, 1,2,3... post
        time_pre   = if_else(m24 < 13, m24, 12L),  # 1,2,3... post is 12
        month_fe   = factor(((m24 - 1) %% 12) + 1L)       # month-of-year FE
      )
    #Sample after exclusions:
    print(table(df_its_ips_out$treatment))
    
    df_its_ips_dept <- df_test_ips_dept %>%
      filter(!outreach_name %in% c("Colorectal Surgery", "Urology Surgery",
                                   "Vascular Surgery")) %>%  
      filter(year>2023 & year <2026) %>% # limiting to 24 months 2024 to 2025 
      mutate(
        treated    = as.integer(treatment == "Treatment"),
        post       = as.integer(m24 >= 13),
        time_since = if_else(m24 >= 13, m24 - 13L, 0L),  # 0 pre, 1,2,3... post
        time_pre   = if_else(m24 < 13, m24, 12L),  # 1,2,3... post is 12
        month_fe   = factor(((m24 - 1) %% 12) + 1L)       # month-of-year FE
      )
    #Sample after exclusions:
    print(table(df_its_ips_dept$treatment))
   
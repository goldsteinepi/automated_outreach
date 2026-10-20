#Set-up analysis file 

library(tidyverse) 
library(janitor)
library(readxl)
#install.packages("tidytext") 
#install.packages("writexl") 

library(flextable)
library(officer)
library(tidytext)
library(gt)
library(writexl)

setwd("C:/Users/kkilburn/OneDrive - Geisinger System Services/GPE and BIT Shared Directory-GPE shared data (PHI) - Documents/Automated Discharge/Data/")

############################################################################
# 1. Descriptive statistics - by outreach group and department
############################################################################

#library(dplyr)

# Restrict to pre-period only 
# df_pre <- df_did %>%
#   filter(month_index<25)


df_pre <- df_did %>%
  filter(month_index <if_else( outreach_name %in% c("Urology Surgery", "Colorectal Surgery"),
        17L, 25L))


# 1. Counts by Treatment vs Control
# counts_by_group <- df_pre %>%
#   group_by(treatment) %>%
#   summarise(
#     n_rows = n(),
#     total_discharges = sum(n, na.rm = TRUE),
#     .groups = "drop"
#   ) %>%
#   arrange(treatment)

counts_by_treatment <- df_pre %>%
  filter(
    (outreach_name %in% c("Urology Surgery", "Colorectal Surgery") &
       month_index >= 5 & month_index <= 16) |
      (!outreach_name %in% c("Urology Surgery", "Colorectal Surgery") &
         month_index >= 13 & month_index <= 24)
  ) %>%
  group_by(treatment) %>%
  summarise(
    n_rows = n(),
    total_discharges = sum(n, na.rm = TRUE),
    .groups = "drop"
  ) %>%
 arrange(treatment)


# 2. Average LOS and Unplanned Readmissions by Outreach Group
summary_by_outreach <- df_pre %>%
  filter((outreach_name %in% c("Urology Surgery", "Colorectal Surgery") &
       month_index >= 5 & month_index <= 16) |
      (!outreach_name %in% c("Urology Surgery", "Colorectal Surgery") &
         month_index >= 13 & month_index <= 24)) %>%
  group_by(outreach_name, treatment) %>%
  summarise(
    avg_los = mean(los, na.rm = TRUE),
    sd_los = sd(los, na.rm = TRUE),
    avg_unplanned_reads = mean(unplanned_reads, na.rm = TRUE),
    sd_unplanned_reads = sd(unplanned_reads, na.rm = TRUE),
    avg_edvisits = mean(edvisits, na.rm = TRUE),
    sd_edvisits = sd(edvisits, na.rm = TRUE),
    total_discharges = sum(n, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  arrange(treatment)

# View results
counts_by_treatment
summary_by_outreach



############################################################################
# 2. Plots for descriptive statistics - by outreach group and department
############################################################################

# 
# table<-summary_by_outreach %>%
#   arrange(treatment, outreach_name) %>%
#   rename(
#     `Outreach Group` = outreach_name,
#     `Treatment Group` = treatment,
#     `Avg Length of Stay` = avg_los,
#     `Avg Unplanned Readmissions` = avg_unplanned_reads,
#     `Avg ED Visits` = avg_edvisits,
#     `Discharges in prior 12 months` = total_discharges
#   ) %>%
#   gt() %>%
#   fmt_number(
#     columns = c(`Avg Length of Stay`, `Avg Unplanned Readmissions`, `Avg ED Visits`),
#     decimals = 2
#   ) %>%
#   tab_header(
#     title = "Pre-Period Descriptive Statistics",
#     subtitle = "By Outreach Group and Treatment Status"
#   ) %>%
#   opt_table_outline()
# 
# gtsave(table, "R/group analysis/outputs/desctable.html")

# Create table-ready data frame
table_df <- summary_by_outreach %>%
  arrange(treatment, outreach_name) %>%
  mutate(
    los = sprintf("%.2f (%.2f)", avg_los, sd_los),
    unplanned_reads = sprintf("%.2f (%.2f)", avg_unplanned_reads, sd_unplanned_reads),
    ed_visits = sprintf("%.2f (%.2f)", avg_edvisits, sd_edvisits)
  ) %>%
  select(
    outreach_name,
    treatment,
    los,
    unplanned_reads,
    ed_visits,
    total_discharges
  ) %>%
  rename(
    `Outreach Group` = outreach_name,
    `Treatment Group` = treatment,
    `Length of Stay Mean (SD)` = los,
    `Unplanned Readmissions Mean (SD)` = unplanned_reads,
    `ED Visits Mean (SD)` = ed_visits,
    `Discharges in prior 12 months` = total_discharges
  )

# Export as CSV
readr::write_csv(
  table_df,
  "R/group analysis/outputs/desctable.csv"
)

# Create gt table for HTML export
table <- table_df %>%
  gt() %>%
  tab_header(
    title = "Pre-Period Descriptive Statistics",
    subtitle = "By Outreach Group and Treatment Status"
  ) %>%
  opt_table_outline()

# Export as HTML
gtsave(table, "R/group analysis/outputs/desctable.html")



###############
# Remove all unnecessary Data and Objects
rm(list = setdiff(ls(), c("df_test", "df_analysis", "df_main", 
                          "df_matches", "df_pdc", "df_did",
                          "df_its", "df_sa", "df_cs")))


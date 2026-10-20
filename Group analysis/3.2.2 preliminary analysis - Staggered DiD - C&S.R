# =============================================================================
# Staggered Difference-in-Differences / Event-Study Analysis
# Estimator: Calloway & Sant'anna (did)
# Outcomes:  readmissions, unplanned readmissions, edvisits
# =============================================================================

required_pkgs <- c("fixest", "ggplot2", "dplyr", "writexl", "patchwork")

to_install <- required_pkgs[!sapply(required_pkgs, requireNamespace, quietly = TRUE)]
if (length(to_install)) install.packages(to_install)

library(fixest)
library(did)
library(ggplot2)
library(writexl)
library(patchwork)
library(dplyr)
library(modelsummary)


# =============================================================================
# Summary
# =============================================================================

# Design: Staggered adoption DiD, dept x month panel (Jan 2024-Dec 2025) - Sun & Abraham
# Cohorts: g=5  (Colorectal/Urology pilot)
#          g=13 (Jan-25 roll-out, all other treated)
#          g=19 (Vascular, Jul-25)
#          g=0  (never-treated controls)
# 
# MODEL SPEC
#   did::att_gt()
#   Treatment timing: g (cohort)
#   Time index: month_index
#   Unit id: dept_id
#   Estimator: doubly robust
#   Covariates: avg_age + avg_cci + pct_female
#   SE: clustered at dept_id

# MAIN OUTPUTS
#   Plots : 
#   Tables: 
# 
# CAVEATS
#   - Colorectal cohort has only 1 dept; cohort-level ATT imprecise.
#   - Vascular cohort (g=19) has short post-period; interpret cautiously.
#   - Many zero outcomes visible in sample; consider other 
#     logit regression models as a robustness check.
#   - Parallel trends assumed with review of pre-trend plots and Wald test.

# =============================================================================
# 1. Data preparation & Review
# =============================================================================

  
  df_cs <- df_sa # NEED to run 3.2 preliminary analysis - Staggered DiD_S&A.R

# =============================================================================
# 2. Calloway and Sant'anna models
# =============================================================================
  ## Callaway & Sant'Anna estimators using did ##
  table(df_cs$g)
    ##NOTES: not including covariates as they can introduce huge CIs, control group is nevertreated
  run_attgt <- function(df,
                        outcome,
                        tname = "month_index",
                        idname = "dept_id",
                        gname = "g",
                        cluster = "dept_id",
                        est_method = "dr",
                        allow_unbalanced = TRUE) {
  
  # Estimate group-time ATTs
  att <- att_gt(
    yname = outcome,
    tname = tname,
    idname = idname,
    gname = gname,
    data = df,
    clustervars = cluster,
    est_method = est_method,
    allow_unbalanced_panel = allow_unbalanced
  )
  
  # Aggregate
  att_group   <- aggte(att, type = "group")
  att_overall <- aggte(att, type = "simple")
  att_dynamic <- aggte(att, type = "dynamic")

  
  # Return all objects in a list
  return(list(
    att_raw = att,
    group   = att_group,
    overall = att_overall,
    dynamic = att_dynamic
  ))
}



#========== Results ==========
 

  results_reads <- run_attgt(df = df_cs,
                             outcome = "pct_reads")
   #att_group <- aggte(results_reads, type = "group") 
   
   #ggdid(results_reads)
  
  results_ureads <- run_attgt(df = df_cs,
                             outcome = "pct_unplanned_reads")
  results_ed <- run_attgt(df = df_cs,
                          outcome = "pct_edvisits")
  
  results_ureads$att_raw
  results_ureads$group
  results_ureads$overall
  results_ureads$dynamic
  
  ggdid(results_ureads$dynamic)
# =============================================================================
# 2. Exporting Results
# =============================================================================
# Function for exporting to Excel 
export_attgt_excel <- function(res,
                               path = "attgt_results.xlsx") {
                               #scale = 100) {
  
  wb <- openxlsx::createWorkbook()
  
  # Overall
  ov <- data.frame(
    att     = res$overall$overall.att , 
    se      = res$overall$overall.se , 
    ci_low  = (res$overall$overall.att - 1.96 * res$overall$overall.se) , 
    ci_high = (res$overall$overall.att + 1.96 * res$overall$overall.se)  # * scale
  )
  openxlsx::addWorksheet(wb, "overall")
  openxlsx::writeData(wb, "overall", ov)
  
  # Group
  grp <- data.frame(
    group   = res$group$group,
    att     = res$group$att.group, 
    se      = res$group$se.group, 
    ci_low  = (res$group$att.group - 1.96 * res$group$se.group), 
    ci_high = (res$group$att.group + 1.96 * res$group$se.group) 
  )
  openxlsx::addWorksheet(wb, "group")
  openxlsx::writeData(wb, "group", grp)
  
  # Dynamic
  dyn <- data.frame(
    event_time = res$dynamic$egt,
    att        = res$dynamic$att.egt , 
    se         = res$dynamic$se.egt , 
    ci_low     = (res$dynamic$att.egt - 1.96 * res$dynamic$se.egt) , 
    ci_high    = (res$dynamic$att.egt + 1.96 * res$dynamic$se.egt)  
  )
  openxlsx::addWorksheet(wb, "dynamic")
  openxlsx::writeData(wb, "dynamic", dyn)
  
  openxlsx::saveWorkbook(wb, path, overwrite = TRUE)
  message("Saved: ", path)
}

#========== Export to Excel ==========
  export_attgt_excel(results_reads, "R/group analysis/outputs/pct_reads_results.xlsx")
  export_attgt_excel(results_ureads, "R/group analysis/outputs/pct_unplanned_reads_results.xlsx")
  export_attgt_excel(results_ureads, "R/group analysis/outputs/pct_edvisits_results.xlsx")
  
  
  
#Functions for Plotting group and dynamic results 

  plot_dynamic_atts <- function(res_dynamic,
                                title_suffix = NULL,
                                file = NULL,
                                width = 7, height = 4, dpi = 300) {
    p <- ggdid(res_dynamic)
    if (!is.null(title_suffix)) {
      p <- p + ggtitle(paste0("Dynamic ATTs (Event Study)", " ", title_suffix))
    }
    if (!is.null(file)) {
      ggsave(filename = file, plot = p, width = width, height = height, dpi = dpi)
    }
    p
  }
  
  plot_group_atts <- function(res_group,
                              title_suffix = NULL,
                              file = NULL,
                              width = 7, height = 4, dpi = 300) {
    p <- ggdid(res_group)
    if (!is.null(title_suffix)) {
      p <- p + ggtitle(paste0("Group ATTs", " ", title_suffix))
    }
    if (!is.null(file)) {
      ggsave(filename = file, plot = p, width = width, height = height, dpi = dpi)
    }
    p
  }
  

#========== Save Plots ==========

#--------Readmissions
p_group_reads <- plot_group_atts(
  results_reads$group,
  title_suffix = "– % readmissions",
  file = "R/group analysis/outputs/group_reads.png"
)

p_dynamic_reads <- plot_dynamic_atts(
  results_reads$dynamic,
  title_suffix = "– % readmissions",
  file = "R/group analysis/outputs/dynamic_reads.png"
)

#---------Unplanned Readmissions
p_group_ureads <- plot_group_atts(
  results_ureads$group,
  title_suffix = "– % unplanned readmissions",
  file = "R/group analysis/outputs/group_ureads.png"
)

p_dynamic_ureads <- plot_dynamic_atts(
  results_ureads$dynamic,
  title_suffix = "– % unplanned readmissions",
  file = "R/group analysis/outputs/dynamic_ureads.png"
)
#---------ED visits
p_group_ed <- plot_group_atts(
  results_ed$group,
  title_suffix = "– % ED Visits",
  file = "R/group analysis/outputs/group_ed.png"
)

p_dynamic_ed <- plot_dynamic_atts(
  results_ed$dynamic,
  title_suffix = "– % ED Visits",
  file = "R/group analysis/outputs/dynamic_ed.png"
)

#=============================================================
# Heterogeneity Results #
#==============================================================

#--------- Separate results by outreach level #
  df_cs <- df_cs %>%
    mutate(
      g_int = ifelse(
        trt_outreach_group == "Treatment, No High Outreach",
        g,
        ifelse(treatment == "Treatment", NA, 0)
      )
    )
  table(df_cs$g_int)
  
  het_run_attgt <- function(df,
                            outcome,
                            tname = "month_index",
                            idname = "dept_id",
                            gname = "g_int",
                            cluster = "dept_id",
                            est_method = "dr",
                            allow_unbalanced = TRUE) {
    
    # Estimate group-time ATTs
    att <- att_gt(
      yname = outcome,
      tname = tname,
      idname = idname,
      gname = gname,
      data = df,
      clustervars = cluster,
      est_method = est_method,
      allow_unbalanced_panel = allow_unbalanced
    )
    
    # Aggregate
    att_group   <- aggte(att, type = "group")
    att_overall <- aggte(att, type = "simple")
    att_dynamic <- aggte(att, type = "dynamic")
    
    # Return all objects in a list
    return(list(
      att_raw = att,
      group   = att_group,
      overall = att_overall,
      dynamic = att_dynamic
    ))
  }
  
  hetoutreach_results_ureads <- het_run_attgt(df = df_cs,
                                      outcome = "pct_unplanned_reads")
  
  #--------- Separate results by CCI median #
    #Just changing g_int to CCI
  df_cs <- df_cs %>%
    mutate(
      g_int = ifelse(
        scenario_label == "Scenario 1: Treatment, CCI avg above median",
        g,
        ifelse(treatment == "Treatment", NA, 0)  # removing treatment groups from analysis
      )
    )
  table(df_cs$g_int)
  
  
  hetcci_results_ureads <- het_run_attgt(df= df_cs,
                                      outcome = "pct_unplanned_reads")
  
  #Pretty table
  att_group1   <- att_gt(data = df_cs,
                           yname  = "pct_unplanned_reads",
                           tname = "month_index",
                           idname = "dept_id",
                           gname = "g",
                           cluster = "dept_id",
                           est_method = "dr",
                           allow_unbalanced = TRUE)
  att_group1 <- aggte(att_group1, type = "group")
  att_group2   <- att_gt(data = df_cs,
                         yname  = "pct_unplanned_reads",
                         tname = "month_index",
                         idname = "dept_id",
                         gname = "g_int_outreach",
                         cluster = "dept_id",
                         est_method = "dr",
                         allow_unbalanced = TRUE)
  att_group2   <- aggte(att_group2, type = "group")
  att_group3   <- att_gt(data = df_cs,
                         yname  = "pct_unplanned_reads",
                         tname = "month_index",
                         idname = "dept_id",
                         gname = "g_int_cci",
                         cluster = "dept_id",
                         est_method = "dr",
                         allow_unbalanced = TRUE)
  att_group3   <- aggte(att_group3, type = "group")
  modelsummary(list(att_group1, att_group2, att_group3), 
               conf_level = 0.90, statistic = "conf.int", 
               stars = c("*" = .10, "**" = .05, "***" = .01), 
               fmt = 2)
  

# Remove all unnecessary Data and Objects
  rm(list = setdiff(ls(), c("df_test", "df_analysis", "df_main", 
                            "df_matches", "df_pdc", "df_did",
                            "df_its", "df_sa", "df_cs",
                            "df_its_ips_out", "df_test_ips_out",
                            "df_its_ips_dept", "df_test_ips_dept")))
  
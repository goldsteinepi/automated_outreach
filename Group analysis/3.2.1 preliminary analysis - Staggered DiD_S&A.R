# =============================================================================
# Staggered Difference-in-Differences / Event-Study Analysis
# Estimator: Sun & Abraham (fixest::sunab)
# Outcomes:  readmissions, unplanned readmissions, edvisits
# =============================================================================

required_pkgs <- c("fixest", "ggplot2", "dplyr", "writexl", "patchwork", "modelsummary")

to_install <- required_pkgs[!sapply(required_pkgs, requireNamespace, quietly = TRUE)]
if (length(to_install)) install.packages(to_install)

library(fixest)
library(ggplot2)
library(writexl)
library(patchwork)
library(dplyr)
library(did)
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
#   fixest::sunab(g, month_index, ref.p = -1)
#   FE: dept_id + month_index (calendar period) + month_fe (seasonality)
#   SE: clustered at dept_id
# 
# OUTPUTS
#   Plots : sunab_pct_reads.png, sunab_pct_edvisits.png, sunab_overview.png
#   Tables: sunab_results.xlsx (5 sheets)
#   Cov   : covariate_trends.png
# 
# CAVEATS
#   - Colorectal cohort has only 1 dept; cohort-level ATT imprecise.
#   - Vascular cohort (g=19) has short post-period; interpret cautiously.
#   - Many zero outcomes visible in sample; consider beta/fractional
#     logit regression as a robustness check.
#   - Parallel trends assumed; validate via pre-trend plots and Wald p.

# -----------------------------------------------------------------------------
# 1. Data preparation & Review
# -----------------------------------------------------------------------------
#--------Setting up staggered DiD groups

  df_sa <- df_did #%>%
    #filter(is.na(smalln_flag))  # GETTING rid of dept-group <10 obs
         
   
    # Data check: no units with missing g
    stopifnot(!anyNA(df_sa$g))
    
    #Unit counts by cohort (g)
    print(table(df_sa$g))

#----Adding in outreach level groupings & CCI groupings
  avgoutreach <- df_sa %>%
    group_by(department) %>%
    summarize(post_outreach = mean(pct_outreach[rel_time >= 0], na.rm = TRUE),
    totavg_cci = mean(avg_cci, na.rm = TRUE)) %>%
    ungroup()
  
  #What's the median cci average?
    summary(df_sa$avg_cci)  
    
    df_sa <- merge(df_sa, avgoutreach, by = "department") %>%
      group_by(department) %>%
      mutate(high_outreach = ifelse(is.na(post_outreach), 0,
                                    ifelse(post_outreach > 80, 1, 0)), # binary variable to identify groups that received high outreach on avg
             cci_median = ifelse(totavg_cci > 3.7, 1, 0)) %>%  # binary variable for groups with avg cci above median
      ungroup() %>%
  
  #create a new grouping to distinguish treatment and control groups by any high outreach in the post-period
  mutate(
    trt_outreach_group = case_when(
      treatment == "Treatment" & high_outreach == 1 ~ "Treatment, High Outreach",
      treatment == "Treatment" & high_outreach == 0 ~ "Treatment, No High Outreach",
      treatment == "Control"  ~ "Control, No High Outreach"
    ),
    cci_median_scenario = case_when(cci_median ==1 & treatment == "Treatment" ~ 1L, 
                                cci_median ==0 & treatment =="Treatment" ~ 2L,
                                 treatment == "Control"  ~ 3L),
           scenario_label = factor(
             cci_median_scenario,
             levels = c(1, 2, 3),
             labels = c(
               "Scenario 1: Treatment, CCI avg above median",
               "Scenario 2: Treatment, CCI avg below median",
               "Scenario 3: Control" )))

# -----Time groups check (event-time count per bin)
    support_tbl <- df_sa %>%
      filter(!is.na(rel_time)) %>%
      count(rel_time, name = "n_obs") %>%
      arrange(rel_time)
    
    write_xlsx(list(support = support_tbl),
               "R/group analysis/outputs/support_check.xlsx")
    
    p_support <- ggplot(support_tbl, aes(rel_time, n_obs)) +
      geom_col(fill = "steelblue") +
      labs(title = "Observation Count by Event Time",
           x = "Event Time (months relative to adoption)",
           y = "N observations") 
    ggsave("R/group analysis/outputs/support_check.png", p_support, width=8, height=4, dpi=150)

# =============================================================================
# 2. Sun & Abraham models
# =============================================================================
#----OLS - Unadjusted
run_sunab <- function(outcome, data = df_sa) {
  fml <- as.formula(
    paste0(outcome, " ~ sunab(g, month_index, ref.p = -1) | dept_id +  month_fe  ")
  )
  feols(fml, data = data, cluster = ~dept_id)
}
#----OLS - Adjusted
run_asunab <- function(outcome, data = df_sa) {
  fml <- as.formula(
    paste0(outcome, " ~ sunab(g, month_index, ref.p = -1)", 
    " + avg_cci + avg_age + pct_female", 
   # "+ pct_outreach", --> only include if interested in mediation by fidelity to treatment
    " | dept_id +  month_fe"
  ))
  feols(fml, data = data, cluster = ~dept_id)
}

#----Poisson - Unadjusted
run_sunab_poisson <- function(outcome_count, denominator= "n", data = df_sa) {
  data$log_denom <- log(data[[denominator]])
  fml <- as.formula(
    paste0(outcome_count, 
           "~ sunab(g, month_index, ref.p = -1)" ,
           "+ offset(log_denom)",
           "| dept_id +  month_fe")
  )
  feglm(fml, data = data, cluster = ~dept_id, family = poisson())
}

#---Poisson - Adjusted
run_sunab_apoisson <- function(outcome_count, denominator= "n", data = df_sa) {
  data$log_denom <- log(data[[denominator]])
  fml <- as.formula(
    paste0(outcome_count, 
           "~ sunab(g, month_index, ref.p = -1)" ,
           " + avg_cci + avg_age + pct_female", 
          # "+ pct_outreach", --> only include if interested in mediation by fidelity to treatment
          "+ offset(log_denom)",
           "| dept_id +  month_fe")
  )
  feglm(fml, data = data, cluster = ~dept_id, family = poisson())
}
  #========== Readmissions ==========
  sa_reads <- run_sunab("pct_reads")
  sa_areads <- run_asunab("pct_reads")
  sa_reads_pois <- run_sunab_poisson("reads")
  sa_areads_pois <- run_sunab_apoisson("reads")
  
  
  #========== Unplanned Readmissions ==========
  sa_u_reads <- run_sunab("pct_unplanned_reads")
  sa_u_areads <- run_asunab("pct_unplanned_reads")
  sa_u_reads_pois <- run_sunab_poisson("unplanned_reads")
  sa_u_areads_pois <- run_sunab_apoisson("unplanned_reads")
  
  #========== ED visits ==========
  sa_edvisits <- run_sunab("pct_edvisits")
  sa_aedvisits <- run_asunab("pct_edvisits")
  # sa_edvisits_pois <- run_sunab_poisson("edvisits")
  # sa_aedvisits_pois <- run_sunab_apoisson("edvisits")
  # 
  #       #see which FE levels were actually estimated 
  #       fe_est <- fixest::fixef(sa_edvisits_pois)
  #       lapply(fe_est, function(x) head(names(x)))
  
# -- 3a. Overall ATT - OLS ---------------------------------------------
  att_sa_reads    <- summary(sa_reads,    agg = "att")
  att_sa_u_reads    <- summary(sa_u_reads,    agg = "att")
  att_sa_edvisits <- summary(sa_edvisits, agg = "att")
  att_sa_areads    <- summary(sa_areads,    agg = "att")
  att_sa_u_areads    <- summary(sa_u_areads,    agg = "att")
  att_sa_aedvisits <- summary(sa_aedvisits, agg = "att")
  
   print(att_sa_reads)
   print(att_sa_u_reads)
   print(att_sa_edvisits)
   print(att_sa_areads)
   print(att_sa_u_areads)
   print(att_sa_aedvisits)
   
   #Pretty table
   modelsummary(list(att_sa_u_reads, att_sa_u_areads), 
                conf_level = 0.90, statistic = "conf.int", 
                stars = c("*" = .10, "**" = .05, "***" = .01), 
                fmt = 2)
   
# -- 3b. Overall ATT - POISSON -----------------------------------
   patt_sa_reads    <- summary(sa_reads_pois,    agg = "att")
   patt_sa_u_reads    <- summary(sa_u_reads_pois,    agg = "att")
   #patt_sa_edvisits <- summary(sa_edvisits_pois, agg = "att")
   patt_sa_areads    <- summary(sa_areads_pois,    agg = "att")
   patt_sa_u_areads    <- summary(sa_u_areads_pois,    agg = "att")
   # patt_sa_aedvisits <- summary(sa_aedvisits_pois, agg = "att")
   
   print(patt_sa_reads)
   print(patt_sa_u_reads)
   # print(patt_sa_edvisits)
   print(patt_sa_areads)
   print(patt_sa_u_areads)
   # print(patt_sa_aedvisits)
   
# -- 3c. Cohort-specific ATTs - OLS --------------------
  att_cohort_reads    <- summary(sa_reads,    agg = "cohort")
   att_cohort_u_reads    <- summary(sa_u_reads,    agg = "cohort")
   att_cohort_edvisits <- summary(sa_edvisits, agg = "cohort")
   att_cohort_areads    <- summary(sa_areads,    agg = "cohort")
   att_cohort_u_areads    <- summary(sa_u_areads,    agg = "cohort")
   att_cohort_aedvisits <- summary(sa_aedvisits, agg = "cohort")
   
    print(att_cohort_reads)
    print(att_cohort_u_reads)
    print(att_cohort_edvisits)
    print(att_cohort_areads)
    print(att_cohort_u_areads)
    print(att_cohort_aedvisits)
    
# -- 3d. Cohort-specific ATTs - Poisson --------------------
    patt_cohort_reads    <- summary(sa_reads_pois,    agg = "cohort")
    patt_cohort_u_reads    <- summary(sa_u_reads_pois,    agg = "cohort")
    # patt_cohort_edvisits <- summary(sa_edvisits_pois, agg = "cohort")
    patt_cohort_areads    <- summary(sa_areads_pois,    agg = "cohort")
    patt_cohort_u_areads    <- summary(sa_u_areads_pois,    agg = "cohort")
    # patt_cohort_aedvisits <- summary(sa_aedvisits_pois, agg = "cohort")
    
    print(patt_cohort_reads)
    print(patt_cohort_u_reads)
    # print(patt_cohort_edvisits)
    print(patt_cohort_areads)
    print(patt_cohort_u_areads)
    # print(patt_cohort_aedvisits)
    
##---- Tables for presentation -----------------------##
    # etable(att_sa_u_reads, coefstat = "confint", ci = 0.90,
    #        signif.code = c("***" = 0.01, "**" = 0.05, "*" = 0.10),
    #        notes = "Confidence intervals shown are 90%.")

    # etable(att_sa_u_areads, coefstat = "confint", ci = 0.90,
    #        signif.code = c("***" = 0.01, "**" = 0.05, "*" = 0.10),
    #        notes = "Confidence intervals shown are 90%.")
    # irr_model <- broom::tidy(patt_sa_u_reads) %>%
    #   mutate(
    #     irr = exp(estimate),
    #     irr_low = exp(estimate - 1.96 * std.error),
    #     irr_high = exp(estimate + 1.96 * std.error)
    #   )
    # etable(patt_sa_u_reads, coefstat = "confint", ci = 0.90,
    #        signif.code = c("***" = 0.01, "**" = 0.05, "*" = 0.10),
    #        notes = "Confidence intervals shown are 90%.")
    # etable(att_cohort_u_reads, coefstat = "confint", ci = 0.90,
    #        signif.code = c("***" = 0.01, "**" = 0.05, "*" = 0.10),
    #        notes = "Confidence intervals shown are 90%.")
    # etable(att_cohort_u_areads, coefstat = "confint", ci = 0.90,
    #        signif.code = c("***" = 0.01, "**" = 0.05, "*" = 0.10),
    #        notes = "Confidence intervals shown are 90%.")
    irr_model_cohort <- broom::tidy(patt_cohort_u_reads) %>%
      mutate(
        irr = exp(estimate),
        irr_low = exp(estimate - 1.96 * std.error),
        irr_high = exp(estimate + 1.96 * std.error)
      )
    # etable(patt_cohort_u_reads, coefstat = "confint", ci = 0.90,
    #        signif.code = c("***" = 0.01, "**" = 0.05, "*" = 0.10),
    #        notes = "Confidence intervals shown are 90%.")

# -- 3e. Plots via fixest::iplot ----------------------------------------
# iplot() is  built-in event-study plotter for sunab models, it aggregates cohort-specific estimates to one point per period.

# png("R/group analysis/outputs/sunab_pct_reads.png", width = 900, height = 500, res = 100)
# iplot(sa_reads,
#       main  = "Sun & Abraham Event Study: % Readmissions",
#       xlab  = "Periods relative to treatment",
#       ylab  = "Coefficient (95% CI)")
# dev.off()
# 
# png("R/group analysis/outputs/sunab_pct_u_reads.png", width = 900, height = 500, res = 100)
# iplot(sa_u_reads,
#       main  = "Sun & Abraham Event Study: % Readmissions",
#       xlab  = "Periods relative to treatment",
#       ylab  = "Coefficient (95% CI)")
# dev.off()
# 
# png("R/group analysis/outputs/sunab_pct_edvisits.png", width = 900, height = 500, res = 100)
# iplot(sa_edvisits,
#       main  = "Sun & Abraham Event Study: % ED Visits",
#       xlab  = "Periods relative to treatment",
#       ylab  = "Coefficient (95% CI)")
# dev.off()
# 
# # Combined overview (both panels in one PNG)
# png("R/group analysis/outputs/sunab_overview.png", width = 900, height = 1000, res = 100)
# par(mfrow = c(3, 1), mar = c(4, 4, 3, 1))
# iplot(sa_reads,
#       main = "Sun & Abraham: % Readmissions",
#       xlab = "Periods relative to treatment",
#       ylab = "Coefficient (95% CI)")
# iplot(sa_u_reads,
#       main = "Sun & Abraham: % Unplanned Readmissions",
#       xlab = "Periods relative to treatment",
#       ylab = "Coefficient (95% CI)")
# iplot(sa_edvisits,
#       main = "Sun & Abraham: % ED Visits",
#       xlab = "Periods relative to treatment",
#       ylab = "Coefficient (95% CI)")
# dev.off()

    
# =============================================================================
# 4. Pre-trend tests
# =============================================================================
    
      b    <- coef(sa_reads)
    V  <- vcov(sa_reads)
    
    
    # keep Sun–Abraham coefficients only
    keep <- grepl("^sunab::", names(b))
    b    <- b[keep]
    
      se   <- sqrt(diag(sa_reads$cov.scaled))
      keep <- grepl("^month_index::", names(b))
      b    <- b[keep]
      se   <- se[keep]
      
      periods   <- suppressWarnings(as.integer(sub("^month_index::([-0-9]+):.*", "\\1", names(b))))
      pre_idx   <- which(!is.na(periods) & periods < 0)
      pre_names <- names(b)[pre_idx]
      
      wt <- wald(sa_reads, keep = pre_names)
      wt <- tryCatch(
        wald(sa_reads, keep = pre_names),
        error = function(e) { message("Wald test skipped: ", e$message); NULL }
      )
    
    pretrend_test <- function(mod, outcome_label) {
      b    <- mod$coefficients
      se   <- sqrt(diag(mod$cov.scaled))
      keep <- grepl("^month_index::", names(b))
      b    <- b[keep]
      se   <- se[keep]
      
      periods   <- suppressWarnings(as.integer(sub("^month_index::([-0-9]+):.*", "\\1", names(b))))
      pre_idx   <- which(!is.na(periods) & periods < 0)
      pre_names <- names(b)[pre_idx]
      
      wt <- tryCatch(
        wald(mod, keep = pre_names),
        error = function(e) { message("Wald test skipped: ", e$message); NULL }
      )
      
      ind_tbl <- data.frame(
        outcome  = outcome_label,
        period   = periods[pre_idx],
        estimate = unname(b[pre_idx]),
        se       = unname(se[pre_idx]),
        stringsAsFactors = FALSE
      ) %>%
        group_by(outcome, period) %>%
        summarise(
          estimate = mean(estimate),
          se       = sqrt(mean(se^2)),
          .groups  = "drop"
        ) %>%
        mutate(
          ci_lo   = estimate - 1.96 * se,
          ci_hi   = estimate + 1.96 * se,
          p_value = 2 * pnorm(-abs(estimate / se))
        ) %>%
        arrange(period)
      
      list(wald = wt, individual = ind_tbl)
    }
    
    pt_reads    <- pretrend_test(sa_reads,    "pct_reads")
    pt_u_reads    <- pretrend_test(sa_u_reads,    "pct_unplanned_reads")
    pt_edvisits <- pretrend_test(sa_edvisits, "pct_edvisits")
    
    #=== Pre-trend Joint Wald Test - pct_reads 
    print(pt_reads$wald)
    print(pt_u_reads$wald)
    
    #=== Pre-trend Joint Wald Test - pct_edvisits 
    print(pt_edvisits$wald)
    
    #==== Individual Pre-period Leads - pct_reads 
    print(pt_reads$individual)
    print(pt_u_reads$individual)
    
    #=== Individual Pre-period Leads - pct_edvisits 
    print(pt_edvisits$individual)
    
    
    # Export pre-trend results
    wald_summary <- data.frame(
      outcome = c("pct_reads", "pct_unplanned_reads", "pct_edvisits"),
      F_stat  = c(pt_reads$wald$stat, pt_u_reads$wald$stat, pt_edvisits$wald$stat),
      p_value = c(pt_reads$wald$p,  pt_u_reads$wald$p,   pt_edvisits$wald$p)
    )
    
    lead_table <- bind_rows(pt_reads$individual, pt_u_reads$individual, pt_edvisits$individual)
    
    write_xlsx(
      list(
        "PreTrend_Wald"       = wald_summary,
        "PreTrend_Individual" = lead_table
      ),
      path = "R/group analysis/outputs/pretrend_results.xlsx"
    )

    
# =============================================================================
# 5. Event Study Plots
# =============================================================================    
    
# ---- 5a. Event-study plot - OLS --------------------------------

  #Extract coefficients from Sunab models first
  extract_es <- function(mod) {
    broom::tidy(mod, conf.int = TRUE) |>
      dplyr::filter(stringr::str_detect(term, "^month_index::\\s*-?\\d+$")) |>
      dplyr::mutate(rel_time = as.numeric(stringr::str_extract(term, "-?\\d+$"))) |>
      dplyr::select(rel_time, estimate, std.error, statistic, p.value, conf.low, conf.high) |>
      dplyr::arrange(rel_time)
  }
  
  
  es_reads <- extract_es(sa_reads) %>% mutate(outcome = "pct_reads")
  es_u_reads <- extract_es(sa_u_reads) %>% mutate(outcome = "pct_uplanned_reads")
  es_ed    <- extract_es(sa_edvisits)    %>% mutate(outcome = "pct_edvisits")
  es_all   <- bind_rows(es_reads, es_u_reads, es_ed)
  es_areads <- extract_es(sa_areads) %>% mutate(outcome = "pct_reads")
  es_u_areads <- extract_es(sa_u_areads) %>% mutate(outcome = "pct_unplanned_reads")
  es_aed    <- extract_es(sa_aedvisits)    %>% mutate(outcome = "pct_edvisits")
  es_adjusted_all   <- bind_rows(es_areads, es_u_areads, es_aed)

  #Plot coefficients for event study visual
  plot_es_compare <- function(es_unadj, es_adj, title_str) {
  # Add a model label to each input and combine
  es_unadj2 <- es_unadj %>% mutate(model = "Unadjusted")
  es_adj2   <- es_adj   %>% mutate(model = "Adjusted")
  
  es_both <- bind_rows(es_unadj2, es_adj2)
  
  ggplot(es_both, aes(rel_time, estimate, colour = model, fill = model, linetype = model)) +
    geom_hline(yintercept = 0, linetype = "dashed", colour = "grey50") +
    geom_vline(xintercept = 0, linetype = "dotted", colour = "red") +
    # CI ribbon per model
    geom_ribbon(aes(ymin = conf.low, ymax = conf.high), alpha = 0.15, colour = NA) +
    # Line & points per model
    geom_line(linewidth = 1) +
    geom_point(size = 2) +
    # Highlight significant points per model
    geom_point(
      data = es_both %>% filter(p.value < 0.05),
      size = 3, shape = 21, colour = "black"
    ) +
    scale_x_continuous(breaks = scales::pretty_breaks(n = 12)) +
    scale_colour_manual(values = c("Unadjusted" = "steelblue", "Adjusted" = "coral4")) +
    scale_fill_manual(values   = c("Unadjusted" = "steelblue", "Adjusted" = "coral4")) +
    scale_linetype_manual(values= c("Unadjusted" = "solid",     "Adjusted" = "dotdash")) +
    labs(
      title    = title_str,
      subtitle = "Sun & Abraham event study | FEOLS | 95% CI | Clustered SE (dept_id)",
      x        = "Months relative to adoption",
      y        = "Coefficient estimate",
      colour   = NULL, fill = NULL, linetype = NULL,
      caption  = "Black-rimmed points: p < 0.05. Vertical dotted line: adoption."
    ) }

  p_es_reads <- plot_es_compare(es_reads, es_areads, "Event Study – % Readmissions")
  p_es_u_reads <- plot_es_compare(es_u_reads, es_u_areads, "Event Study – % Unplanned Readmissions")
  p_es_ed    <- plot_es_compare(es_ed, es_aed,   "Event Study – % ED Visits")
  ggsave("R/group analysis/outputs/event_study_pct_reads.png",    p_es_reads, width=10, height=5, dpi=150)
  ggsave("R/group analysis/outputs/event_study_pct_ureads.png",    p_es_u_reads, width=10, height=5, dpi=150)
  ggsave("R/group analysis/outputs/event_study_pct_edvisits.png", p_es_ed,    width=10, height=5, dpi=150)
  ggsave("R/group analysis/outputs/event_study_combined.png",
         p_es_reads / p_es_u_reads / p_es_ed, width=10, height=10, dpi=150)

  

 

# ---- 5c. Event-study plot - Poisson --------------------------------

#Extract coefficients from Sunab models first
  es_reads_pois <- extract_es(sa_reads_pois) %>% mutate(outcome = "reads")
  es_u_reads_pois <- extract_es(sa_u_reads_pois) %>% mutate(outcome = "unplanned_reads")
  #es_ed_pois    <- extract_es(sa_edvisits_pois)    %>% mutate(outcome = "edvisits")
  es_all_pois   <- bind_rows(es_reads_pois, es_u_reads_pois)
  es_areads_pois <- extract_es(sa_areads_pois) %>% mutate(outcome = "reads")
  es_u_areads_pois <- extract_es(sa_u_areads_pois) %>% mutate(outcome = "unplanned_reads")
  #es_aed_pois    <- extract_es(sa_aedvisits_pois)    %>% mutate(outcome = "edvisits")
  es_adjusted_all_pois   <- bind_rows(es_areads_pois, es_u_areads_pois)

#Plot coefficients for event study visual
  plot_es_compare_pois <- function(es_unadj, es_adj, title_str) {
  # Add a model label to each input and combine
  es_unadj2 <- es_unadj %>% mutate(model = "Unadjusted")
  es_adj2   <- es_adj   %>% mutate(model = "Adjusted")
  
  es_both <- bind_rows(es_unadj2, es_adj2)
  
  ggplot(es_both, aes(rel_time, estimate, colour = model, fill = model, linetype = model)) +
    geom_hline(yintercept = 0, linetype = "dashed", colour = "grey50") +
    geom_vline(xintercept = 0, linetype = "dotted", colour = "red") +
    # CI ribbon per model
    geom_ribbon(aes(ymin = conf.low, ymax = conf.high), alpha = 0.15, colour = NA) +
    # Line & points per model
    geom_line(linewidth = 1) +
    geom_point(size = 2) +
    # Highlight significant points per model
    geom_point(
      data = es_both %>% filter(p.value < 0.05),
      size = 3, shape = 21, colour = "black"
    ) +
    scale_x_continuous(breaks = scales::pretty_breaks(n = 12)) +
    scale_colour_manual(values = c("Unadjusted" = "steelblue", "Adjusted" = "coral4")) +
    scale_fill_manual(values   = c("Unadjusted" = "steelblue", "Adjusted" = "coral4")) +
    scale_linetype_manual(values= c("Unadjusted" = "solid",     "Adjusted" = "dotdash")) +
    labs(
      title    = title_str,
      subtitle = "Sun & Abraham event study | FEGLM | 95% CI | Clustered SE (dept_id)",
      x        = "Months relative to adoption",
      y        = "Coefficient estimate",
      colour   = NULL, fill = NULL, linetype = NULL,
      caption  = "Black-rimmed points: p < 0.05. Vertical dotted line: adoption."
    ) }
  
  p_es_reads_pois <- plot_es_compare_pois(es_reads_pois, es_areads_pois, "Event Study – Readmission rate")
  p_es_u_reads_pois <- plot_es_compare_pois(es_u_reads_pois, es_u_areads_pois, "Event Study – Unplanned Readmission rate")
  #p_es_ed_pois    <- plot_es_compare_pois(es_ed_pois, es_aed_pois,   "Event Study – ED Visit rate")
  ggsave("R/group analysis/outputs/event_study_reads_pois.png",    p_es_reads_pois, width=10, height=5, dpi=150)
  ggsave("R/group analysis/outputs/event_study_ureads_pois.png",    p_es_u_reads_pois, width=10, height=5, dpi=150)
  #ggsave("R/group analysis/outputs/event_study_edvisits_pois.png", p_es_ed_pois,    width=10, height=5, dpi=150)
  ggsave("R/group analysis/outputs/event_study_combined_pois.png",
         p_es_reads_pois / p_es_u_reads_pois, width=10, height=10, dpi=150)


# =============================================================================
# 6. Export results to Excel
# =============================================================================

# Overall ATT
  extract_att <- function(agg_mod, outcome) {
    data.frame(
      outcome  = outcome,
      ATT      = agg_mod["ATT", "Estimate"],
      SE       = agg_mod["ATT", "Std. Error"],
      ci_lo    = agg_mod["ATT", "Estimate"] - 1.96 * agg_mod["ATT", "Std. Error"],
      ci_hi    = agg_mod["ATT", "Estimate"] + 1.96 * agg_mod["ATT", "Std. Error"],
      p_value  = agg_mod["ATT", "Pr(>|t|)"]
    )
  }
  # att_summary <- bind_rows(
  #   extract_att(att_sa_reads,    "pct_reads"),
  #   extract_att(att_sa_u_reads,    "pct_unplanned_reads"),
  #   extract_att(att_sa_edvisits, "pct_edvisits")
  # )

# Cohort ATTs
  extract_cohort_att <- function(agg_mod, outcome) {
    data.frame(
      outcome      = outcome,
      cohort       = rownames(agg_mod),
      ATT          = agg_mod[, "Estimate"],
      SE           = agg_mod[, "Std. Error"],
      ci_lo        = agg_mod[, "Estimate"] - 1.96 * agg_mod[, "Std. Error"],
      ci_hi        = agg_mod[, "Estimate"] + 1.96 * agg_mod[, "Std. Error"],
      p_value      = agg_mod[, "Pr(>|t|)"],
      row.names    = NULL
    )
  }

  # cohort_atts <- bind_rows(
  #   extract_cohort_att(att_cohort_reads,    "pct_reads"),
  #   extract_cohort_att(att_cohort_u_reads,    "pct_unplanned_reads"),
  #   extract_cohort_att(att_cohort_edvisits, "pct_edvisits")
  # )

# Wald test summary
  wald_summary <- data.frame(
    outcome = c("pct_reads", "pct_unplanned_reads"  , "pct_edvisits"),
    F_stat  = c(pt_reads$wald$stat, pt_u_reads$wald$stat, pt_edvisits$wald$stat),
    p_value = c(pt_reads$wald$stat,  pt_u_reads$wald$stat,   pt_edvisits$wald$p)
  )

extract_es_coefs <- function(mod, outcome) {
  b    <- mod$coefficients
  se   <- sqrt(diag(mod$cov.scaled))
  keep <- grepl("^month_index::", names(b))
  b    <- b[keep]
  se   <- se[keep]
  
  periods <- suppressWarnings(as.integer(sub("^month_index::([-0-9]+):.*", "\\1", names(b))))
  
  data.frame(
    outcome  = outcome,
    period   = periods,
    estimate = unname(b),
    se       = unname(se),
    stringsAsFactors = FALSE
  ) %>%
    filter(!is.na(period)) %>%
    group_by(outcome, period) %>%
    summarise(
      estimate = mean(estimate),
      se       = sqrt(mean(se^2)),
      .groups  = "drop"
    ) %>%
    mutate(
      ci_lo = estimate - 1.96 * se,
      ci_hi = estimate + 1.96 * se
    ) %>%
    arrange(period)
}

es_coefs <- bind_rows(
  extract_es_coefs(sa_reads,    "pct_reads"),
  extract_es_coefs(sa_u_reads,    "pct_unplanned_reads"),
  extract_es_coefs(sa_edvisits, "pct_edvisits")
)

write_xlsx(
  list(
#ATT_Overall"         = att_summary,
 #   "ATT_By_Cohort"       = cohort_atts,
    "PreTrend_Wald"       = wald_summary,
    "PreTrend_Individual" = lead_table,
    "EventStudy_Coefs"    = es_coefs
  ),
  path = "R/group analysis/outputs/sunab_results.xlsx"
)

# =============================================================================
# 7. Testing heterogeneity in Sun & Abraham models
# =============================================================================

#---Heterogeneity by high outreach

    table(df_sa$g)
    table(df_sa$trt_outreach_group)
    
    df_sa <- df_sa %>%
      mutate(
        g_int_outreach = ifelse(
          trt_outreach_group == "Treatment, No High Outreach",
          g,
          ifelse(treatment == "Treatment", NA, 0)))
    table(df_sa$g_int_outreach)

#---Heterogeneity by CCI median
    df_sa <- df_sa %>%
      mutate(
        g_int_cci = ifelse(
          scenario_label == "Scenario 1: Treatment, CCI avg above median",
          g,
          ifelse(treatment == "Treatment", NA, 0)))  # removing treatment groups from analysis
    table(df_sa$g_int_cci)

    run_sunab_het <- function(outcome, g, data = df_sa) {
      fml <- as.formula(
        paste0(outcome, " ~ sunab(", g, ",month_index, ref.p = -1)  | dept_id +  month_fe  ") #change to CCI model
      )
      feols(fml, data = data, cluster = ~dept_id)
    }
    
      hetout_sa_u_reads <- run_sunab_het("pct_unplanned_reads", "g_int_outreach")
      hetcci_sa_u_reads <- run_sunab_het("pct_unplanned_reads", "g_int_cci")
        #Outreach interaction
        hetout_att_u_reads    <- summary(hetout_sa_u_reads,    agg = "att")
          etable(hetout_att_u_reads, coefstat = "confint", ci = 0.90,
                 signif.code = c("***" = 0.01, "**" = 0.05, "*" = 0.10),
                 notes = "Confidence intervals shown are 90%.")
        hetout_att_cohort_u_reads    <- summary(hetout_sa_u_reads,    agg = "cohort")
          etable(hetout_att_cohort_u_reads, coefstat = "confint", ci = 0.90,
                 signif.code = c("***" = 0.01, "**" = 0.05, "*" = 0.10),
                 notes = "Confidence intervals shown are 90%.")
        #CCI interaction  
        hetcci_att_u_reads    <- summary(hetcci_sa_u_reads,    agg = "att")
        etable(hetcci_att_u_reads, coefstat = "confint", ci = 0.90,
               signif.code = c("***" = 0.01, "**" = 0.05, "*" = 0.10),
               notes = "Confidence intervals shown are 90%.")
        hetcci_att_cohort_u_reads    <- summary(hetcci_sa_u_reads,    agg = "cohort")
        etable(hetcci_att_cohort_u_reads, coefstat = "confint", ci = 0.90,
               signif.code = c("***" = 0.01, "**" = 0.05, "*" = 0.10),
               notes = "Confidence intervals shown are 90%.")


#Extract coefficients from Sunab models first
  extract_es <- function(mod) {
    broom::tidy(mod, conf.int = TRUE) |>
      dplyr::filter(stringr::str_detect(term, "^month_index::\\s*-?\\d+$")) |>
      dplyr::mutate(rel_time = as.numeric(stringr::str_extract(term, "-?\\d+$"))) |>
      dplyr::select(rel_time, estimate, std.error, statistic, p.value, conf.low, conf.high) |>
      dplyr::arrange(rel_time)
  }


  #het_es_u_reads <- extract_es(het_sa_u_reads) %>% mutate(outcome = "pct_uplanned_reads")

#-----Plot coefficients for event study visual  -- unadjusted vs het model
    plot_es_compare <- function(es_unadj, es_het, title_str) {
      # Add a model label to each input and combine
      es_unadj2 <- es_unadj %>% mutate(model = "Unadjusted")
      es_het2   <- es_het   %>% mutate(model = "Heterogeneity")
      
      es_both <- bind_rows(es_unadj2, es_het2)
      
      ggplot(es_both, aes(rel_time, estimate, colour = model, fill = model, linetype = model)) +
        geom_hline(yintercept = 0, linetype = "dashed", colour = "grey50") +
        geom_vline(xintercept = 0, linetype = "dotted", colour = "red") +
        # CI ribbon per model
        geom_ribbon(aes(ymin = conf.low, ymax = conf.high), alpha = 0.15, colour = NA) +
        # Line & points per model
        geom_line(linewidth = 1) +
        geom_point(size = 2) +
        # Highlight significant points per model
        geom_point(
          data = es_both %>% filter(p.value < 0.05),
          size = 3, shape = 21, colour = "black"
        ) +
        scale_x_continuous(breaks = scales::pretty_breaks(n = 12)) +
        scale_colour_manual(values = c("Unadjusted" = "steelblue", "Heterogeneity" = "coral4")) +
        scale_fill_manual(values   = c("Unadjusted" = "steelblue", "Heterogeneity" = "coral4")) +
        scale_linetype_manual(values= c("Unadjusted" = "solid",     "Heterogeneity" = "dotdash")) +
        labs(
          title    = title_str,
          subtitle = "Sun & Abraham event study | FEOLS | 95% CI | Clustered SE (dept_id)",
          x        = "Months relative to adoption",
          y        = "Coefficient estimate",
          colour   = NULL, fill = NULL, linetype = NULL,
          caption  = "Black-rimmed points: p < 0.05. Vertical dotted line: adoption."
        ) }


  p_es_u_reads <- plot_es_compare(es_u_reads, es_u_reads, "Event Study – % Unplanned Readmissions")


ggsave("R/group analysis/outputs/het_event_study_pct_ureads.png",    p_es_u_reads, width=10, height=5, dpi=150)


##------- Comparing models ---------------  


#Pretty table
modelsummary(list(att_sa_u_reads, att_cohort_u_reads, hetout_att_u_reads , hetout_att_cohort_u_reads, hetcci_att_u_reads , hetcci_att_cohort_u_reads), 
             conf_level = 0.90, statistic = "conf.int", 
             stars = c("*" = .10, "**" = .05, "***" = .01), 
             fmt = 2)


########
rm(list = setdiff(ls(), c("df_test", "df_analysis", "df_main", 
                          "df_matches", "df_pdc", "df_did",
                          "df_its", "df_sa", "df_cs",
                          "df_its_ips_out", "df_test_ips_out",
                          "df_its_ips_dept", "df_test_ips_dept")))
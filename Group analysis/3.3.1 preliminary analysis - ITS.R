# =============================================================================
# Interrupted Time Series (ITS) Analysis
# Single intervention point: m24 = 13 (Jan 2025)
# Outcomes: readmissions, unplanned readmissions, edvisits
# Estimator: CITS via fixest
# =============================================================================
required_pkgs <- c("fixest", "ggplot2", "dplyr", "writexl", "patchwork")

to_install <- required_pkgs[!sapply(required_pkgs, requireNamespace, quietly = TRUE)]
if (length(to_install)) install.packages(to_install)

library(fixest)
library(ggplot2)
library(writexl)
library(patchwork)
library(dplyr)
# -----------------------------------------------------------------------------
# Summary
# -----------------------------------------------------------------------------

# Design: Single-cutpoint ITS, dept x month panel (Jan 2024-Dec 2025)
# Excluded: Colorectal Surgery, Urology Surgery, Vascular Surgery
# Intervention: m24 = 13 (Jan 2025) for all treated departments\n
# 
# MODEL SPEC
#   feols(outcome ~ ~ treated * time_pre + treated * post + treated * time_since
#         | dept_id + month_fe)
#   treated:post      = differential post-period level shift for treated vs control
#   treated:time_since = differential slope change post-intervention for treated vs control
#   post = immediate level change for control group
#   time_pre = slope for control group in pre-period
#   treated:time_pre = differential pre-period slope for treated vs control
#      ->non zero values suggest parallel trends violation
#   FE: dept_id +  month_fe (seasonality)
#   SE: clustered at dept_id\n
# 
# OUTPUTS
#   Plots : its_pct_reads.png, its_pct_edvisits.png, its_overview.png
#           its_poisson_pct_reads.png, its_poisson_pct_edvisits.png, its_poisson_overview.png
#   Tables: its_results.xlsx (5 sheets)
#     - ITS_Coefficients : OLS level shift and slope change estimates
#     - Poisson_IRR      : Incidence rate ratios from Poisson model
#     - Overdispersion   : Pearson overdispersion ratios (>2 = use negbin)
#     - PreTrend_Wald    : Joint F-test on pre-period parallel trends
#     - Obs_vs_Fitted    : Observed and fitted means by group and period\n
# 
# CAVEATS
#   - Parallel trends assumed pre-intervention; check Wald p-value.
#   - Level shift and slope change are jointly identified; interpret together.
#   - Many zero outcomes in sample; consider fractional logit as robustness check.
#   - dept_id and calendar month FE absorb unit and time effects; interactions
#     identify treatment effects of within-unit variation.


# -----------------------------------------------------------------------------
# 1. Data preparation
# -----------------------------------------------------------------------------

  #Pull df_its from data setup script
  table(df_its$treated)
  table(df_its$m24) 
  #Sample after exclusions:
  print(table(df_its$treatment))
  #Departments retained:
  print(sort(unique(df_its$outreach_name)))

# -----------------------------------------------------------------------------
# 2. ITS model
# -----------------------------------------------------------------------------
  # Level shift:  treated x post        => immediate jump at intervention
  # Slope change: treated x post x time  => differential post-intervention trend
  # Pre-trend slope differences: treated x time 
  # control group level change: treated
  
  run_its <- function(outcome, data = df_its) {
    fml <- as.formula(paste0(
      outcome,
      " ~ treated * time_pre + treated * post + treated * time_since ",
      " | dept_id +  month_fe"
    ))
    feols(fml, data = data, cluster = ~dept_id)
  }
  
  run_its_adj <- function(outcome, data = df_its) {
    fml <- as.formula(paste0(
      outcome,
      " ~ treated * time_pre + treated * post + treated * time_since ", 
      " + avg_cci + avg_age + pct_female",
      " | dept_id +  month_fe"
      
    ))
    feols(fml, data = data, cluster = ~dept_id)
}

#========== UNADJUSTED ITS MODEL: READMISSIONS AND ED VISITS ==========
its_reads <- run_its("pct_reads")
summary(its_reads)

its_u_reads <- run_its("pct_unplanned_reads")
summary(its_u_reads)

    etable(
      its_u_reads, coefstat = "confint", ci = 0.90,
      signif.code = c("***" = 0.01, "**" = 0.05, "*" = 0.10),
      notes = "Confidence intervals shown are 90%."
    )

its_edvisits <- run_its("pct_edvisits")
summary(its_edvisits)

#========== ADJUSTED ITS MODEL: READMISSIONS AND ED VISITS ==========
aits_reads <- run_its_adj("pct_reads")
summary(its_reads)

aits_u_reads <- run_its_adj("pct_unplanned_reads")
summary(its_u_reads)


    etable(
      aits_u_reads, coefstat = "confint", ci = 0.90,
      signif.code = c("***" = 0.01, "**" = 0.05, "*" = 0.10),
      notes = "Confidence intervals shown are 90%."
    )
aits_edvisits <- run_its_adj("pct_edvisits")
summary(its_edvisits)

# -----------------------------------------------------------------------------
# 3. Poisson model with offset (count ~ n)
# -----------------------------------------------------------------------------
# Uses raw numerator (reads, edvisits) and denominator (n = total discharges)
# Coefficients are on the log scale; exponentiate for incidence rate ratios (IRR)

run_poisson <- function(outcome_count, denominator = "n", data = df_its) {
  data$log_denom <- log(data[[denominator]])
  fml <- as.formula(paste0(
    outcome_count,
    " ~ treated * time_pre + treated * post + treated * time_since + offset(log_denom) | dept_id +  month_fe"
  ))
  feglm(fml, data = data, cluster = ~dept_id, family = poisson())
}
run_poisson_adj <- function(outcome_count, denominator = "n", data = df_its) {
  data$log_denom <- log(data[[denominator]])
  fml <- as.formula(paste0(
    outcome_count,
    "~ treated * time_pre + treated * post + treated * time_since + offset(log_denom) + avg_cci + avg_age + pct_female", 
    "| dept_id +  month_fe"
  ))
  feglm(fml, data = data, cluster = ~dept_id, family = poisson())
}

#========== Unadjusted POISSON MODEL: reads (offset = n) ==========
pois_reads <- run_poisson("reads")
summary(pois_reads)
pois_u_reads <- run_poisson("unplanned_reads")
summary(pois_u_reads)
      etable(
        pois_u_reads, coefstat = "confint", ci = 0.90,
        signif.code = c("***" = 0.01, "**" = 0.05, "*" = 0.10),
        notes = "Confidence intervals shown are 90%."
      )
pois_edvisits <- run_poisson("edvisits")
summary(pois_edvisits)

#========== Adjusted POISSON MODEL: edvisits (offset = n) ==========
apois_reads <- run_poisson_adj("reads")
summary(apois_reads)
apois_u_reads <- run_poisson_adj("unplanned_reads")
summary(apois_u_reads)
apois_edvisits <- run_poisson_adj("edvisits")
summary(apois_edvisits)

# -- Overdispersion check -----------------------------------------------------
# Ratio should be close to 1; values >> 1 suggest overdispersion -> use negbin
overdisp <- function(mod, label) {
  ratio <- sum(residuals(mod, type = "pearson")^2) / mod$nobs
  cat(sprintf("Overdispersion ratio (%s): %.3f %s\n", label, ratio,
              ifelse(ratio > 2, "  WARNING: consider negative binomial", "")))
  ratio
}

od_reads    <- overdisp(pois_reads,    "reads")
od_u_reads    <- overdisp(pois_u_reads,    "unplanned_reads")
od_edvisits <- overdisp(pois_edvisits, "edvisits")

# -- IRR table (exponentiated coefficients) -----------------------------------
extract_irr <- function(mod, outcome) {
  ct <- coeftable(mod)
  cf <- confint(mod, level = 0.95)
  keep <- grepl("treated:post|treated:time_since", rownames(ct))
  data.frame(
    outcome  = outcome,
    term     = rownames(ct)[keep],
    IRR      = exp(ct[keep, "Estimate"]),
    ci_lo    = exp(cf[keep, 1]),
    ci_hi    = exp(cf[keep, 2]),
    p_value  = ct[keep, "Pr(>|z|)"],
    row.names = NULL
  )
}

irr_summary <- bind_rows(
  extract_irr(pois_reads,    "reads"),
  extract_irr(pois_u_reads,    "unplanned_reads"),
  extract_irr(pois_edvisits, "edvisits"),
  extract_irr(apois_reads,    "reads"),
  extract_irr(apois_u_reads,    "unplanned_reads"),
  extract_irr(apois_edvisits, "edvisits")
)

#Poisson IRR summary:
print(irr_summary)



# -----------------------------------------------------------------------------
# 4. Pre-trend test (joint F-test on pre-period event-study leads)
# -----------------------------------------------------------------------------
# Fit an event-study version in the pre-period to test parallel trends

run_pretrend <- function(outcome, data = df_its) {
  fml <- as.formula(paste0(
    outcome,
    " ~ i(m24, treated, ref = 6) | dept_id "
  ))
  feols(fml, data = data[data$m24 < 13, ], cluster = ~dept_id)
}

run_pretrend_adj <- function(outcome, data = df_its) {
  fml <- as.formula(paste0(
    outcome,
    " ~ i(m24, treated, ref = 6) | dept_id +  month_fe"
  ))
  feols(fml, data = data[data$m24 < 13, ], cluster = ~dept_id)
}

#========== PRE-TREND TEST: readmissions, unplanned readmissions, ed visits ==========
#Unadjusted model
pre_reads <- run_pretrend("pct_reads")
summary(pre_reads)

pre_u_reads <- run_pretrend("pct_unplanned_reads")
summary(pre_u_reads)

pre_edvisits <- run_pretrend("pct_edvisits")
summary(pre_edvisits)

#*Adjusted model (conditional trends assumption)
pre_areads <- run_pretrend_adj("pct_reads")
summary(pre_areads)

pre_u_areads <- run_pretrend_adj("pct_unplanned_reads")
summary(pre_u_areads)

pre_aedvisits <- run_pretrend_adj("pct_edvisits")
summary(pre_aedvisits)


# Joint Wald F-test on all pre-period interaction terms
wald_reads        <- wald(pre_reads)

wald_u_reads      <- wald(pre_u_reads)

wald_edvisits     <- wald(pre_edvisits)

wald_areads        <- wald(pre_areads)

wald_u_areads      <- wald(pre_u_areads)

wald_aedvisits     <- wald(pre_aedvisits)

# -----------------------------------------------------------------------------
# 4. ITS plots: observed means + fitted regression lines
# -----------------------------------------------------------------------------

##--------------- FEOLS First------------------------------------

# Compute observed group means per period
obs_means <- df_its %>%
  group_by(m24, treatment) %>%
  summarise(
    pct_reads    = mean(pct_reads,    na.rm = TRUE),
    pct_u_reads    = mean(pct_unplanned_reads,    na.rm = TRUE),
    pct_edvisits = mean(pct_edvisits, na.rm = TRUE),
    .groups = "drop"
  )

# Compute fitted model means
  #Fitted values from model-UNADJUSTED
  df_its$fitted_reads <- predict(its_reads, newdata = df_its)
  df_its$fitted_u_reads <- predict(its_u_reads, newdata = df_its)
  df_its$fitted_edvisits <- predict(its_edvisits, newdata = df_its)
  
  #Fitted values from model-ADJUSTED
  df_its$afitted_reads <- predict(aits_reads, newdata = df_its)
  df_its$afitted_u_reads <- predict(aits_u_reads, newdata = df_its)
  df_its$afitted_edvisits <- predict(aits_edvisits, newdata = df_its)


  fitted_means <- df_its %>%
    group_by(m24, treatment) %>%
    summarise(
      fit_reads    = mean(fitted_reads,    na.rm = TRUE),
      fit_u_reads    = mean(fitted_u_reads,    na.rm = TRUE),
      fit_edvisits = mean(fitted_edvisits, na.rm = TRUE),
      afit_reads    = mean(afitted_reads,    na.rm = TRUE),
      afit_u_reads    = mean(afitted_u_reads,    na.rm = TRUE),
      afit_edvisits = mean(afitted_edvisits, na.rm = TRUE),
      .groups = "drop"
    )



# #plot fitted vs observed means for each outcome -function for unadjusted and adjusted
# plot_its <- function(obs, fit, obs_var, fit_var, afit_var, title, ylab) {
#   obs_df <- obs %>% dplyr::select(m24, treatment, value = all_of(obs_var))
#   #fit_df <- fit %>% dplyr::select(m24, treatment, value = all_of(c(fit_var, afit_var)))
#   
#   fit_df <- fit %>%
#     dplyr::select(m24, treatment, any_of(c(fit_var, afit_var))) %>%
#     pivot_longer(cols = any_of(c(fit_var, afit_var)),
#                  names_to = "model",
#                  values_to = "value")
#   
#   series_present <- unique(fit_df$model) # Determine all series that are present
#   
#   # Build linetype map (values) for only those present
#   lt_full <- c( "dotdash", "solid")
#   names(lt_full) <- c(fit_var, afit_var)
#   lt_map <- lt_full[series_present]
#   
#   # Labels for lines
#   lab_full <- c("Unadjusted", "Adjusted")
#   names(lab_full) <- c(fit_var, afit_var)
#   lab_map <- lab_full[series_present]
#   
#   ggplot() +
#     geom_vline(
#       xintercept = 12.5, linetype = "dashed",
#       colour = "grey40", linewidth = 0.8
#     ) +
#     geom_point(
#       data = obs_df,
#       aes(x = m24, y = value, colour = treatment),
#       alpha = 0.5, size = 1.8
#     ) +
#     # Lines for one or both fitted series (works if only one exists too)
#     geom_line(
#       data = fit_df,
#       aes(x = m24, y = value, colour = treatment, linetype = model),
#       linewidth = 1
#     ) +
#     scale_colour_manual(values = c("Treatment" = "steelblue",
#                                    "Control"   = "coral")) +
#     #scale_y_continuous(limits=c(0,25))+
#     scale_linetype_manual(values = lt_map, labels = lab_map) +
#     labs(
#       title   = title,
#       x       = "Month (m24)",
#       y       = ylab,
#       colour  = NULL,
#       caption = "Points: observed group means  |  Lines: model fitted values  |  Dashed: intervention (m24=13)"
#     ) +
#     theme_minimal(base_size = 12) +
#     theme(
#       plot.title     = element_text(face = "bold"),
#       legend.position = "top")
# }

#plot fitted vs observed means for each outcome - function for just unadjusted
plot_its <- function(obs, fit, obs_var, fit_var, title, ylab) {
  obs_df <- obs %>% dplyr::select(m24, treatment, value = all_of(obs_var))
  #fit_df <- fit %>% dplyr::select(m24, treatment, value = all_of(c(fit_var, afit_var)))

  fit_df <- fit %>%
    dplyr::select(m24, treatment, any_of(c(fit_var))) %>%
    pivot_longer(cols = any_of(c(fit_var)),
                 names_to = "model",
                 values_to = "value")

  series_present <- unique(fit_df$model) # Determine all series that are present

  # Build linetype map (values) for only those present
  lt_full <- c( "solid")
  names(lt_full) <- c(fit_var)
  lt_map <- lt_full[series_present]

  # Labels for lines
  lab_full <- c("Unadjusted")
  names(lab_full) <- c(fit_var)
  lab_map <- lab_full[series_present]

  ggplot() +
    geom_vline(
      xintercept = 12.5, linetype = "dashed",
      colour = "grey40", linewidth = 0.8) +
    geom_point(
      data = obs_df,
      aes(x = m24, y = value, colour = treatment),
      alpha = 0.5, size = 1.8) +
    # Lines for one or both fitted series (works if only one exists too)
    geom_line(
      data = fit_df,
      aes(x = m24, y = value, colour = treatment),
      linewidth = 1 ) +
    scale_y_continuous(limits = c(6,14)) +
    scale_x_continuous(breaks = seq(0,24, by=24),
                       labels = c("Jan '24", "Dec '25")) +
    scale_colour_manual(values = c("Treatment" = "steelblue",
                                   "Control"   = "coral")) +
    #scale_y_continuous(limits=c(0,25))+
    scale_linetype_manual(values = lt_map, labels = lab_map) +
    labs(
      title   = title,
      x       = "Month",
      y       = ylab,
      colour  = NULL,
     # caption = "Points: observed group means  |  Lines: model fitted values  |  Dashed: intervention (Jan '25)"
    ) +
    theme_minimal(base_size = 16) +
    theme(
      plot.title     = element_text(face = "bold"),
      legend.position = "top")
}


p_u_reads <- plot_its(obs_means, fitted_means,
                    "pct_u_reads", "fit_u_reads", 
                    "% Unplanned Readmissions",
                    "% Unplanned Readmissions")
p_edvisits <- plot_its(obs_means, fitted_means,
                       "pct_edvisits", "fit_edvisits", 
                       "% ED Visits",
                       "% ED Visits")

ggsave("R/group analysis/outputs/its_pct_u_reads.png",    p_u_reads,    width = 9, height = 5, dpi = 300)
ggsave("R/group analysis/outputs/its_pct_edvisits.png", p_edvisits, width = 9, height = 5, dpi = 300)

# Combined overview plot
p_overview <-  p_u_reads | p_edvisits +
  plot_annotation(
    title    = "Interrupted Time Series Results",
    subtitle = "Single intervention time  |  dept_id + calendar-period FEs  | FEOLS, SE clustered at dept_id",
    theme    = theme(plot.title = element_text(face = "bold", size = 16)),
  )

ggsave("R/group analysis/outputs/its_overview.png", p_overview, width = 18, height = 10, dpi = 300)

##---------------------- Poisson second ----------------------------------------------

  # Poisson fitted values (on response/rate scale) -UNADJUSTED
  df_its$fitted_pois_reads    <- fitted(pois_reads)    / df_its$n
  df_its$fitted_pois_u_reads    <- fitted(pois_u_reads)    / df_its$n
  df_its$fitted_pois_edvisits <- fitted(pois_edvisits) / df_its$n
  
  #Fitted values from model-ADJUSTED
      # add offset used in estimation exists in newdata
      df_its$log_denom <- log(df_its$n)
  
  df_its$afitted_pois_reads    <- predict(apois_reads, newdata = df_its, type = "response")/ df_its$n # response = expected counts for Poisson)#
  df_its$afitted_pois_u_reads    <- predict(apois_u_reads, newdata = df_its, type = "response" )/ df_its$n
  df_its$afitted_pois_edvisits <- predict(apois_edvisits, newdata = df_its, type = "response" ) / df_its$n
  

  
    fitted_pois_means <- df_its %>%
    group_by(m24, treatment) %>%
    summarise(
      fit_pois_reads    = mean(fitted_pois_reads,    na.rm = TRUE),
      fit_pois_u_reads    = mean(fitted_pois_u_reads,    na.rm = TRUE),
      fit_pois_edvisits = mean(fitted_pois_edvisits, na.rm = TRUE),
      afit_pois_reads    = mean(afitted_pois_reads,    na.rm = TRUE),
      afit_pois_u_reads    = mean(afitted_pois_u_reads,    na.rm = TRUE),
      afit_pois_edvisits = mean(afitted_pois_edvisits, na.rm = TRUE),
      .groups = "drop"
    )

#plot fitted vs observed means for each outcome
  plot_its_fitted <- function(fit, fit_var, afit_var, title, ylab) {
 # fit_df <- fit %>% dplyr:: select(m24, treatment, value = all_of(fit_var))
  
  fit_df <- fit %>%
    dplyr::select(m24, treatment, any_of(c(fit_var, afit_var))) %>%
    pivot_longer(cols = any_of(c(fit_var, afit_var)),
                 names_to = "model",
                 values_to = "value")
  # Determine all series that are present
  series_present <- unique(fit_df$model)
  
  # Build linetype map (values) for only those present
  lt_full <- c( "dotdash", "solid")
  names(lt_full) <- c(fit_var, afit_var)
  lt_map <- lt_full[series_present]
  
  # Labels 
  lab_full <- c("Unadjusted", "Adjusted")
  names(lab_full) <- c(fit_var, afit_var)
  lab_map <- lab_full[series_present]
  
  ggplot() +
    geom_vline(
      xintercept = 12.5, linetype = "dashed",
      colour = "grey40", linewidth = 0.8
    ) +
    
    # Lines for one or both fitted series (works if only one exists too)
    geom_line(
      data = fit_df,
      aes(x = m24, y = value, colour = treatment, linetype = model),
      linewidth = 1
    ) +
    scale_colour_manual(values = c("Treatment" = "steelblue",
                                   "Control"   = "coral")) +
    scale_linetype_manual(values = lt_map, labels = lab_map) +
   # scale_y_continuous(limits=c(0,0.25))+
    labs(
      title   = title,
      x       = "Month (m24)",
      y       = ylab,
      colour  = NULL,
      caption = "Points: observed group means  |  Lines: model fitted values  |  Dashed: intervention (m24=13)"
    ) +
    theme_minimal(base_size = 12) +
    theme(
      plot.title     = element_text(face = "bold"),
      legend.position = "top")
}

p_pois_reads <- plot_its_fitted(fitted_pois_means, "fit_pois_reads", "afit_pois_reads",
                                "Poisson ITS Fitted Trends: % Readmissions",
                                "% Readmissions")
p_pois_u_reads <- plot_its_fitted(fitted_pois_means, "fit_pois_u_reads", "afit_pois_u_reads",
                                "Poisson ITS Fitted Trends: % Unplanned Readmissions",
                                "% Uplanned Readmissions")
p_pois_edvisits <- plot_its_fitted(fitted_pois_means, "fit_pois_edvisits", "afit_pois_edvisits",
                                   "Poisson ITS Fitted Trends: % ED Visits",
                                   "% ED Visits")

ggsave("R/group analysis/outputs/its_poisson_pct_reads.png",    p_pois_reads,    width = 9, height = 5, dpi = 300)
ggsave("R/group analysis/outputs/its_poisson_pct_u_reads.png",    p_pois_reads,    width = 9, height = 5, dpi = 300)
ggsave("R/group analysis/outputs/its_poisson_pct_edvisits.png", p_pois_edvisits, width = 9, height = 5, dpi = 300)

# Combined overview plot
p_pois_overview <- p_pois_reads / p_pois_u_reads / p_pois_edvisits +
  plot_annotation(
    title    = "Poisson ITS Fitted Trend Results",
    subtitle = "Fitted values on rate scale (counts / n)  |  dept_id + calendar-period FEs  | FEGLM, SE clustered at dept_id",
    theme    = theme(plot.title = element_text(face = "bold", size = 14))
  )

ggsave("R/group analysis/outputs/its_poisson_overview.png", p_pois_overview, width = 10, height = 10, dpi = 300)

# -----------------------------------------------------------------------------
# 5. Export results to Excel
# -----------------------------------------------------------------------------

  # Extract key ITS coefficients from model
  extract_its_coefs <- function(mod, outcome) {
    ct <- coeftable(mod)
    cf <- confint(mod, level = 0.95)
    # Keep only the two treatment interaction terms
    keep <- grepl("treated:post|treated:time_since", rownames(ct))
    data.frame(
      outcome    = outcome,
      term       = rownames(ct)[keep],
      estimate   = ct[keep, "Estimate"],
      se         = ct[keep, "Std. Error"],
      ci_lo      = cf[keep, 1],
      ci_hi      = cf[keep, 2],
      p_value    = ct[keep, "Pr(>|t|)"],
      row.names  = NULL
    )
  }
  
  coef_summary <- bind_rows(
    #unadjusted
    extract_its_coefs(its_reads,    "pct_reads"),
    extract_its_coefs(aits_reads,    "pct_reads"),
    extract_its_coefs(its_u_reads,    "pct_unplanned_reads"),
    #adjusted
    extract_its_coefs(aits_u_reads,    "pct_unplanned_reads"),
    extract_its_coefs(its_edvisits, "pct_edvisits"),
    extract_its_coefs(aits_edvisits, "pct_edvisits")
    
  )

  # Wald test summary
  wald_summary <- data.frame(
    outcome = c("pct_reads", "pct_unplanned_reads", "pct_edvisits"),
    F_stat  = c(if (!is.null(wald_reads))    wald_reads$stat    else NA,
                if (!is.null(wald_u_reads))    wald_u_reads$stat    else NA,
                if (!is.null(wald_edvisits)) wald_edvisits$stat else NA,
                if (!is.null(wald_areads))    wald_areads$stat    else NA,
                if (!is.null(wald_u_areads))    wald_u_areads$stat    else NA,
                if (!is.null(wald_aedvisits)) wald_aedvisits$stat else NA),
    p_value = c(if (!is.null(wald_reads))    wald_reads$p       else NA,
                if (!is.null(wald_u_reads))    wald_u_reads$p       else NA,
                if (!is.null(wald_edvisits)) wald_edvisits$p    else NA,
                if (!is.null(wald_areads))    wald_areads$p       else NA,
                if (!is.null(wald_u_areads))    wald_u_areads$p       else NA,
                if (!is.null(wald_aedvisits)) wald_aedvisits$p    else NA)
  )

  # Observed vs fitted means for export
  obs_fit_export <- left_join(obs_means, fitted_means, by = c("m24", "treatment"))
  
  write_xlsx(
    list(
      "ITS_Coefficients" = coef_summary,
      "Poisson_IRR"      = irr_summary,
      "Overdispersion"   = data.frame(
                             outcome = c("reads", "unplanned_reads", "edvisits"),
                             overdisp_ratio = c(od_reads, od_u_reads, od_edvisits)
                           ),
      "PreTrend_Wald"    = wald_summary,
      "Obs_vs_Fitted"    = obs_fit_export
    ),
    path = "R/group analysis/outputs/its_results.xlsx"
  )


  rm(list = setdiff(ls(), c("df_test", "df_analysis", "df_main", 
                            "df_matches", "df_pdc", "df_did",
                            "df_its", "df_sa", "df_cs",
                            "df_its_ips_out", "df_test_ips_out",
                            "df_its_ips_dept", "df_test_ips_dept")))
  

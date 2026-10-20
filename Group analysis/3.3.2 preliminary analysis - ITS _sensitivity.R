# =============================================================================
# Interrupted Time Series (ITS) Analysis
# Single intervention point: m24 = 13 (Jan 2025)
# Outcomes: readmissions, unplanned readmissions, edvisits
# Estimator: TWFE via fixest
# =============================================================================
required_pkgs <- c("fixest", "ggplot2", "dplyr", "writexl", "patchwork")

to_install <- required_pkgs[!sapply(required_pkgs, requireNamespace, quietly = TRUE)]
if (length(to_install)) install.packages(to_install)

library(fixest)
library(ggplot2)
library(writexl)
library(patchwork)
library(dplyr)
library(modelsummary)
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
  avgoutreach <- df_its %>%
  group_by(department) %>%
  summarize(post_outreach = mean(pct_outreach[post == 1]), #, na.rm = TRUE),
            cci_avggroup = mean(avg_cci), #na.rm= TRUE),
            rr_avggroup = mean(dis_rrs)) %>%
                  
  ungroup()
  summary(avgoutreach$rr_avggroup)
  
  #develop some subgroups for risk / implementation related indicators
  df_its <- merge(df_its, avgoutreach, by = "department") %>%
    group_by(department) %>%
    mutate( high_rr = ifelse(is.na(rr_avggroup) | is.nan(rr_avggroup), 0, 
                             ifelse(rr_avggroup>10, 1, 0)),
            high_outreach = ifelse(is.na(post_outreach), 0,
                                  ifelse(post_outreach > 80, 1, 0)),
           cci_median = ifelse(cci_avggroup > 3.5, 1, 0)) %>%  # binary variable to identify groups that received high outreach at any point
    ungroup() %>%
    
    #create a new grouping to distinguish treatment and control groups by any high outreach in the post-period
    mutate(
      trt_outreach_group = case_when(
        treatment == "Treatment" & high_outreach == 1 ~ "Treatment, High Outreach",
        treatment == "Treatment" & high_outreach == 0 ~ "Treatment, No High Outreach",
        # treatment == "Control" & ever_high_outreach == 1 ~ "Control, High Outreach",  
        #     # from current data there is one dept (Neurosurgery Scranton) that had high outreach bc it was only 1 patient and they received outreach in a month
        #     # this department only has pre-data from one month so I'm going to leave this one out. 
        treatment == "Control" & high_outreach == 0 ~ "Control, No High Outreach"
      )
    )
   summary(df_its$cci_median)

   
  #Sample after exclusions:
  print(table(df_its$treatment))
  #Departments retained:
  print(sort(unique(df_its$outreach_name)))
  #Departments in subgroups
  table(df_its$department, df_its$trt_outreach_group)
  table(df_its$department, df_its$high_rr)
  
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
    feols(fml, data = data, cluster = ~outreach_name)
  }
  
  run_its_adj <- function(outcome, data = df_its) {
    fml <- as.formula(paste0(
      outcome,
      " ~ treated * time_pre + treated * post + treated * time_since ", 
      " + avg_cci + avg_age + pct_female",
      " | dept_id + m24 + month_fe"
      
    ))
    feols(fml, data = data, cluster = ~dept_id)
}

#========== UNADJUSTED ITS MODEL: READMISSIONS AND ED VISITS ==========
its_reads <- run_its("pct_reads")
summary(its_reads)

its_u_reads <- run_its("pct_unplanned_reads")
summary(its_u_reads)

its_edvisits <- run_its("pct_edvisits")
summary(its_edvisits)

#========== ADJUSTED ITS MODEL: READMISSIONS AND ED VISITS ==========
aits_reads <- run_its_adj("pct_reads")
summary(its_reads)

aits_u_reads <- run_its_adj("pct_unplanned_reads")
summary(its_u_reads)

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
pois_edvisits <- run_poisson("edvisits")
summary(pois_edvisits)

#========== Adjusted POISSON MODEL: edvisits (offset = n) ==========
apois_reads <- run_poisson_adj("reads")
summary(apois_reads)
apois_u_reads <- run_poisson_adj("unplanned_reads")
summary(apois_u_reads)
apois_edvisits <- run_poisson_adj("edvisits")
summary(apois_edvisits)

# -----------------------------------------------------------------------------
# 2. Mediator model
# -----------------------------------------------------------------------------
# Mediator model: pct_outreach as outcome
m_med <- feols(
  pct_outreach ~  treated * time_pre + treated * post + treated * time_since | dept_id + month_fe,
  data = df_its, 
  cluster = ~dept_id
)

summary(m_med)
df_its$Mhat <- predict(m_med)
df_its$M_resid <- df_its$pct_outreach - df_its$Mhat

m_direct <- feols(
  pct_outreach ~  treated * time_pre + treated * post + treated * time_since + M_resid | dept_id + month_fe,
  data = df_its, 
  cluster = ~dept_id
)

summary(m_direct)

beta_TM <- tidy(m_med) %>%
  filter(term != "(Intercept)")   # sunab terms

beta_MY <- tidy(m_direct) %>%
  filter(term == "M_resid") %>%
  pull(estimate)

# Indirect effect for each sunab horizon:
indirect_effects <- beta_TM %>%
  mutate(indirect = estimate * beta_MY)

indirect_effects

direct <- tidy(m_direct) %>%
  filter(grepl("treated:", term)) %>%
  mutate(type = "direct")

m_primary <- its_u_reads

summary(m_primary)

total <- tidy(m_primary) %>%
  filter(grepl("treated:", term)) %>%
  mutate(type = "total")

indirect <- indirect_effects %>%
  mutate(type = "indirect",
         term = beta_TM$term)

all_effects <- bind_rows(total, direct, indirect)

all_effects

# -----------------------------------------------------------------------------
# 2. Effect modification model
# ----------------------------------------------------------------------------
feols(
  pct_unplanned_reads ~ treated * post * pct_outreach +
    treated * time_pre +
    treated * time_since  |
    dept_id + month_fe,
  data = df_its,
  cluster = ~dept_id
)

df_its <- df_its %>% 
  mutate(trt_outreach_group = 
           case_when (
  trt_outreach_group == "Treatment, High Outreach" ~1,
  trt_outreach_group == "Treatment, No High Outreach" ~0,
  trt_outreach_group == "Control, No High Outreach"~ 0 )
  )
  
m<-feols(
  pct_unplanned_reads ~
    treated * time_pre +
    treated * post  +
    treated:post:trt_outreach_group +
    treated * time_since |  
    dept_id + month_fe, 
  data = df_its,
  cluster = ~dept_id
)
#----------Fitted values from model-ADJUSTED
  df_its$fitted_u_reads <- predict(m, newdata = df_its)
  
  # For outreach groups
  fitted_means <- df_its %>%
    mutate(scenario = case_when(trt_outreach_group == 1 ~ 1L, 
                                trt_outreach_group ==0 & treated ==1 ~ 2L, 
                                treated == 0  ~ 3L)) %>%
    group_by(m24, scenario) %>%
    summarise(
      fit_u_reads    = mean(fitted_u_reads,    na.rm = TRUE),
      .groups = "drop"
    )

    # Labels for diff scenarios
    fitted_means <- fitted_means %>%
      mutate(
        scenario_label = factor(
          scenario,
          levels = c(1, 2, 3),
          labels = c(
            "Scenario 1: Treatment, High Outreach",
            "Scenario 2: Treatment, No High Outreach",
            "Scenario 3: Control" )))
   #---------------------------------------------------
   # For CCI groups
    # fitted_means <- df_its %>%
    #   mutate(scenario = case_when(cci_median ==1 & treated == 1 ~ 1L, 
    #                               cci_median ==0 & treated ==1 ~ 2L,
    #                               treated == 0  ~ 3L)) %>%
    #   group_by(m24, scenario) %>%
    #   summarise(
    #     fit_u_reads    = mean(fitted_u_reads,    na.rm = TRUE),
    #     .groups = "drop"
    #   )
    # 
    # # Labels for diff scenarios
    # fitted_means <- fitted_means %>%
    #   mutate(
    #     scenario_label = factor(
    #       scenario,
    #       levels = c(1, 2, 3),
    #       labels = c(
    #         "Scenario 1: Treatment, CCI avg above median",
    #         "Scenario 2: Treatment, CCI avg below median",
    #         "Scenario 3: Control" )))
    

ggplot(fitted_means, aes(x = m24, y = fit_u_reads, color = scenario_label)) +
  geom_line(size = 1.3) +
  geom_point(size = 2) +
 #scale_y_continuous(limits= c(8,13)) +
  scale_color_manual(
    values = c(
      "Scenario 1: Treatment, High Outreach" = "#1b9e77",
      "Scenario 2: Treatment, No High Outreach" = "#d95f02",
      "Scenario 3: Control" = "#7570b3"
    )) +
  #   scale_color_manual(
  #     values = c(
  #     "Scenario 1: Treatment, CCI avg above median"= "#1b9e77",
  #     "Scenario 2: Treatment, CCI avg below median"= "#d95f02",
  #     "Scenario 3: Control"  = "#7570b3"
  # )) +
  labs(
    title = "Model-Adjusted Fitted pct_unplanned_reads Over Time",
    subtitle = "Averaged fitted values by scenario",
    x = "Time (m24)",
    y = "Adjusted fitted pct_unplanned_reads",
    color = "Scenario"
  ) +
  geom_vline(
    xintercept = 12.5, linetype = "dashed",
    colour = "grey40", linewidth = 0.8
  ) +
  theme_minimal(base_size = 14) +
  theme(
    plot.title = element_text(face = "bold"),
    legend.position = "top",
    legend.title = element_text(size = 12),
    legend.text = element_text(size = 11)
  )

# -----------------------------------------------------------------------------
# 2. Heterogenous effect model
# ----------------------------------------------------------------------------
# Interaction model for high outreach treatment group (OLS)
  het_model <-feols(
    pct_unplanned_reads ~
      treated * time_pre +
      treated * post  +
      treated:post:trt_outreach_group +
      treated * time_since |  
      dept_id + month_fe, 
    data = df_its,
    cluster = ~dept_id
  )


  etable(
    het_model, coefstat = "confint", ci = 0.90,
    signif.code = c("***" = 0.01, "**" = 0.05, "*" = 0.10),
    notes = "Confidence intervals shown are 90%."
  )

# # Interaction model for high outreach treatment group (Poisson)
#   run_poisson_het <- function(outcome_count, denominator = "n", data = df_its) {
#     data$log_denom <- log(data[[denominator]])
#     fml <- as.formula(paste0(
#       outcome_count,
#       " ~ treated * time_pre + treated * post + treated * time_since",
#       "treated:post:trt_outreach_group + offset(log_denom) | dept_id +  month_fe"
#     ))
#     feglm(fml, data = data, cluster = ~dept_id, family = poisson())
#   }
#   
#   het_pois_model <-run_poisson_het("unplanned_reads")
#     
#   etable(
#     het_pois_model, coefstat = "confint", ci = 0.90,
#     signif.code = c("***" = 0.01, "**" = 0.05, "*" = 0.10),
#     notes = "Confidence intervals shown are 90%."
#   )
  
# Interaction model for CCI (indicator for above median of 3.5) (OLS)
  het_model_cci <-feols(
    pct_unplanned_reads ~
      treated * time_pre +
      treated * post  +
      treated:post:cci_median +
      treated * time_since |  
      dept_id + month_fe, 
    data = df_its, #subset(df_its, cci_median > 0),
    cluster = ~dept_id
  )
  
  
  etable(
    het_model_cci, coefstat = "confint", ci = 0.90,
    signif.code = c("***" = 0.01, "**" = 0.05, "*" = 0.10),
    notes = "Confidence intervals shown are 90%."
  )
  
  
##------- Comparing models ---------------  
  
  
  #Pretty table
  modelsummary(list(its_u_reads, het_model), 
               conf_level = 0.90, statistic = "conf.int", 
               stars = c("*" = .10, "**" = .05, "***" = .01), 
               fmt = 2)
  
  ###
  rm(list = setdiff(ls(), c("df_test", "df_analysis", "df_main", "df_matches", "df_pdc", "df_did","df_its", "df_sa")))
  
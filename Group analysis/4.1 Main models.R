# =============================================================================
# Main models from each estimator
# Estimators: Sun & Abraham (fixest::sunab), Calloway & Sant'anna (did), CITS via fixest
# Outcomes:  unplanned readmissions
# =============================================================================

# =============================================================================
# 0. Run individual model scripts to ensure correct data frames are used
# =============================================================================
source("R/group analysis/3.2.1 preliminary analysis - Staggered DiD_S&A.R")
source("R/group analysis/3.2.2 preliminary analysis - Staggered DiD - C&S.R")
source("R/group analysis/3.3.1 preliminary analysis - ITS.R")


# =============================================================================
# 1. Sun & Abraham models
# =============================================================================
#----OLS - Unadjusted
sunab <- feols(pct_unplanned_reads ~ sunab(g, month_index, ref.p = -1) | dept_id +  month_fe, 
        data = df_sa, 
        cluster = ~dept_id)

sunab_group <- summary(sunab, agg = "cohort")
sunab_overall <- summary(sunab, agg = "att")

#------Event study plot 

    png("R/group analysis/outputs/sunab_iplot.png", width = 800, height = 600)
     sunab_iplot <- iplot(
       sunab,
       ref.line = -1,
       xlab = "Months Relative to Intervention",
       ylab = "Effect on % Unplanned Readmissions",
       main = "Event Study: Sun & Abraham Model"
     )
     dev.off()
   
     
     sunab_ed <- feols(pct_edvisits ~ sunab(g, month_index, ref.p = -1) | dept_id +  month_fe, 
                    data = df_sa, 
                    cluster = ~dept_id)
     
     sunab_ed_overall <- summary(sunab_ed, agg = "att")
# =============================================================================
# 2. Calloway and Sant'anna models
# =============================================================================

casa_read <- att_gt(data =df_cs,
                 yname = "pct_unplanned_reads",
                      tname = "month_index",
                      idname = "dept_id",
                      gname = "g",
                      cluster = "dept_id",
                      est_method = "dr",
                      allow_unbalanced = TRUE)

casa_group   <- aggte(casa_read, type = "group")
casa_overall <- aggte(casa_read, type = "simple")
casa_dynamic <- aggte(casa_read, type = "dynamic")

#------Event study plot 

    ggdid(casa_read)
    ggdid(casa_dynamic)
    
    casa_ed <- att_gt(data =df_cs,
                        yname = "pct_edvisits",
                        tname = "month_index",
                        idname = "dept_id",
                        gname = "g",
                        cluster = "dept_id",
                        est_method = "dr",
                        allow_unbalanced = TRUE)
    casa_ed_overall<- aggte(casa_ed, type = "simple")
    casa_ed_group   <- aggte(casa_ed, type = "group")
    
# =============================================================================
# 3. ITS model
# =============================================================================


    its_r <-  feols(pct_unplanned_reads ~ treated * time_pre + treated * post +
                treated * time_since | dept_id +  month_fe,
        data = df_its, 
        cluster = ~dept_id)
    
    its_ed <-  feols(pct_edvisits ~ treated * time_pre + treated * post +
                    treated * time_since | dept_id +  month_fe,
                  data = df_its, 
                  cluster = ~dept_id)
    

#------Checking how an event study plot looks with ITS model
      
        df_its <- df_its |>
          mutate(
            event_time = m24 - 13
          )
        
        
        es_mod <- feols(
          pct_unplanned_reads ~ 
            i(event_time, treated, ref = -1) |
            dept_id + month_fe,
          data = df_its,
          cluster = ~dept_id
        )
      
        png("R/group analysis/outputs/its_iplot.png", width = 800, height = 600)    
      its_iplot <-iplot(
          es_mod,
          ref.line = -1,
          xlab = "Months Relative to Intervention",
          ylab = "Effect on % Unplanned Readmissions",
          main = "Event Study: Multi-time period DiD"
        )
      dev.off()

#------Creating a better looking joint event study plot looks with Sunab and ITS model
      
      show<-ggiplot(list("ITS" = es_mod, "Sun & Abraham" = sunab),
                    main = "Effect on % Unplanned Readmissions: Split mutli-model",
                    xlab = "Months Relative to Intervention",
                    multi_style = "facet",
                    geom_style = "ribbon",
                    aggr_es = "none",
                    facet_args = list(
                      labeller = labeller(id = \(x) gsub(".*: ", "", x))) +
                      theme_minimal() +
                      theme(
                        text = element_text(family = "HersheySans"),
                        plot.title = element_text(hjust = 0.5),
                        legend.position = "none"
                      ))
      ggsave(filename = "R/group analysis/outputs/joint_iplot.png",
             plot = show)           
# =============================================================================
# Results Tables
# =============================================================================
  #Tables for staggered DiD
    modelsummary(
      list(sunab_overall, sunab_group, casa_overall, casa_group),
      conf_level = 0.90,
      coef_rename = c(
        "ATT" = "ATE",
        "ATT(Average)" = "ATE",
        "ATT(simple average)" = "ATE",
        "ATT(17)" = "cohort::17",
        "ATT(25)" = "cohort::25",
        "ATT(31)" = "cohort::31"
      ),
      statistic = "conf.int",
      gof_omit = 'DF|Deviance|R2|AIC|BIC',
      stars = c("*" = .10, "**" = .05, "***" = .01),
      fmt = 2
    )

      modelsummary(
        list(casa_ed_overall, sunab_ed_overall),
        conf_level = 0.90,
        # coef_rename = c(
        #   "ATT" = "ATE",
        #   "ATT(Average)" = "ATE",
        #   "ATT(simple average)" = "ATE",
        # ),
        statistic = "conf.int",
        gof_omit = 'DF|Deviance|R2|AIC|BIC',
        stars = c("*" = .10, "**" = .05, "***" = .01),
        fmt = 2
      )
      
      
  #Tables for ITS results
    modelsummary(
      list(its_r, its_ed),
      coef_map = c(
        "treated:post" = "ITT (level shift)",
        "treated:time_since" = "ITT (slope change)"
      ),
      conf_level = 0.90,
      statistic = "conf.int",
      gof_omit = 'DF|Deviance|R2|AIC|BIC',
      stars = c("*" = .10, "**" = .05, "***" = .01),
      fmt = 2
    )
# =============================================================================
# Forest Plots
# =============================================================================    
 ##Sunab estimates (overall and group)
    ci_r <- confint(sunab_overall, level = 0.95)
    sunab_overall_est <- data.frame(
      term      = names(coef(sunab_overall)),
      estimate  = as.numeric(coef(sunab_overall)),
      conf.low  = ci_r[, 1],
      conf.high = ci_r[, 2],
      model     = "B. Sun–Abraham estimator",
      outcome = "Readmissions"
    )
    
    ci_group_r <- confint(sunab_group, level = 0.95)
    sunab_group_est <- data.frame(
      term      = names(coef(sunab_group)),
      estimate  = as.numeric(coef(sunab_group)),
      conf.low  = ci_group_r[, 1],
      conf.high = ci_group_r[, 2],
      model     = "Sun–Abraham (Group)",
      outcome = "Readmissions"
    )
    
    ci_ed <- confint(sunab_ed_overall, level = 0.95)
    sunab_ed_overall_est <- data.frame(
      term      = names(coef(sunab_ed_overall)),
      estimate  = as.numeric(coef(sunab_ed_overall)),
      conf.low  = ci_ed[, 1],
      conf.high = ci_ed[, 2],
      model     = "B. Sun–Abraham estimator",
      outcome = "ED Visits"
    )
    
 ##ITS estimates 
    b  <- coef(its_r)
    se <- sqrt(diag(vcov(its_r)))
    z95 <- qnorm(0.975)
    
    its_est_r <- data.frame(
      term      = names(b),
      estimate  = as.numeric(b),
      conf.low  = as.numeric(b - z95 * se),
      conf.high = as.numeric(b + z95 * se),
      model     = "A. DiD with multiple time periods",
      outcome = "Readmissions"
    )
    
    b  <- coef(its_ed)
    se <- sqrt(diag(vcov(its_ed)))
    z95 <- qnorm(0.975)
    
  ##ITS estimates   
    its_est_ed <- data.frame(
      term      = names(b),
      estimate  = as.numeric(b),
      conf.low  = as.numeric(b - z95 * se),
      conf.high = as.numeric(b + z95 * se),
      model     = "A. DiD with multiple time periods",
      outcome = "ED Visits"
    )
 
  ##Casa estimates (overall and group)
    casa_overall_est <- data.frame(
      term      = paste0("ATT"),
      estimate  = casa_overall$overall.att,
      conf.low  = casa_overall$overall.att - z95 * casa_overall$overall.se,
      conf.high = casa_overall$overall.att + z95 * casa_overall$overall.se,
      model     = "C. Callaway-Sant'anna estimator",
      outcome = "Readmissions"
    )
    
    casa_group_est <- data.frame(
      term      = paste0("ATT(", casa_group$egt, ")"),
      estimate  = casa_group$att.egt,
      conf.low  = casa_group$att.egt - casa_group$crit.val * casa_group$se.egt,
      conf.high = casa_group$att.egt + casa_group$crit.val * casa_group$se.egt,
      model     = "D. Callaway-Sant'anna group estimator",
      outcome = "Readmissions"
    )
    
    casa_ed_overall_est <- data.frame(
      term      = paste0("ATT"),
      estimate  = casa_ed_overall$overall.att,
      conf.low  = casa_ed_overall$overall.att - z95 * casa_ed_overall$overall.se,
      conf.high = casa_ed_overall$overall.att + z95 * casa_ed_overall$overall.se,
      model     = "C. Callaway-Sant'anna estimator",
      outcome = "ED Visits"
    )
    
    casa_ed_group_est <- data.frame(
      term      = paste0("ATT(", casa_ed_group$egt, ")"),
      estimate  = casa_ed_group$att.egt,
      conf.low  = casa_ed_group$att.egt - casa_ed_group$crit.val * casa_ed_group$se.egt,
      conf.high = casa_ed_group$att.egt + casa_ed_group$crit.val * casa_ed_group$se.egt,
      model     = "D. Callaway-Sant'anna group estimator",
      outcome = "ED Visits"
    )
    
 #combine all estimates into one df
    est_all <- dplyr::bind_rows(
      its_est_r,
      its_est_ed,
      sunab_overall_est,
      sunab_ed_overall_est,
      casa_overall_est,
      casa_ed_overall_est,
      casa_group_est,
      casa_ed_group_est,
      .id = "source"
    ) 
    
    est_all <- est_all %>%
      mutate(
        model = factor(model, levels = c(
            "A. DiD with multiple time periods",
            "B. Sun–Abraham estimator",
            "C. Callaway-Sant'anna estimator",
            "D. Callaway-Sant'anna group estimator"
          )
        )
      ) %>%
      mutate(
        effect = case_when(
          term == "treated:post"                ~ "Immediate change",
          term == "treated:time_since"          ~ "Slope change",
          term == "ATT" ~ "Overall effect",
          term %in% c("cohort::17", "ATT(17)") ~ "Pilot",
          term %in% c("cohort::25", "ATT(25)") ~ "Main launch",
          term %in% c("cohort::31", "ATT(31)") ~ "Vascular",
          
          TRUE                                  ~ NA_character_
        )
      ) %>%
      filter(!is.na(effect)) %>%
      mutate(
        effect = factor(
          effect,
          levels = c("Vascular","Main launch","Pilot", "Overall effect", 
                     "Slope change" ,"Immediate change")))
    
    
# ##Forest plot of ITS - readmissions and ed visits only
#     
#         
#       forest_its<-ggplot(est_all %>% filter(model %in% c("ITS")),
#                      aes(x = estimate, y = effect, color = outcome)) +
#       geom_point(position = position_dodge(0.6), size = 2) +
#       geom_errorbarh(
#         aes(xmin = conf.low, xmax = conf.high),
#         position = position_dodge(0.6),
#         height = 0.2
#       ) +
#       geom_vline(xintercept = 0, linetype = "dashed") +
#       labs(
#         x = "% Outcome (95% CI)",
#         title = "Treatment effects from ITS models",
#         color ="Outcome",
#         y = NULL
#       ) +
#         theme_minimal(base_size = 16)
#     
#     print(forest_its)
#     ggsave("R/group analysis/outputs/forest_its.png", forest_its, width=8, height=4, dpi=300)
#     

    ##Forest plot of all effects together
    
    
    #est_all <- est_all %>% filter(source != 6)
    
    forest<-ggplot(est_all,
                   aes(x = estimate, y = effect, color = outcome, group = outcome)) +
      geom_point(position = position_dodge(0.6), size = 2) +
      geom_errorbarh(
        aes(xmin = conf.low, xmax = conf.high),
        position = position_dodge(0.6),
        height = 0.2
      ) +
      geom_vline(xintercept = 0, linetype = "dashed") +
      labs(
        x = "Effect size (pp with 95% CIs)",
        y = NULL,
       # title = "Treatment effects from staggered adoption estimators",
        #subtitle = "Percent unplanned readmissions",
        color ="Outcome"
      ) +
      theme_minimal(base_size = 14) +
      facet_wrap(~ model, ncol = 1, scales = "free_y") 
     # theme(legend.position = "none")
    
    print(forest)
    ggsave("R/group analysis/outputs/forest.png", forest, width=8, height=4, dpi=300)
    
    
# =============================================================================
#  Subgroup Analysis - ITS model
# =============================================================================
    
    #Model using surgical group aggregation
    its_ips_out_r <-  feols(pct_unplanned_reads ~ treated * time_pre +
                        treated * post + treated * time_since | out_id + month_fe,
                  data = df_its_ips_out, 
                  cluster = ~out_id)
    its_ips_out_ed <-  feols(pct_edvisits ~ treated * time_pre +
                            treated * post + treated * time_since | out_id + month_fe,
                          data = df_its_ips_out, 
                          cluster = ~out_id)
    
    #model using department aggregation
    its_ips_dept_r <-  feols(pct_unplanned_reads ~ treated * time_pre +
                       treated * post + treated * time_since | dept_id + month_fe,
                     data = df_its_ips_dept, 
                     cluster = ~dept_id)
    its_ips_dept_ed <-  feols(pct_edvisits ~ treated * time_pre +
                            treated * post + treated * time_since | dept_id + month_fe,
                          data = df_its_ips_dept, 
                          cluster = ~dept_id)
    
    #Tables for ITS results
    modelsummary(
      list(its_r, its_ips_dept_r, its_ips_out_r),
      coef_map = c(
        "treated:post" = "ITT (level shift)",
        "treated:time_since" = "ITT (slope change)"
      ),
      conf_level = 0.95,
      statistic = "conf.int",
      gof_omit = 'DF|Deviance|R2|AIC|BIC',
      stars = c("*" = .10, "**" = .05, "***" = .01),
      fmt = 2
    )
    modelsummary(
      list(its_ed, its_ips_dept_ed, its_ips_out_ed),
      coef_map = c(
        "treated:post" = "ITT (level shift)",
        "treated:time_since" = "ITT (slope change)"
      ),
      conf_level = 0.95,
      statistic = "conf.int",
      gof_omit = 'DF|Deviance|R2|AIC|BIC',
      stars = c("*" = .10, "**" = .05, "***" = .01),
      fmt = 2
    )


    # =============================================================================
    # Forest Plots for ITS subgroup analysis
    # =============================================================================     
    
    ##Main model - all patients -Unplanned readmissions and ED visits
    b  <- coef(its_r)
    se <- sqrt(diag(vcov(its_r)))
    z95 <- qnorm(0.975)
    
    its_est_r <- data.frame(
      term      = names(b),
      estimate  = as.numeric(b),
      conf.low  = as.numeric(b - z95 * se),
      conf.high = as.numeric(b + z95 * se),
      model     = "All patients",
      outcome = "Readmissions"
    )
    
    b  <- coef(its_ed)
    se <- sqrt(diag(vcov(its_ed)))
    z95 <- qnorm(0.975)
    
    its_est_ed <- data.frame(
      term      = names(b),
      estimate  = as.numeric(b),
      conf.low  = as.numeric(b - z95 * se),
      conf.high = as.numeric(b + z95 * se),
      model     = "All patients",
      outcome = "ED visits"
    )
    
    ##Surgical group  
    
    b  <- coef(its_ips_out_r)
    se <- sqrt(diag(vcov(its_ips_out_r)))
    z95 <- qnorm(0.975)
    
    its_est_out_r <- data.frame(
      term      = names(b),
      estimate  = as.numeric(b),
      conf.low  = as.numeric(b - z95 * se),
      conf.high = as.numeric(b + z95 * se),
      model     = "Inpatients (surgery group level)",
      outcome = "Readmissions"
    )   
    
    b  <- coef(its_ips_out_ed)
    se <- sqrt(diag(vcov(its_ips_out_ed)))
    z95 <- qnorm(0.975)
    
    its_est_out_ed <- data.frame(
      term      = names(b),
      estimate  = as.numeric(b),
      conf.low  = as.numeric(b - z95 * se),
      conf.high = as.numeric(b + z95 * se),
      model     = "Inpatients (surgery group level)",
      outcome = "ED visits"
    )   
    
    ##Department group  
    
    b  <- coef(its_ips_dept_r)
    se <- sqrt(diag(vcov(its_ips_dept_r)))
    z95 <- qnorm(0.975)
    
    its_est_dept_r <- data.frame(
      term      = names(b),
      estimate  = as.numeric(b),
      conf.low  = as.numeric(b - z95 * se),
      conf.high = as.numeric(b + z95 * se),
      model     = "Inpatients",
      outcome = "Readmissions"
    )    
    
    b  <- coef(its_ips_dept_ed)
    se <- sqrt(diag(vcov(its_ips_dept_ed)))
    z95 <- qnorm(0.975)
    
    its_est_dept_ed <- data.frame(
      term      = names(b),
      estimate  = as.numeric(b),
      conf.low  = as.numeric(b - z95 * se),
      conf.high = as.numeric(b + z95 * se),
      model     = "Inpatients",
      outcome = "ED visits"
    )    
    
    
    #combine all estimates into one df
    #est_all_its <- dplyr::bind_rows(its_est_r,its_est_dept_r, its_est_out_r,
                                    #its_est_ed,its_est_dept_ed, its_est_out_ed) 
    est_all_its <- dplyr::bind_rows(its_est_r,its_est_dept_r, 
                                    its_est_ed,its_est_dept_ed) 
    
    
    est_all_its <- est_all_its %>%
      mutate( model = factor(model, levels = c(
           "All patients",
           "Inpatients",
           "Inpatients (surgery group level)"  ))) %>%
      mutate(
        effect = case_when(
          term == "treated:post"                ~ "Immediate change",
          term == "treated:time_since"          ~ "Slope change",
          TRUE                                  ~ NA_character_
        )
      ) %>%
      filter(!is.na(effect)) %>%
      mutate(
        effect = factor(
          effect,
          levels = c("Slope change", "Immediate change" )))
    
    
    ##Forest plot of all effects together
    
    
    forest_subgroup <- ggplot(
      est_all_its,
      aes(
        x = estimate,
        y = effect,
        color = outcome,
        group = outcome
      )
    ) +
      geom_point(position = position_dodge(0.6), size = 2) +
      geom_errorbarh(
        aes(xmin = conf.low, xmax = conf.high),
        position = position_dodge(0.6),
        height = 0.2
      ) +
      facet_wrap(~ model, ncol = 1, scales = "free_y") +
      coord_cartesian(xlim = c(-10, 10)) +
      geom_vline(xintercept = 0, linetype = "dashed") +
      labs(
        x = "Effect size (pp with 95% CI)",
        y = NULL,
        #title = "Treatment effects for ITS models by sample",
        #subtitle = "Percent ED Visits and Readmissions",
        color = "Outcome"
      ) +
      theme_minimal(base_size = 12)
    
    print(forest_subgroup)
    
    ggsave("R/group analysis/outputs/forest_subgroup.png", forest_subgroup, width=8, height=4, dpi=300)
    

# =============================================================================
# Changing clustering at service level as a sensitivity analysis
# =============================================================================
    library(fixest)
    library(clubSandwich)
    
    #----Sun & Abraham
    sunab_service <- feols(
      pct_unplanned_reads ~ sunab(g, month_index, ref.p = -1) |
        dept_id + month_fe,
      data = df_sa)
    
    #Manual Wild bootstrap clustering
    B <- 999
    
    boot_coef <- replicate(B, {
      
      # cluster-level weights
      w <- sample(
        c(-1, 1),
        size = length(unique(df_sa$outreach_name)),
        replace = TRUE)
      
      # map weights to observations
      df_boot <- df_sa
      
      # create bootstrap outcome
      df_boot$y_boot <- fitted(sunab_service) +
        resid(sunab_service) * w[match(df_boot$outreach_name,
                             unique(df_boot$outreach_name))]
      
      boot_mod <- feols(
        y_boot ~ sunab(g, month_index, ref.p = -1) |
          dept_id + month_fe,
        data = df_boot
      )
      
      coef(boot_mod)
    })
    
    sunab_service_overall <- summary(sunab_service, agg = "att")   #----Sun and Abraham
  
    
    #----Callaway & Santanna
    casa_service <- att_gt(data =df_cs,
                        yname = "pct_unplanned_reads",
                        tname = "month_index",
                        idname = "dept_id",
                        gname = "g",
                        est_method = "dr",
                        allow_unbalanced = TRUE, 
                        bstrap = TRUE,
                        biters = 999,
                        cband = TRUE,
                        clustervars = "outreach_name")
    
    
    casa_group_service   <- aggte(casa_service, type = "group")
    casa_service_overall <- aggte(casa_service, type = "simple")

    
    
    #Tables for staggered DiD
    modelsummary(
      list(sunab_overall, sunab_service_overall, casa_overall, casa_service_overall, casa_group, casa_group_service),
      conf_level = 0.95,
      coef_rename = c(
        "ATT" = "ATE",
        "ATT(Average)" = "ATE",
        "ATT(simple average)" = "ATE",
        "ATT(17)" = "cohort::17",
        "ATT(25)" = "cohort::25",
        "ATT(31)" = "cohort::31"
      ),
      statistic = "conf.int",
      gof_omit = 'DF|Deviance|R2|AIC|BIC',
      stars = c("*" = .10, "**" = .05, "***" = .01),
      fmt = 2
    )
    
    ###############
    # Remove all unnecessary Data and Objects
    rm(list = setdiff(ls(), c("df_test", "df_analysis", "df_main", 
                              "df_matches", "df_pdc", "df_did",
                              "df_its", "df_sa", "df_cs",
                              "df_its_ips_out", "df_test_ips_out",
                              "df_its_ips_dept", "df_test_ips_dept")))
    
    
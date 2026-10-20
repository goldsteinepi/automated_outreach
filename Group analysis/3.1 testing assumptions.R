library(dplyr)
library(ggplot2)
library(scales)
library(did)
library(plm)
#install.packages("omnibus")
#library(omnibus)

# =============================================================================
### Departmental plots - pre-trends
# =============================================================================


# Aggregate to monthly department x treatment means (and SE for ribbons)
  dept_month <- df_test %>%
        filter(
          (outreach_name %in% c("Colorectal Surgery", "Urology Surgery") & month_index < 17) |  # pilot groups started earlier
            !(outreach_name %in% c("Colorectal Surgery", "Urology Surgery") 
            )) %>%
        filter(month_index<25) %>%
    group_by(outreach_name, treatment, month_index) %>%
    summarise(
      n         = n(),
      mean_pct  = mean(pct_unplanned_reads, na.rm = TRUE),
      sd_pct    = sd(pct_unplanned_reads, na.rm = TRUE),
      se_pct    = sd_pct / sqrt(n),
      .groups = "drop"
  )

plot_dept <- ggplot(dept_month, aes(x = month_index, y = mean_pct, color = treatment, group = treatment)) +
      geom_ribbon(aes(ymin = mean_pct - 1.96*se_pct,
                      ymax = mean_pct + 1.96*se_pct,
                      fill = treatment),
                  alpha = 0.15, color = NA, show.legend = FALSE) +
      geom_line(aes(color=treatment),linewidth = 1) +
      geom_point(size = 1.8) +
      scale_y_continuous(limits = c(0,50)) +
      scale_x_continuous(
        limits = c(0.5, 24.5),
        breaks = c(1, 24),
        labels = c("Jan '23",  "Dec '24"),
        expand = c(0, 0)
      ) + 
      scale_color_manual(
        values = c("Control" = "#1f77b4",  # blue
                   "Treatment" = "#d62728"), # red
        name = "Group"
      ) +
      labs(
        title = "Pre-trend Check (2024): Treatment vs Control by Department",
        subtitle = "Mean monthly unplanned readmission rate with 95% CI",
        y = "Unplanned Readmission Rate",
        color = "Group"
      ) +
      theme_minimal(base_size = 10) +
      theme(panel.grid = element_blank()) +
      facet_wrap(~ outreach_name)

ggsave("R/group analysis/outputs/bysurgery_pretrend.png", plot_dept, width=10, height=5, dpi=150)


############################################################################
#  Pre- trend plots for different samples
############################################################################

#------#Line plot unplanned readmissions by  treatment and control - Pre Trends
#A. Sample = staggered DiD sample 
#-----------------------------------------------------------------------------

df_line <- df_did %>%
  filter(year<2025)%>% 
  filter(
    !(month_index > 16 & outreach_name %in% c("Colorectal Surgery", "Urology Surgery"))) %>%
  group_by(treatment, month_index) %>%
  summarize(
    #pct_unplanned_reads = mean(unplanned_reads/n),
    totunplanned_reads = sum(unplanned_reads),
    totaln = sum(n),
    totedvisits = sum(edvisits),
    #avg_ureads = mean(unplanned_reads, na.rm = TRUE)
  ) %>%
  mutate(pct_unplanned_reads = (totunplanned_reads/totaln) *100, 
         pct_edvisits = (totedvisits/totaln)*100) %>%
  ungroup()


avg_lines <- df_line %>%
  group_by(treatment) %>%
  summarize(
    avg_pct_r = mean(pct_unplanned_reads, na.rm = TRUE),
    avg_ureads = mean(totunplanned_reads, na.rm = TRUE),
    avg_pct_ed = mean(pct_edvisits, na.rm = TRUE),
    avg_ed = mean(totedvisits, na.rm = TRUE)
  )

# -----Plotting -----
plot_pre_did_read<- ggplot(df_line, aes(x = month_index)) +
  geom_line(aes(y = pct_unplanned_reads, color = treatment), linewidth = 1.2) +
  geom_point(aes(y = pct_unplanned_reads, color = treatment), size = 2) +
  
  # Average % line
  geom_hline(data = avg_lines,
             aes(yintercept = avg_pct_r, color = treatment),
             linetype = "dashed", linewidth = 0.8) +
  
  scale_x_continuous(
    limits = c(0.5, 24.5),
    breaks = c(1, 6, 12, 18, 24),
    labels = c("Jan '23", "Jun '23", "Dec '23", "Jun '24", "Dec '24")
  ) +
  scale_y_continuous(
    name = "% Unplanned Readmissions",
    limits = c(0, 15)
  ) +
  scale_colour_manual(values = c("Treatment" = "steelblue",
                                 "Control"   = "coral")) +
  labs(
    title = "A. Unplanned Readmissions Trends",
    x = "Month",
    y = "Pct Unplanned Readmissions",
    color = "Group"
  ) +
  theme_minimal(base_size = 12) +
  theme(
    panel.grid.minor = element_blank(),
    legend.position = "bottom",
    axis.title.x = element_text(margin = margin(t = 10)),
    axis.title.y  = element_text(margin = margin(r = 10)),
  )
plot_pre_did_read
ggsave("R/group analysis/outputs/plot_pre_did_read.png", plot_pre_did_read, width=10, height=5, dpi=150)

plot_pre_did_ed<- ggplot(df_line, aes(x = month_index)) +
  geom_line(aes(y = pct_edvisits, color = treatment), linewidth = 1.2) +
  geom_point(aes(y = pct_edvisits, color = treatment), size = 2) +
  
  # Average % line
  geom_hline(data = avg_lines,
             aes(yintercept = avg_pct_ed, color = treatment),
             linetype = "dashed", linewidth = 0.8) +
  
  scale_x_continuous(
    limits = c(0.5, 24.5),
    breaks = c(1, 6, 12, 18, 24),
    labels = c("Jan '23", "Jun '23", "Dec '23", "Jun '24", "Dec '24")
  ) +
  scale_y_continuous(
    name = "% ED Visits",
    limits = c(0, 15)
  ) +
  scale_colour_manual(values = c("Treatment" = "steelblue",
                                 "Control"   = "coral")) +
labs(
  title = "B. ED Visits Trends",
  x = "Month",
  y = "Pct ED Visits",
  color = "Group"
) +
  theme_minimal(base_size = 12) +
  theme(
    panel.grid.minor = element_blank(),
    legend.position = "bottom",
    axis.title.x = element_text(margin = margin(t = 10)),
    axis.title.y  = element_text(margin = margin(r = 10)),
  )
plot_pre_did_ed
ggsave("R/group analysis/outputs/plot_pre_did_ed.png", plot_pre_did_ed, width=10, height=5, dpi=150)

#-----------------------------------------------------------------------------
#B. Sample = ITS 
#-----------------------------------------------------------------------------
df_line6 <- df_its %>% 
  filter(year==2024)%>% 
  filter(!outreach_name %in% c("Colorectal Surgery", "Urology Surgery",
                               "Vascular Surgery")) %>%  
  group_by(treatment, month_index) %>%
  summarize(
    #pct_unplanned_reads = mean(unplanned_reads/n),
    totunplanned_reads = sum(unplanned_reads),
    totaln = sum(n),
    totedvisits = sum(edvisits),
    #avg_ureads = mean(unplanned_reads, na.rm = TRUE)
  ) %>%
  mutate(pct_unplanned_reads = totunplanned_reads/totaln, 
         pct_edvisits = totedvisits/totaln) %>%
  ungroup()


avg_lines6 <- df_line6 %>%
  group_by(treatment) %>%
  summarize(
    avg_pct_r = mean(pct_unplanned_reads, na.rm = TRUE),
    avg_ureads = mean(totunplanned_reads, na.rm = TRUE),
    avg_pct_ed = mean(pct_edvisits, na.rm = TRUE),
    avg_ed = mean(totedvisits, na.rm = TRUE)
  )

# -----Plotting -----

#Readmissions
plot_pre_its_r<-ggplot(df_line6, aes(x = month_index)) +
  geom_line(aes(y = pct_unplanned_reads, color = treatment), linewidth = 1.2) +
  geom_point(aes(y = pct_unplanned_reads, color = treatment), size = 2) +
  
  # Average % line
  geom_hline(data = avg_lines6,
             aes(yintercept = avg_pct_r, color = treatment),
             linetype = "dashed", linewidth = 0.8) +
  scale_x_continuous(
    limits = c(12.5, 24.5),
    breaks = c(13, 18, 24),
    labels = c( "Jan '24", "Jun '24", "Dec '24")
  ) +
  scale_colour_manual(values = c("Treatment" = "steelblue",
                                 "Control"   = "coral")) +
  scale_y_continuous(
    name = "% Unplanned Readmissions",
    limits = c(0, 0.15)
  ) +
  labs(
    title = "Trend in % Unplanned Readmissions (Treatment vs Control) ",
    subtitle = "Excludes Pilot and Vascular Departments",
    x = "Month",
    y = "Pct Unplanned Readmissions",
    color = "Group"
  ) +
  theme_minimal(base_size = 16) +
  theme(
    panel.grid.minor = element_blank(),
    legend.position = "bottom"
  )

plot_pre_its_r
ggsave("R/group analysis/outputs/plot_pre_its_r.png", plot_pre_its_r, width=10, height=5, dpi=150)

#ED visits
plot_pre_its_ed<-ggplot(df_line6, aes(x = month_index)) +
  geom_line(aes(y = pct_edvisits, color = treatment), linewidth = 1.2) +
  geom_point(aes(y = pct_edvisits, color = treatment), size = 2) +
  
  # Average % line
  geom_hline(data = avg_lines6,
             aes(yintercept = avg_pct_ed, color = treatment),
             linetype = "dashed", linewidth = 0.8) +
  scale_x_continuous(
    limits = c(12.5, 24.5),
    breaks = c(13, 18, 24),
    labels = c( "Jan '24", "Jun '24", "Dec '24")
  ) +
  scale_y_continuous(
    name = "% ED Visits",
    limits = c(0, 0.15)
  ) +
  scale_colour_manual(values = c("Treatment" = "steelblue",
                                 "Control"   = "coral")) +
  labs(
    title = "Trend in % ED Visits (Treatment vs Control) ",
    subtitle = "Excludes Pilot and Vascular Departments",
    x = "Month",
    y = "Pct ED Visits",
    color = "Group"
  ) +
  theme_minimal(base_size = 16) +
  theme(
    panel.grid.minor = element_blank(),
    legend.position = "bottom"
  )

plot_pre_its_ed
ggsave("R/group analysis/outputs/plot_pre_its_ed.png", plot_pre_its_ed, width=10, height=5, dpi=150)

#-----------------------------------------------------------------------------     
#C. Sample = Staggered DiD sample -- individual cohorts
#-----------------------------------------------------------------------------

df_line7 <- df_did %>%
  filter(year<2025)%>% 
  filter(
    !(month_index > 16 & outreach_name %in% c("Colorectal Surgery", "Urology Surgery"))
  ) %>%
  group_by(g, month_index) %>%
  summarize(
    pct_unplanned_reads = mean(unplanned_reads/n),
    unplanned_reads = sum(unplanned_reads)
  ) %>%
  ungroup() %>%
  mutate(g_label = factor(
    g,
    levels = c(0, 17, 25, 31),
    labels = c(
      "Control",
      "Pilot",
      "Jan '25",
      "Vascular" )))


avg_lines7 <- df_line7 %>%
  group_by(g_label) %>%
  summarize(
    avg_pct = mean(pct_unplanned_reads, na.rm = TRUE),
    avg_ureads = mean(unplanned_reads, na.rm = TRUE)
  )

# -----Plotting -----
plot_pre_did_cohort <-ggplot(df_line7, aes(x = month_index)) +
  geom_line(aes(y = pct_unplanned_reads, color = g_label), linewidth = 1.2) +
  geom_point(aes(y = pct_unplanned_reads, color = g_label), size = 2) +
  
  # Average % line
  geom_hline(data = avg_lines7,
             aes(yintercept = avg_pct, color = g_label),
             linetype = "dashed", linewidth = 0.8) +
  
  scale_x_continuous(
    limits = c(0.5, 24.5),
    breaks = c(1, 6, 12, 18, 24),
    labels = c("Jan '23", "Jun '23", "Dec '23", "Jun '24", "Dec '24")
  ) +
  scale_y_continuous(
    name = "% Unplanned Readmissions",
    limits = c(0, 0.25)
  ) +
  labs(
    title = "Trend in % Unplanned Readmissions",
    x = "Month",
    y = "Pct Unplanned Readmissions",
    color = "Group"
  ) +
  theme_minimal(base_size = 11) +
  theme(
    panel.grid.minor = element_blank(),
    legend.position = "bottom"
  )
plot_pre_did_cohort
ggsave("R/group analysis/outputs/plot_pre_did_cohort.png", plot_pre_did_cohort, width=10, height=5, dpi=150)


#-----------------------------------------------------------------------------
#B. Sample = ITS sample & Outreach subgroup
#-----------------------------------------------------------------------------      
#df_test$high_outreach <- ifelse(df_test$pct_outreach >= 80, 1, 0)
avgoutreach <- df_test %>%  mutate(
  # First-treatment period (g): 0 = never treated, all other months are base off starting date of Jan '23
  g = case_when(
    treatment == "Control"                                         ~ 0L,
    outreach_name %in% c("Colorectal Surgery", "Urology Surgery") ~ 17L,
    outreach_name == "Vascular Surgery"                           ~ 31L,
    treatment == "Treatment"                                       ~ 25L,
    TRUE                                                           ~ NA_integer_
  ),
  month_index = (year - 2023) * 12 + month,  # new index variable for all months 2023-2026 we have
  # Relative event time (NA for never-treated)
  rel_time = ifelse(g > 0, month_index - g, NA_real_)) %>%
  group_by(department) %>%
  summarize(post_outreach = mean(pct_outreach[rel_time > 0], na.rm = TRUE)) %>%
  ungroup()



df_line2 <- merge(df_test, avgoutreach, by = "department") %>%
  group_by(department) %>%
  mutate(high_outreach = ifelse(is.na(post_outreach), 0,
                                ifelse(post_outreach > 80, 1, 0))) %>%  # binary variable to identify groups that received high outreach at any point
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
  ) %>%
  #write_xlsx(df_line2, "R/group analysis/df_line2.xlsx")
  
  filter(!outreach_name %in% c("Colorectal Surgery", "Urology Surgery",
                               "Vascular Surgery")) %>%  
  filter(year<2025)%>% 
  group_by(treatment, trt_outreach_group, m24) %>%
  summarize(
    pct_unplanned_reads = mean(unplanned_reads/n),
    unplanned_reads = sum(unplanned_reads)
  ) %>%
  ungroup()

table(df_line2$trt_outreach_group)

avg_lines2 <- df_line2 %>%
  group_by(treatment, trt_outreach_group) %>%
  summarize(
    avg_pct = mean(pct_unplanned_reads, na.rm = TRUE),
    avg_ureads = mean(unplanned_reads, na.rm = TRUE)
  )

# -----Plotting -----

ggplot(df_line2, aes(x = m24)) +
  geom_line(aes(y = pct_unplanned_reads, color = trt_outreach_group), linewidth = 1.2) +
  geom_point(aes(y = pct_unplanned_reads, color = trt_outreach_group), size = 2) +
  
  # Average % line
  geom_hline(data = avg_lines2,
             aes(yintercept = avg_pct, color = trt_outreach_group),
             linetype = "dashed", linewidth = 0.8) +
  scale_x_continuous(
    limits = c(0.5, 12.5),
    breaks = c(1, 6, 12, 18, 24),
    labels = c("Jan '24", "Jun '24", "Dec '24", "Jun '25", "Dec '25")
  ) +
  scale_y_continuous(
    name = "% Unplanned Readmissions",
    limits = c(0, 0.25)
  ) +
  labs(
    title = "Trend in % Unplanned Readmissions (Treatment & high oureach groups vs Control) ",
    subtitle = "Excludes Pilot and Vascular Departments",
    x = "Month",
    y = "Pct Unplanned Readmissions",
    color = "Group"
  ) +
  theme_minimal(base_size = 11) +
  theme(
    panel.grid.minor = element_blank(),
    legend.position = "bottom"
  )

#-----------------------------------------------------------------------------
#C. Sample = Staggered DiD & Outreach subgroup
#-----------------------------------------------------------------------------      
unique(df_sa$department[df_sa$trt_outreach_group == "Treatment, High Outreach"])
unique(df_sa$department[df_sa$trt_outreach_group == "Treatment, No High Outreach"])

df_line3 <- df_sa %>% #from script 3.2 preliminary analysis - Staggered DiD_S&A.R
  filter(year<2025)%>% 
  filter(
    !(month_index > 16 & outreach_name %in% c("Colorectal Surgery", "Urology Surgery"))
  ) %>%
  group_by(treatment, trt_outreach_group, month_index) %>%
  summarize(
    pct_unplanned_reads = mean(unplanned_reads/n),
    unplanned_reads = sum(unplanned_reads)
  ) %>%
  ungroup()
table(df_line3$trt_outreach_group)

avg_lines3 <- df_line3 %>%
  group_by(treatment, trt_outreach_group) %>%
  summarize(
    avg_pct = mean(pct_unplanned_reads, na.rm = TRUE),
    avg_ureads = mean(unplanned_reads, na.rm = TRUE)
  )

# -----Plotting -----

ggplot(df_line3, aes(x = month_index)) +
  geom_line(aes(y = pct_unplanned_reads, color = trt_outreach_group), linewidth = 1.2) +
  geom_point(aes(y = pct_unplanned_reads, color = trt_outreach_group), size = 2) +
  
  # Average % line
  geom_hline(data = avg_lines3,
             aes(yintercept = avg_pct, color = trt_outreach_group),
             linetype = "dashed", linewidth = 0.8) +
  
  scale_x_continuous(
    limits = c(0.5, 24.5),
    breaks = c(1, 6, 12, 18, 24),
    labels = c("Jan '23", "Jun '23", "Dec '23", "Jun '24", "Dec '24")
  ) +
  scale_y_continuous(
    name = "% Unplanned Readmissions",
    limits = c(0, 0.25)
  ) +
  labs(
    title = "Trend in % Unplanned Readmissions (Treatment & high oureach groups vs Control) ",
    x = "Month",
    y = "Pct Unplanned Readmissions",
    color = "Group"
  ) +
  theme_minimal(base_size = 11) +
  theme(
    panel.grid.minor = element_blank(),
    legend.position = "bottom"
  )

#-----------------------------------------------------------------------------
#C2. Sample = Staggered DiD & CCI subgroup
#-----------------------------------------------------------------------------      
unique(df_sa$department[df_sa$scenario_label == "Scenario 1: Treatment, CCI avg above median"])
unique(df_sa$department[df_sa$scenario_label == "Scenario 2: Treatment, CCI avg below median"])

df_line4 <- df_sa %>% #from script 3.2 preliminary analysis - Staggered DiD_S&A.R
  filter(year<2025)%>% 
  filter(
    !(month_index > 16 & outreach_name %in% c("Colorectal Surgery", "Urology Surgery"))
  ) %>%
  group_by(treatment, scenario_label, month_index) %>%
  summarize(
    pct_unplanned_reads = mean(unplanned_reads/n),
    unplanned_reads = sum(unplanned_reads)
  ) %>%
  ungroup()
table(df_line4$scenario_label)

avg_lines4 <- df_line4 %>%
  group_by(treatment, scenario_label) %>%
  summarize(
    avg_pct = mean(pct_unplanned_reads, na.rm = TRUE),
    avg_ureads = mean(unplanned_reads, na.rm = TRUE)
  )

# -----Plotting -----

ggplot(df_line4, aes(x = month_index)) +
  geom_line(aes(y = pct_unplanned_reads, color = scenario_label), linewidth = 1.2) +
  geom_point(aes(y = pct_unplanned_reads, color = scenario_label), size = 2) +
  
  # Average % line
  geom_hline(data = avg_lines4,
             aes(yintercept = avg_pct, color = scenario_label),
             linetype = "dashed", linewidth = 0.8) +
  scale_x_continuous(
    limits = c(0.5, 24.5),
    breaks = c(1, 6, 12, 18, 24),
    labels = c("Jan '23", "Jun '23", "Dec '23", "Jun '24", "Dec '24")
  ) +
  scale_y_continuous(
    name = "% Unplanned Readmissions",
    limits = c(0, 0.25)
  ) +
  labs(
    title = "Trend in % Unplanned Readmissions (Treatment & CCI groups vs Control) ",
    x = "Month",
    y = "Pct Unplanned Readmissions",
    color = "Group"
  ) +
  theme_minimal(base_size = 11) +
  theme(
    panel.grid.minor = element_blank(),
    legend.position = "bottom"
  )

#-----------------------------------------------------------------------------
#D. Sample = ITS & CCI subgroup
#-----------------------------------------------------------------------------
df_line5 <- df_sa %>% #from script 3.2 preliminary analysis - Staggered DiD_S&A.R
  filter(year<2025)%>% 
  filter(!outreach_name %in% c("Colorectal Surgery", "Urology Surgery",
                               "Vascular Surgery")) %>%  
  group_by(treatment, scenario_label, month_index) %>%
  summarize(
    pct_unplanned_reads = mean(unplanned_reads/n),
    unplanned_reads = sum(unplanned_reads)
  ) %>%
  ungroup()
table(df_line5$scenario_label)


avg_lines5 <- df_line5 %>%
  group_by(treatment, scenario_label) %>%
  summarize(
    avg_pct = mean(pct_unplanned_reads, na.rm = TRUE),
    avg_ureads = mean(unplanned_reads, na.rm = TRUE)
  )

# -----Plotting -----

ggplot(df_line5, aes(x = month_index)) +
  geom_line(aes(y = pct_unplanned_reads, color = scenario_label), linewidth = 1.2) +
  geom_point(aes(y = pct_unplanned_reads, color = scenario_label), size = 2) +
  
  # Average % line
  geom_hline(data = avg_lines5,
             aes(yintercept = avg_pct, color = scenario_label),
             linetype = "dashed", linewidth = 0.8) +
  scale_x_continuous(
    limits = c(12.5, 24.5),
    breaks = c(13, 18, 24),
    labels = c( "Jan '24", "Jun '24", "Dec '24")
  ) +
  scale_y_continuous(
    name = "% Unplanned Readmissions",
    limits = c(0, 0.25)
  ) +
  labs(
    title = "Trend in % Unplanned Readmissions (Treatment & CCI groups vs Control) ",
    x = "Month",
    y = "Pct Unplanned Readmissions",
    color = "Group"
  ) +
  theme_minimal(base_size = 11) +
  theme(
    panel.grid.minor = element_blank(),
    legend.position = "bottom"
  )





# =============================================================================
## Visualizing covariate trends by treatment and control 
# =============================================================================
      #grouping by outreach name
      covar_trends_name <- df_did %>%
        group_by(month_index, outreach_name) %>%
        summarise(
          mean_pct_ureads = mean(pct_unplanned_reads, na.rm = TRUE),
          mean_cci    = mean(avg_cci,          na.rm = TRUE),
          mean_pct_or = mean(pct_or_dash_case, na.rm = TRUE),
          mean_age    = mean(avg_age,          na.rm = TRUE),
          mean_female = mean(pct_female,       na.rm = TRUE),
          mean_pct_outreach = mean(pct_outreach),
          mean_los = mean(los, na.rm = TRUE),
          mean_disrrs = mean(dis_rrs, na.rm = TRUE),
          .groups = "drop"
        ) 

    #grouping by staggered did cohorts

      covar_trends_did <- df_did %>%
        group_by(month_index, g) %>%
        summarise(
          mean_pct_ureads = mean(pct_unplanned_reads, na.rm = TRUE),
          mean_cci    = mean(avg_cci,          na.rm = TRUE),
          mean_pct_or = mean(pct_or_dash_case, na.rm = TRUE),
          mean_age    = mean(avg_age,          na.rm = TRUE),
          mean_female = mean(pct_female,       na.rm = TRUE),
          mean_pct_outreach = mean(pct_outreach,   na.rm = TRUE),
          mean_los = mean(los, na.rm = TRUE),
          mean_disrrs = mean(dis_rrs, na.rm = TRUE),
          .groups = "drop"
        ) %>%
        mutate(cohort = factor(g,
                               levels = c(0, 17, 25, 31),
                               labels = c("Control (g=0)", "Pilot (g=17)",
                                          "Jan-25 (g=25)", "Vascular (g=31)")))
      # By ITS groupings
      covar_trends_its <- df_its %>%
        group_by(m24, treated) %>%
        summarise(
          mean_pct_ureads = mean(pct_unplanned_reads, na.rm = TRUE),
          mean_cci    = mean(avg_cci,          na.rm = TRUE),
          mean_pct_or = mean(pct_or_dash_case, na.rm = TRUE),
          mean_age    = mean(avg_age,          na.rm = TRUE),
          mean_female = mean(pct_female,       na.rm = TRUE),
          mean_pct_outreach = mean(pct_outreach,   na.rm = TRUE),
          mean_los = mean(los, na.rm = TRUE),
          mean_disrrs = mean(dis_rrs, na.rm = TRUE),
          .groups = "drop"
        )  %>%
        mutate(treatment = factor(treated,
                                  levels = c(0,1),
                                  labels = c("Control", "Treatment")))
      
      
      make_cov_plot <- function(data, month, yvar, trtgroup, ytitle) {
        ggplot(data, aes_string(x = month, y = yvar,
                                colour = trtgroup, group = trtgroup)) +
          
          geom_line() + geom_point(size = 1.5) +
          labs(title = paste("Covariate Trend:", ytitle),
               x = month, y = ytitle, 
               colour = trtgroup ) +
          theme_minimal(base_size = 11)
      }
      
      p_cov_outreach <- make_cov_plot(covar_trends_name, "month_index", "mean_pct_outreach", "outreach_name", "% Automated Outreach" )    

      
      
      p_cov_did <- make_cov_plot(covar_trends_did, "month_index",  "mean_cci",  "cohort",  "Avg CCI")    /
        make_cov_plot(covar_trends_did, "month_index", "mean_age",  "cohort",  "Avg Age")    /
        make_cov_plot(covar_trends_did, "month_index", "mean_female", "cohort", "% Female")
      
      p_cov2_did <- make_cov_plot(covar_trends_did, "month_index","mean_pct_or", "cohort",   "% OR cases")    /
        make_cov_plot(covar_trends_did,"month_index", "mean_pct_outreach",  "cohort",  "% Outreach cases")    
      
      p_cov3_did <- make_cov_plot(covar_trends_did, "month_index", "mean_pct_ureads", "cohort","% Unplanned Readmissions") +
        scale_y_continuous(lim=c(0,45)) 
      
      p_cov4_did <- make_cov_plot(covar_trends_did, "month_index", "mean_los", "cohort","Avg length of Stay") /
        make_cov_plot(covar_trends_did, "month_index", "mean_disrrs", "cohort","Avg discharge readmission risk score") 
      #+    scale_y_continuous(lim=c(7,17)) 
      
      
      p_cov_its <- make_cov_plot(covar_trends_its, "m24",  "mean_cci",  "treatment",  "Avg CCI")    /
        make_cov_plot(covar_trends_its, "m24", "mean_age",  "treatment",  "Avg Age")    /
        make_cov_plot(covar_trends_its, "m24", "mean_female", "treatment", "% Female")
      
      p_cov2_its <- make_cov_plot(covar_trends_its, "m24","mean_pct_or", "treatment",   "% OR cases")    /
        make_cov_plot(covar_trends_its,"m24", "mean_pct_outreach",  "treatment",  "% Outreach cases")    
      
      p_cov3_its <- make_cov_plot(covar_trends_its, "m24", "mean_pct_ureads", "treatment","% Unplanned Readmissions") +
        scale_y_continuous(lim=c(0,45))  
      
      p_cov4_its <- make_cov_plot(covar_trends_its, "m24", "mean_los", "treatment","Avg length of Stay") /
        make_cov_plot(covar_trends_its, "m24", "mean_disrrs", "treatment","Avg discharge readmission risk score") 
      
      ggsave("R/group analysis/outputs/outreachcovariate_trends.png", p_cov_outreach, width = 9, height = 11, dpi = 300)
      
      ggsave("R/group analysis/outputs/didcovariate_trends.png", p_cov_did, width = 9, height = 11, dpi = 300)
      ggsave("R/group analysis/outputs/didcovariate_trends2.png", p_cov2_did, width = 9, height = 11, dpi = 300)
      ggsave("R/group analysis/outputs/didcovariate_trends4.png", p_cov4_did, width = 9, height = 11, dpi = 300)
      
      ggsave("R/group analysis/outputs/itscovariate_trends.png", p_cov_its, width = 9, height = 11, dpi = 300)
      ggsave("R/group analysis/outputs/itscovariate_trends2.png", p_cov2_its, width = 9, height = 11, dpi = 300)
      ggsave("R/group analysis/outputs/itscovariate_trends4.png", p_cov4_its, width = 9, height = 11, dpi = 300)
      
    #########################################
    
      rm(list = setdiff(ls(), c("df_test", "df_analysis", "df_main", "df_matches", "df_pdc", "df_did","df_its", "df_sa")))
      
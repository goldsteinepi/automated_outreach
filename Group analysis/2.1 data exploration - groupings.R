#Set-up analysis file 

library(tidyverse) 
library(janitor)
library(readxl)
#install.packages("tidytext") 
library(tidytext)

setwd("C:/Users/kkilburn/OneDrive - Geisinger System Services/GPE and BIT Shared Directory-GPE shared data (PHI) - Documents/Automated Discharge/Data/")


#Look at numbers of patients in each surgery grouping


df_counts <- df_analysis %>%
  count(year, outreach_name, treatment) %>%
  group_by(year, outreach_name) %>%
  mutate(total_n = sum(n)) %>%
  ungroup() %>%
  arrange(desc(total_n))

ggplot(df_counts, aes(x = fct_reorder(outreach_name, n, .desc = TRUE), y = n,
                      fill = factor(treatment, levels = c(0,1), labels = c("Control", "Treatment")))) +
  geom_col() +
  labs(
    title = "Number of discharges by surgery group (2024-2025)",
    x = "Surgery group",
    y = "Discharges",
    fill = "Treament"
  ) +
  theme_minimal(base_size = 14) +
  theme(
    strip.background = element_rect(fill = "lightgray", color = NA)
  ) +
  coord_flip() +
  facet_wrap(~ year, ncol = 1, scales = "free_y")


### Graphs of discharges by department

df_counts2 <- df_analysis %>%
  count(department, treatment, year) %>%
  mutate(
    treat_lbl = factor(treatment, levels = c(0, 1),
                       labels = c("Control", "Treatment")),
    year = factor(year, levels = c(2024, 2025), labels = c("2024", "2025")),  # 2025 at the bottom of stack
  )

# Build an x-axis order driven by Control totals (then Treatment)
order_tbl <- df_counts2 %>%
  group_by(department, treat_lbl) %>%
  summarise(n_sum = sum(n), .groups = "drop") %>%
  tidyr::pivot_wider(names_from = treat_lbl, values_from = n_sum, values_fill = 0) %>%
  mutate(
    ctrl = `Control`,
    trt  = `Treatment`
  ) %>%
  arrange(desc(ctrl), desc(trt)) %>%
  pull(department)

df_counts <- df_counts %>%
  mutate(outreach_name = factor(outreach_name, levels = unique(order_tbl)))

pos <- position_dodge2(width = 0.8, preserve = "single")

df_counts2$department <- stringr::str_wrap(df_counts2$department, width = 25)

ggplot(df_counts2,
       aes(x = reorder_within(department, n, treat_lbl), y = n, fill = year)) +
  scale_x_reordered()

ggplot(
  df_counts2,
  aes(
    x = department,
    y = n,
    fill = year,
    #color = treat_lbl,
    group = treat_lbl
  )
) +
  geom_col(position = pos, width = 0.7, linewidth = 0.8) +
  scale_fill_manual(values = c( "2024" = "#E69F00", "2025" = "#56B4E9"), drop = FALSE) +
  labs(
    title = "Discharges by department",
    x = "Surgery department",
    y = "Discharges",
    fill = "Year",
    color = "Treatment"
  ) +
  theme_minimal(base_size = 10) +
  theme(panel.grid.major.y = element_blank()) +
  coord_flip() + facet_wrap(~ treat_lbl, ncol = 2, scales = "free_y")
``
### Heat map of discharges by department and ratios by month
df_heat<-df_group %>%
  mutate(m24 = ifelse(year==2025, month+12, month))

ggplot(df_heat, aes(x = m24, y = department, fill = n)) +
  geom_tile(color = "white") +
  scale_x_continuous(
    limits = c(0.5, 24.5),
    breaks = c(1, 6, 12, 18, 24),
    labels = c("Jan '24", "June '24", "Dec '24", "June '25", "Dec '25"),
    expand = c(0, 0)
  ) +
  scale_fill_viridis_c(option = "C") +
  labs(
    title = "Monthly Discharges by Department (2024-2025)",
    x = "Month",
    y = "Department",
    fill = "Discharges"
  ) +
  theme_minimal(base_size = 10) +
  theme(panel.grid = element_blank())  +
  geom_vline(xintercept = 12.5, linetype = "dashed", color = "red", linewidth = 2)

ggplot(df_heat, 
       aes(x = m24, y = department, fill = ratio_dis)) +
  geom_tile(color = "white") +
  scale_fill_viridis_c(option = "C") + 
  scale_x_continuous(
    limits = c(0.5, 24.5),
    breaks = c(1, 6, 12, 18, 24),
    labels = c("Jan '24", "June '24", "Dec '24", "June '25", "Dec '25"),
    expand = c(0, 0)
  ) +
  labs( title = "Monthly Discharge Ratios by Department", 
        x = "Month", y = "Department", fill = "Ratios" ) +
  theme_minimal(base_size = 10) +
  geom_vline(xintercept = 12.5, linetype = "dashed", color = "red", linewidth = 2)


## Looking at the percent of individual discharges that are also found on the OR dashboard


ggplot(
  df_group,
  aes(
    x = fct_reorder(department, pct_or_dash_case, .desc = TRUE),
    y = pct_or_dash_case,
    #color = factor(year, levels = c(2024, 2025), labels = c("2024", "2025"))
    color = factor(treatment, levels = c(0,1), labels = c("Control", "Treatment"))
  )
) +
  geom_point(position = position_dodge(width = 0.6), size = 2) +
  labs(
    title = "Percent of patients overlap in OR dashboard (2024-2025)",
    x = "Surgery group",
    y = "Percent OR dash cases",
    color = "Treatment"
  ) +
  theme_minimal(base_size = 10) +
  theme(strip.background = element_rect(fill = "lightgray", color = NA)) +
  coord_flip() + 
  facet_wrap(~ treatment, ncol = 1, scales = "free_y"
             ,                 # each panel gets its own category set
             space  = "free_y"
  )


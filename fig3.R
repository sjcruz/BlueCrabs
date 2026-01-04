library(dplyr)
library(ggplot2)
library(ggpubr)
source("functions/paper_graph_labels.R")

data = read.csv(paste0(getwd(), "/results/res3.csv"))
#999 codes no recovery 

data = data %>%
  mutate( recovery_delta_months = recovery_time_months - recovery_time_3months_baseline) %>%
  mutate(
    months_saved_bin = case_when(
      recovery_time_months != 999 & recovery_time_3months_baseline == 999 ~ "Recovers with policy\nbut not under status quo",
      recovery_time_months == 999 & recovery_time_3months_baseline != 999 ~ "Recovers with policy\nbut not under status quo",
      recovery_time_months == 999 & recovery_time_3months_baseline == 999 ~ "No recovery with policy\nor status quo",
      abs(recovery_delta_months) < 3 | is.na(recovery_delta_months) ~ "No difference",
      recovery_delta_months <= -36 ~ ">36 months faster",
      recovery_delta_months <= -12 ~ "≤36 months faster",   # (i.e., 12–36 months faster)
      recovery_delta_months <= -3  ~ "≤12 months faster",    # (i.e., 3–12 months faster)
      # Slower recovery (positive deltas)
      recovery_delta_months <= 12  ~ "≤12 months slower",
      recovery_delta_months <= 36  ~ "≤36 months slower",
      recovery_delta_months > 36   ~ ">36 months slower"),
    
    # Legend / plotting order
    months_saved_bin = factor(months_saved_bin,levels = c(
      "No recovery with policy\nor status quo", "Recovers with policy\nbut not under status quo",
      ">36 months faster", "≤36 months faster", "≤12 months faster",
      "No difference",
      "≤12 months slower", "≤36 months slower", ">36 months slower")),
    policy = factor( policy,
                     levels = c("Increase male size limit", "No fishing females", "No fishing males", "Full closure")),
    climate_shock_labels = factor(climate_shock_labels, levels = rev(climate_order)),
    summer = ifelse(climate_shock_labels%in%summer, "Summer", "Other"), 
    summer = factor(summer, levels = c("Summer", "Other")))



ggarrange(data %>%
            subset(type == "stock") %>%
            mutate(months_saved_bin) %>%
            ggplot(aes(x = factor(scenario), y = climate_shock_labels, fill = months_saved_bin)) +
            geom_tile() +
            labs(x = x_temp_lab, y = "Timing of warming", subtitle = "A") +
            theme_classic() +
            facet_grid(summer~policy, switch = "y", scales ="free_y", space="free_y")+
            base_theme+ 
            geom_vline(xintercept = 7.5)+
            theme(plot.title = element_text(face = "bold", size=10), panel.spacing.y = unit(0, "pt")) +
            custom_fill_months_saved(),
          data %>%
            subset(type == "NPV") %>%
            ggplot(aes(x = factor(scenario), y = climate_shock_labels, fill = months_saved_bin)) +
            geom_tile() +
            labs(x = x_temp_lab, y = "Timing of warming", subtitle = "B") +
            theme_classic() +
            base_theme+ 
            facet_grid(summer~policy, switch = "y", scales ="free_y", space="free_y")+
            base_theme+ 
            geom_vline(xintercept = 7.5)+
            theme(plot.title = element_text(face = "bold", size=10),  panel.spacing.y = unit(0, "pt")) +
            custom_fill_months_saved(),
          
          nrow = 2, common.legend=T, legend = "right")

#ggsave(filename = "fig3.jpg", width = 7.25, units="in", dpi=1200)
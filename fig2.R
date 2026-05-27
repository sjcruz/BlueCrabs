
source("functions/paper_graph_labels.R")


difference = fread(paste0(getwd(), "/results/res2_a.csv"))
difference$climate_shock_labels = factor(difference$climate_shock_labels, levels =rev(climate_order))
difference$optimized = factor(difference$optimized, levels = c("No fishing", "Optimized fishing, status quo policy"))
difference$summer = factor(difference$summer, levels = c("Summer", "Other"))
difference$percent_decline_cat = factor(difference$percent_decline_cat, levels = percent_change_levels)


recovery = fread(paste0(getwd(), "/results/res2_b.csv"))
recovery$climate_shock_labels = factor(recovery$climate_shock_labels, levels =rev(climate_order))
recovery$optimized = factor(recovery$optimized, levels = c("No fishing", "Optimized fishing, status quo policy"))
recovery$summer = factor(recovery$summer, levels = c("Summer", "Other"))
recovery$recovery_time_cat = factor(recovery$recovery_time_cat, levels = recovery_levels)

ggarrange(
  difference %>%
    ggplot(aes(x = factor(scenario), y = climate_shock_labels, fill = percent_decline_cat)) +
    geom_tile() +
    labs(x = x_temp_lab, y = "Timing of warming", title="a") +
    facet_grid(summer~optimized, switch = "y", scales ="free_y", space="free_y") +
    theme_classic() +
    geom_vline(xintercept = 7.5)+
    base_theme +
    theme(plot.title = element_text(face = "bold"),  panel.spacing.y = unit(0, "pt")) +
    custom_fill_difference (name=paste0("% Change stock\nrelative to baseline")), 
  recovery %>%
    ggplot(aes(x = factor(scenario), y = climate_shock_labels, fill = recovery_time_cat)) +
    geom_tile() +
    geom_vline(xintercept = 7.5)+
    labs(x = x_temp_lab, y = "Timing of warming", title="b") +
    facet_grid(summer~optimized, switch = "y", scales ="free_y", space="free_y") +
    theme_classic() +
    base_theme +
    theme(plot.title = element_text(face = "bold"),  panel.spacing.y = unit(0, "pt")) +
    custom_fill_recovery (name=paste("Time to stock\nrecovery (months)")),
  nrow=2)


#ggsave(filename = "fig2.pdf", width = 180, units="mm", dpi=1200)
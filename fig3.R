
source("functions/paper_graph_labels.R")

data = fread(paste0(getwd(), "/results/res3.csv"))

data = data %>%
    mutate (months_saved_bin = factor(months_saved_bin,levels = time_to_recovery_levels),
    policy = factor( policy, levels =policy_labels),
    climate_shock_labels = factor(climate_shock_labels, levels = rev(climate_order)),
    summer = factor(summer, levels = c("Summer", "Other")))

ggarrange(data %>%
            subset(type == "stock") %>%
            mutate(months_saved_bin) %>%
            ggplot(aes(x = factor(scenario), y = climate_shock_labels, fill = months_saved_bin)) +
            geom_tile() +
            labs(x = x_temp_lab, y = "Timing of warming", subtitle = "a") +
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
            labs(x = x_temp_lab, y = "Timing of warming", subtitle = "b") +
            theme_classic() +
            base_theme+ 
            facet_grid(summer~policy, switch = "y", scales ="free_y", space="free_y")+
            base_theme+ 
            geom_vline(xintercept = 7.5)+
            theme(plot.title = element_text(face = "bold", size=10),  panel.spacing.y = unit(0, "pt")) +
            custom_fill_months_saved(),
          
          nrow = 2, common.legend=T, legend = "right")

#ggsave(filename = "fig3.pdf", width = 180, units="mm", dpi=1200)
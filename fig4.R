library(dplyr)
library(ggplot2)
library(ggpubr)
library(patchwork)
library(data.table)
source("functions/paper_graph_labels.R")
size =5
dark2_colors = RColorBrewer::brewer.pal(4, "Dark2")
policy_colors = c("Status quo" = "black", "Full closure" = dark2_colors[4], "No fishing males" = dark2_colors[3],
                   "No fishing females" = dark2_colors[2], "Increase male size limit" = dark2_colors[1])

summary = fread(paste0(getwd(), "/results/res4_a_c.csv")) %>%
  mutate(policy = factor(policy, levels = rev(c("Status quo","Increase male size limit", "No fishing females" , "No fishing males", "Full closure"))))


p1 = summary %>%
  ggplot(aes(y = policy, x = abs(mean_cost),  fill=policy)) +
  geom_col(position = "dodge") +
  geom_errorbarh(aes(xmin = abs(mean_cost) - (1.96*se_cost), xmax = abs(mean_cost) + (1.96*se_cost)), height = 0.2) +
  scale_x_continuous(labels = scales::label_percent(accuracy = 1), expand=c(0,0)) +
  labs(y = "Policy type", x = "Net cost of policy") +
  scale_fill_manual(values = policy_colors, name = "Policy") +
  theme_classic(base_size = size) +
  theme( legend.position = "none",
         strip.text.y.right = element_text(angle = 90)) +
  facet_wrap(~temp_band, labeller = label_parsed, nrow = 2, strip.position = "left")

p2 = summary %>%
  mutate(se_recovery_time=ifelse(policy == "Status quo", 0, se_recovery_time))%>%
  ggplot(aes(y = policy, x = mean_recovery_time, fill=policy) )+
  geom_col(position = "dodge") +
  geom_errorbarh(aes(xmin = mean_recovery_time - (1.96*se_recovery_time),
                     xmax = mean_recovery_time + (1.96*se_recovery_time)),
                 height = 0.2) +
  labs(x = "Change in time to\nrecovery (months)") +
  theme_classic(base_size = size) +
  scale_fill_manual(values = policy_colors, name = "Policy") +
  scale_x_continuous(expand=c(0,0)) +
  geom_vline(xintercept = 0)+
  theme( legend.position = "none",
         axis.title.y = element_blank(),
         axis.text.y = element_blank(),
         strip.text = element_blank(),
         strip.background = element_blank()) +
  facet_wrap(~temp_band, labeller = label_parsed, nrow = 2)

p3 = summary %>%
  ggplot(aes(y = policy, x = recovery_proportion, fill=policy)) +
  geom_col(position = "dodge") +
  geom_errorbarh(aes(xmin = recovery_proportion - (1.96*se_recovery_proportion),
                     xmax = recovery_proportion + (1.96*se_recovery_proportion)),
                 height = 0.2) +
  labs(x = "Proporton of scenarios\nthat recover") +
  theme_classic(base_size = size) +
  scale_x_continuous(expand=c(0,0), breaks=c(0,0.25, 0.5, 0.75,1), labels=c(0,0.25, 0.5, 0.75,1)) +
  scale_fill_manual(values = policy_colors, name = "Policy") +
  theme( legend.position = "none",
         axis.title.y = element_blank(),
         axis.text.y = element_blank(),
         strip.text = element_blank(),
         strip.background = element_blank()) +
  facet_wrap(~temp_band, labeller = label_parsed, nrow = 2)


# Combine plots with relative widths

summary_plot = (p1+labs(subtitle = "A") | p2+labs(subtitle = "B")|p3+labs(subtitle = "C")) + plot_layout(widths = c(1, 1, 1))


###########################################################
# Survival curve plot 
data = fread(paste0(getwd(), "/results/res4_d.csv")) %>%
  mutate(policy = factor(policy, levels = rev(c("Status quo","Increase male size limit", "No fishing females" , "No fishing males", "Full closure"))))

recovery_thresholds = seq(6, 12*3, by = 6)
threshold_colors = colorRampPalette(c("black", "lightblue"))(length(recovery_thresholds))

recovery_summary =lapply(recovery_thresholds, function(thresh) {
  data %>%
    filter(recovery_time_months != 999 & policy != "Increase male size limit") %>%
    group_by(policy_length, temp_band) %>%
    summarise(
      threshold = thresh,
      pr_recovery = mean(recovery_time_months < thresh),
      se = sqrt(pr_recovery * (1 - pr_recovery) / n()),
      .groups = "drop"
    )
}) %>%
  bind_rows() %>%
  mutate(threshold = factor(threshold,
                            levels = recovery_thresholds,
                            labels = paste0(recovery_thresholds)))

plot_d = ggplot(recovery_summary,
                aes(x = policy_length, y = pr_recovery, color = threshold)) +
  geom_line(aes(group = threshold)) +
  geom_point(size=0.5) +
  geom_ribbon(aes(ymin = pmax(0, pr_recovery -  se),
                  ymax = pmin(1, pr_recovery +  se),
                  fill = threshold,
                  group = threshold),
              alpha = 0.2, color = NA) +
  scale_color_manual(values = threshold_colors) +
  scale_fill_manual(values = threshold_colors) +
  scale_y_continuous(labels = scales::label_percent(accuracy = 1), limits = c(0, 1),expand = expansion(mult = c(0, 0))) +
  scale_x_continuous(breaks = seq(12, 84, 12), minor_breaks = NULL,
                     expand = expansion(mult = c(0, 0))) +
  labs(x = "Length of policy (months)",
       y = "Proportion of scenarios recovered within threshold",
       color = "Recovery target\n(months)",
       fill  = "Recovery target\n(months)",
       subtitle = "D") +
  theme_classic(base_size =size) +
  theme(legend.position = "bottom", base_size =size) +
  facet_wrap(~ temp_band, labeller = label_parsed)



##################################################

e = read.csv(paste0(getwd(), "/results/res4_e.csv"))%>%
  mutate(policy = factor(policy, levels = (c("Increase male size limit", "No fishing females" , "No fishing males", "Full closure"))))


plot_e = e %>%
  ggplot() +
  geom_smooth(aes(policy_length,
                  cost_effective,
                  color = policy, fill=policy),
              se = TRUE) +
  labs(x = "Policy length (months)", y = "Cost of one month faster recovery (000 2020 USD)",
       color = "Policy", linetype = "", fill="Policy") +
  facet_wrap(~ temp_band, scales = "free",labeller = label_parsed) +
  scale_y_continuous(expand=c(0,0))+
  scale_x_continuous(limits=c(6,84), expand=c(0,0), breaks = seq(12, 84, 12))+
  scale_color_manual(values=policy_colors)+
  scale_fill_manual(values=policy_colors)+
  theme_classic(base_size = size) +
  theme(legend.position = "bottom")+
  geom_hline(yintercept =0)


##################################################
odds = read.csv(paste0(getwd(), "/results/res4_f.csv"))%>%
  mutate(policy = factor(policy, levels = (c("Increase male size limit", "No fishing females" , "No fishing males", "Full closure"))))

plot_f = odds %>%
  ggplot(aes(x = policy_length, y = (abs(percent_decline)/(pr_recovery))/1000, color = policy, fill=policy)) +
  geom_smooth(method = NULL, se = T)+
  scale_colour_brewer(palette = "Dark2", name = "Policy") +
  scale_fill_brewer(palette = "Dark2", name = "Policy") +
  labs(y="Cost of increasing odds of\nrecovery (000 2020 USD)", x= "Policy length (months)") +
  scale_y_continuous(expand=c(0,0))+
  scale_x_continuous(limits=c(6,84), expand=c(0,0), breaks = seq(12, 84, 12))+
  theme_classic(base_size = size)+
  geom_hline(yintercept = 0)+
  theme(legend.position = "bottom")+
  facet_wrap(~temp_band, labeller = label_parsed, nrow = 2)


ggarrange(summary_plot,
          ggarrange(plot_d, ggarrange(plot_e+labs(subtitle = "E")+ guides(color = guide_legend(nrow = 2, byrow = TRUE)),
          plot_f+labs(subtitle = "F")+ guides(color = guide_legend(nrow = 2, byrow = TRUE)), nrow = 1, common.legend = T,  legend = "bottom"), 
          nrow=1) , nrow=2)

#ggsave(filename = "fig4.jpg", width = 7.25, units="in", dpi=1200)

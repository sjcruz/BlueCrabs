
source("functions/paper_graph_labels.R")

size = fread(paste0(getwd(), "/results/res1_d.csv"))

size %>%
  mutate( month = factor(month), scenario = factor(scenario, levels = c(0,5), labels = c(0,5))) %>%
  ggplot(aes(x = CW,y = month,fill = scenario,color = scenario,group = interaction(month, scenario))) +
  geom_density_ridges(alpha = 0.2, scale = 1.5, rel_min_height = 0.01, linewidth = 0.4, bandwidth =2.5) +
  scale_fill_manual(values = c('0' = "blue", '5' = "orange")) +
  scale_color_manual(values = c('0' = "blue", '5' = "orange")) +
  labs( x = "Carapace width (mm)", y = "Month", fill = expression(Delta*' T ('*degree*C*')'), color = expression(Delta*' T ('*degree*C*')')) +
  scale_y_discrete(expand = c(0, 0)) +
  scale_x_continuous(expand = c(0, 0) ) +
  theme_bw() +
  theme(legend.position = "top") +base_theme

#ggsave(filename = "Fig1_d.png", height = 100, width = 70, units="mm")

data = fread(paste0(getwd(), "/results/res1_c.csv"))

data %>%
  mutate(market_cat= factor(market_cat, levels = c("Large males (>5.5in)", "Small males (5-5.5in)", "Females (dredge)",
                                                   "Females (pots)", "Peelers"))) %>%
  ggplot()+
  geom_line(aes(x=month, y=mean/1000, color=scenario), size=1)+
  geom_point(aes(x=month, y=mean/1000, color=scenario), size=1)+
  geom_ribbon(aes(ymin=min/1000, ymax=max/1000, x=month, fill=scenario), alpha=0.3)+
  scale_color_manual(values = c('actual' = "black", 'model' = "orange"), labels = c("actual" = "Actual", "model" = "Model"))+
  scale_fill_manual(values = c('actual' = "black", 'model' = "orange"), labels = c("actual" = "Actual", "model" = "Model"))+
  
  theme_classic() + year+
  labs(y="Harvest (000 pounds)",x="Month", color= " ", fill=" ")+
  facet_wrap(~market_cat, scales="free", nrow=2)+
  theme(panel.grid = element_blank(),
        panel.border = element_blank(), legend.position = c(0.9, 0.2), legend.justification = c(1,0)) +
  base_theme



#ggsave(filename = "Fig1_c.png", height = 60, width = 110, units="mm")

## Loads packages
library(dplyr)
library(ggplot2)
library(ggpubr)
library(patchwork)
library(data.table)
library(ggridges)

mon_let = c('J','F','M','A','M','J','J','A','S','O','N','D')
year = scale_x_continuous(breaks=seq(1,12,1), labels=c(rep(mon_let,1)), expand=(c(0,0)))

### Labels 
x_temp_lab = expression('Warming intensity ('*Delta*' T '*degree*C*')')
policy_labels = c("Increase male size limit", "No fishing females", "No fishing males", "Full closure")
dark2_colors = RColorBrewer::brewer.pal(4, "Dark2")
policy_colors = c("Status quo" = "black", "Full closure" = dark2_colors[4], "No fishing males" = dark2_colors[3],
                  "No fishing females" = dark2_colors[2], "Increase male size limit" = dark2_colors[1])

warm_labels = data.frame(climate_shock = c("all","year5", "year","jan_aug","feb_aug",  "jun_dec", "mar_aug", "apr_aug", "may_aug", "jun_sep","jun_oct","jun_nov", 
                                           "winter", "spring", "summer", "fall", "june", "july", "august", "jun_jul", "jul_aug"),
                         climate_shock_labels = c("persistent","Jan-Dec", "Jan-Dec", "Jan-Aug", "Feb-Aug", "Jun-Dec", "Mar-Aug", "Apr-Aug", "May-Aug", "Jun-Sep", "Jun-Oct", "Jun-Nov",
                                                  "Dec-Feb", "Mar-May", "Jun-Aug", "Sep-Nov",  "June", "July", "August", "Jun-Jul", "Jul-Aug"))


climate_order = c("persistent","Jan-Dec", "Jan-Aug","Feb-Aug", "Jun-Dec", "Mar-Aug","Jun-Nov", "Apr-Aug","May-Aug","Jun-Sep","Jun-Oct",
                  "Jun-Aug", "Jun-Jul","Jul-Aug","June","July","August", "Mar-May","Sep-Nov","Dec-Feb")

summer = c("persistent","Jan-Dec", "Jan-Aug","Feb-Aug", "Jun-Dec", "Mar-Aug","Jun-Nov", "Apr-Aug","May-Aug","Jun-Sep","Jun-Oct",
           "Jun-Aug", "Jun-Jul", "June","July","August")

percent_change_levels = c("< -90%", "-90% to -70%", "-70% to -50%", "-50% to -30%", "-30% to -10%", "-10% to 0%", "0% to 10%", "> 10%")
categorize_percent_decline <- function(x) {
  cut( x,
    breaks = c(-100, -0.9, -0.7, -0.5, -0.3, -0.1, 0, 0.1, 1),
    labels = percent_change_levels,
    include.lowest = TRUE)}

custom_fill_difference <- function(name = "% change\nrelative to baseline") {
  scale_fill_manual(
    values = c(
      "< -90%"      = "#7F0000",  # dark burgundy
      "-90% to -70%" = "#B30000",
      "-70% to -50%" = "#D73027",
      "-50% to -30%" = "#F46D43",
      "-30% to -10%" = "#FDAE61",
      "-10% to 0%"   = "#FFFFBF", # pale yellow
      "0% to 10%"    = "#66C2A5", # teal
      "> 10%"        = "#3288BD"  # blue-teal
    ),
    name = name
  )
}

recovery_levels = c("< 3", "3–12", "12–24", "24–36", "36–60", "60–100", "No recovery")
categorize_recovery_time <- function(x) {
  cut(
    x,
    breaks = c(0, 3, 12, 24, 36, 60, 100, 999),
    labels = recovery_levels,
    levels = recovery_levels,
    include.lowest = TRUE,
    right = FALSE
  )
}

custom_fill_recovery <- function(name = "Recovery time (months)") {
  scale_fill_manual(
    values = c("< 3" = "#3288BD",
               "3–12" = "#66C2A5",
               "12–24" = "#FFFFBF",
               "24–36" = "#FDAE61",
               "36–60" = "#F46D43",
               "60–100" = "#D73027",
               "No recovery" = "#7F0000"),
    name = name
  )
}

text_size = 7
base_theme <- theme(
  text = element_text(size = text_size),
  axis.text = element_text(size = text_size),
  axis.title = element_text(size = text_size),
  legend.text = element_text(size = text_size),
  legend.title = element_text(size = text_size),
  legend.key.size = unit(0.5, "cm"),  # Adjust as needed
  strip.text = element_text(size = text_size) # for facet labels
)

time_to_recovery_levels = c("No recovery with policy\nor status quo", "Recovers with policy\nbut not under status quo",
                            ">36 months faster", "≤36 months faster", "≤12 months faster", "No difference",
                            "≤12 months slower", "≤36 months slower", ">36 months slower")

custom_fill_months_saved <- function(name = "Change in time to recovery") {
  scale_fill_manual(
    values = c(
      "Recovers with policy\nbut not under status quo" = "#542788", # blue
      "No recovery with policy\nor status quo"         = "#000000", # black
      ">36 months faster"  = "#2166AC", # dark blue
      "≤36 months faster"  = "#67A9CF", # medium blue
      "≤12 months faster"  = "#D1E5F0", # light blue
      "No difference"      = "#BDBDBD", # gray
      "≤12 months slower"  = "#FDDBC7", # light orange
      "≤36 months slower"  = "#EF8A62", # orange-red
      ">36 months slower"  = "#B2182B"  # dark red
    ),
    name = name,
    drop = FALSE
  )
}


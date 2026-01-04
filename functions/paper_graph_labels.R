### Labels 
x_temp_lab = expression('Warming intensity ('*Delta*' T '*degree*C*')')
model_cats = c("all",  "males", "females", "large_males", "small_males", "females_pot", "peeler", "females_dredge")
optimize_cats = c('large males'  , 'pot females' , 'dredge females')
categories = c("large males","small males", "pot females", "peelers", "dredge females")
model_cat_labels = c("all", "males", "females", categories)
policy = c("status quo", "males_size", "females_closed", "males_closed", "all_closed")
policy_labels = c("Status quo","Increase male size limit", "No fishing females", "No fishing males", "Full closure")

warm_labels = data.frame(climate_shock = c("all","year5", "year","jan_aug","feb_aug",  "jun_dec", "mar_aug", "apr_aug", "may_aug", "jun_sep","jun_oct","jun_nov", 
                                           "winter", "spring", "summer", "fall", "june", "july", "august", "jun_jul", "jul_aug"),
                         climate_shock_labels = c("persistent","Jan-Dec", "Jan-Dec", "Jan-Aug", "Feb-Aug", "Jun-Dec", "Mar-Aug", "Apr-Aug", "May-Aug", "Jun-Sep", "Jun-Oct", "Jun-Nov",
                                                  "Dec-Feb", "Mar-May", "Jun-Aug", "Sep-Nov",  "June", "July", "August", "Jun-Jul", "Jul-Aug"))


climate_order = c("persistent","Jan-Dec", "Jan-Aug","Feb-Aug", "Jun-Dec", "Mar-Aug","Jun-Nov", "Apr-Aug","May-Aug","Jun-Sep","Jun-Oct",
                  "Jun-Aug", "Jun-Jul","Jul-Aug","June","July","August", "Mar-May","Sep-Nov","Dec-Feb")

summer = c("persistent","Jan-Dec", "Jan-Aug","Feb-Aug", "Jun-Dec", "Mar-Aug","Jun-Nov", "Apr-Aug","May-Aug","Jun-Sep","Jun-Oct",
           "Jun-Aug", "Jun-Jul", "June","July","August")

#climate_order <- c("persistent","Jan-Dec","Jan-Aug", "Feb-Aug", "Mar-May", "Mar-Aug", "Apr-Aug", "May-Aug", 
#                  "Jun-Aug", "Jun-Sep", "June", "July", "August", "Jun-Jul", "Jul-Aug" , "Jun-Oct", "Jun-Nov", "Jun-Dec",
#                 "Sep-Nov", "Dec-Feb")

categorize_percent_decline <- function(x) {
  cut(
    x,
    breaks = c(-100, -0.9, -0.7, -0.5, -0.3, -0.1, 0, 0.1, 1),
    labels = c("< -90%", "-90% to -70%", "-70% to -50%", "-50% to -30%", "-30% to -10%", "-10% to 0%", 
               "0% to 10%", "> 10%"),
    include.lowest = TRUE
  )
}

custom_fill_difference <- function(name = "% change\nrelative to baseline") {
  scale_fill_manual(
    values = c("< -90%" = "darkred",
               "-90% to -70%" = "red3",
               "-70% to -50%" = "red",
               "-50% to -30%" = "chocolate3",
               "-30% to -10%" = "orange3",
               "-10% to 0%"  = "yellow",
               "0% to 10%"   = "lightgreen",
               "> 10%"       = "springgreen3"),
    name = name
  )
}

categorize_recovery_time <- function(x) {
  cut(
    x,
    breaks = c(0, 3, 12, 24, 36, 60, 100, 999),
    labels = c("< 3", "3–12", "12–24", 
               "24–36", "36–60", "60–100", "No recovery"),
    include.lowest = TRUE,
    right = FALSE
  )
}

custom_fill_recovery <- function(name = "Recovery time (months)") {
  scale_fill_manual(
    values = c("< 3" = "darkgreen",
               "3–12" = "lightgreen",
               "12–24" = "yellow",
               "24–36" = "orange",
               "36–60" = "chocolate",
               "60–100" = "red",
               "No recovery" = "darkred"),
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

custom_fill_months_saved <- function(name = "Change in time to recovery") {
  scale_fill_manual(
    values = c(
      "Recovers with policy\nbut not under status quo"  = "dodgerblue3",
      "No recovery with policy\nor status quo" = "black",
      ">36 months faster"       = "darkgreen",
      "≤36 months faster"       = "green3",
      "≤12 months faster"       = "lightgreen",
      "No difference"           = "gray70",
      "≤12 months slower"       = "gold",
      "≤36 months slower"       = "orange",
      ">36 months slower"       = "darkred"
    ),
    name = name,
    drop = FALSE
  )
}


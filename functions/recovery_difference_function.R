

recovery_difference_function = function(data) {
  groups = c('market_cat','type', "scenario", "optimized", "policy", "climate_shock", "discount_rate", "start_policy_month", "policy_length")
 
  data= data.table(data)
  data = data[month > 24]  # removed first two years when population stabillizes 

  
  # Aggregate all data categories before calculating the difference
  all_data = data[, .(mean = sum(mean), base_mean = sum(base_mean)), by = setdiff(groups, "market_cat")][, market_cat := "all"]
  market_cat = data[, .(mean = sum(mean), base_mean = sum(base_mean)), by = groups]

  # Modify market categories for sex-based aggregation
  # Aggregate without market_cat first
  sex_data = data[, sex := fifelse(market_cat %in% c("small_males", "large_males"), "males",
                               fifelse(market_cat == "peeler", "peeler", "females"))]
  
  sex_data =  data[, .(mean = sum(mean), base_mean = sum(base_mean)), by = c("sex", setdiff(groups, "market_cat"))][, market_cat := sex][, sex := NULL][] 

  # Merge the three datasets
  difference = rbindlist(list(all_data, market_cat, sex_data), use.names = TRUE)
  
  # Compute group-level difference relative to baseline
  difference[, `:=` (difference = mean - base_mean)]
  difference[, `:=` (percent_decline = (mean - base_mean)/base_mean)]
  

  ####################### Compute recovery times
  # Merge climate shock last month data
  ### 1. Aggregate Data by `month` and `groups`
  all_data = data[, .(mean = sum(mean), base_mean = sum(base_mean)), by = c("month", setdiff(groups, "market_cat"))][, market_cat := "all"]  # Exclude market_cat from grouping
  
  market_cat_data = data[, .(mean = sum(mean), base_mean = sum(base_mean)), by = c("month", groups)]
  
  sex_data = data[, sex := fifelse(market_cat %in% c("small_males", "large_males"), "males",
                                   fifelse(market_cat == "peeler", "peeler", "females"))]
  
  sex_data =  data[, .(mean = sum(mean), base_mean = sum(base_mean)), by = c("sex", "month", setdiff(groups, "market_cat"))][, market_cat := sex][, sex := NULL][] 
  
  
  ### 2. Compute Recovery Time for Each Category
  compute_recovery_time = function(dt) {
    # 1. Subset to post-climate shock period
    dt = dt[month >= start_policy_month]
    
    # 2. Flag months where scenario mean exceeds baseline
    dt[, above_base := mean >= base_mean]
    
    dt_recovery = unique(dt[, .SD[1], by = groups])
    
    # now doing the same for 3,6,12 months 
    for (w in c(1, 3,6,12)){
      dt[, rolling_sum := rollapply(above_base, width = w, FUN = sum, align = "left", fill = NA), by = groups]
      
      # Get the w-th month in the window (i.e., last month of recovery window) ie recovery happens onthe 3rd motnh in a 3 month widnow 
      dt[, recovery_month := shift(month, type = "lead", n = w - 1), by = groups]
      
      # Extract earliest recovery month where condition is met
      recov = dt[rolling_sum == w, .(recovery = min(recovery_month, na.rm = TRUE)), by = groups]
      
      # Calculate recovery time relative to shock
      recov[, paste0("recovery_time_months_window", w) := recovery - start_policy_month]
      recov[, recovery := NULL]  # Drop intermediate
      
      # Merge into dt_recovery
      dt_recovery = merge(dt_recovery, recov, by = groups, all.x = TRUE)
    }
   
    return(dt_recovery)
  }
  
  all_recovery = compute_recovery_time(all_data)
  market_cat_recovery = compute_recovery_time(market_cat_data)
  sex_recovery = compute_recovery_time(sex_data)
  
  ### 3. Ensure All Groups are Included
  # Extract unique group combinations for each recovery dataset
  all_groups_recovery = unique(all_data[, ..groups])
  market_groups_recovery = unique(market_cat_data[, ..groups])
  sex_groups_recovery = unique(sex_data[, ..groups])
  
  # Merge with the appropriate unique groups and fill missing values with 999
  all_recovery = merge(all_groups_recovery, all_recovery, by = groups, all.x = TRUE)
  all_recovery[, recovery_time_months_window3 := fifelse(is.na(recovery_time_months_window3), 999, recovery_time_months_window3)]
  
  market_cat_recovery = merge(market_groups_recovery, market_cat_recovery, by = groups, all.x = TRUE)
  market_cat_recovery[, recovery_time_months_window3 := fifelse(is.na(recovery_time_months_window3), 999, recovery_time_months_window3)]
  
  sex_recovery = merge(sex_groups_recovery, sex_recovery, by = groups, all.x = TRUE)
  sex_recovery[, recovery_time_months_window3 := fifelse(is.na(recovery_time_months_window3), 999, recovery_time_months_window3)]
  
  
  ### 4. Merge All Recovery Data
  recovery = rbindlist(list(all_recovery, market_cat_recovery, sex_recovery), use.names = TRUE, fill = TRUE)
  
  return(list(difference, recovery))
}
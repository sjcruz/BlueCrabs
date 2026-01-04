################################################################################
# discount function 
################################################################################
discount = function (r){
  delta = (1+r)^(1/12)-1  #discount rate for each time step
  
  R = rep(0,months)
  
  for (tt in 1:months){       #each time period discount factor
    R[tt] = (1/(1+delta))^(tt-1) # 1/(1+delta)^t
  }
  
  return(R)
}

################################################################################
# Price function 
# Inputs 
# 1. a = choke price 
# 2. b = slope of demand 
# 3. Har = harvest 
# 4. n = # of fishers 
# Output is a matrix timeXi of prices
################################################################################

p = function (a, b, har, proportion_DE_CB, endo){
  if(endo ==1){
    row_sum = rowSums(har)
    
    LM = (((row_sum/proportion_DE_CB)+row_sum)/a)^(b)
    LM[LM<0] = 0
    LM[LM>choke[1]] = choke[1]
    LM = cbind(LM, LM, LM, LM, LM)
    pr = price_ratios * LM
    
  } else {
    pr=prices
  }
  return(pr)
  
}

################################################################################
# policy_variables
# updates a policy variables based on scenario being simulated including 
# closed seasons, size limtis, and catchability coefficients 
################################################################################


policy_variables = function(policy, start_policy_month, end_policy_month, policy_length){
  #############################  Status quo ##################################
  closed_season = function(closed_months, years){
    if(years == 1){
      return(closed_months)
    } else {
      out = closed_months
      for(i in 2:years-1){
        out = c(out, closed_months+(12*i))
      }
      return(out)
    }}
  
  closed_pot = closed_season(closed_months = c(1, 2, 3, 4,11,12), years=years)
  closed_pot_f = closed_season(closed_months = c(1, 2, 3,12), years=years)
  closed_dredge = closed_season(closed_months = c(4:11), years=years)
  closed_seasons = c(closed_pot, closed_pot_f+months, closed_dredge+months+months)
  closed_seasons_no_policy = closed_seasons
  
  if(policy_length>0){
    policy_closure_males = union(closed_pot, start_policy_month:end_policy_month)
    policy_closure_females = union(closed_pot_f, start_policy_month:end_policy_month)
    policy_closure_dredge = union(closed_dredge, start_policy_month:end_policy_month)
  }
  
  hardshell_size = rep(127, months) 
  
  #catchability coefficients including q tidas for each market category
  q_dredge =  as.matrix(fread(paste0(inputs_dir, "dredge_catchability.csv"), drop="V1"))
  q_dredge[,5] = 180
  q_dredge = do.call(rbind, replicate(years, q_dredge, simplify=FALSE))
  
  q_pot_males = c(0.45,0.65,0.5,0.1,0.0)
  q_pot_males = do.call(rbind, replicate(months, q_pot_males, simplify=FALSE))
  
  q_pot_females = c(0.05,0.2,1.3,0.1,0)
  q_pot_females = do.call(rbind, replicate(months, q_pot_females, simplify=FALSE))
  
  if (policy == "males_closed" & policy_length>0){
    
    closed_seasons = c(policy_closure_males, closed_pot_f+months, closed_dredge+months+months)
    q_pot_males[start_policy_month:end_policy_month,1:2] = 0 
    q_pot_females[start_policy_month:end_policy_month, 1:2] = 0
    
  } else if (policy == "females_closed" & policy_length>0){
    
    closed_seasons = c(closed_pot, policy_closure_females+months, policy_closure_dredge+months+months) 
    q_pot_males[start_policy_month:end_policy_month,3] = 0 
    q_pot_females[start_policy_month:end_policy_month,] = 0
    q_dredge[start_policy_month:end_policy_month,] = 0
    
  } else if (policy == "all_closed" & policy_length>0){
    # ensures no effort during females_closed - we assume perfect compliance 
    closed_seasons = c(policy_closure_males, policy_closure_females+months, policy_closure_dredge+months+months) 
    q_pot_males[start_policy_month:end_policy_month,] = 0 
    q_pot_females[start_policy_month:end_policy_month,] = 0
    q_dredge[start_policy_month:end_policy_month,] = 0
    
  } else if (policy == "males_size" & policy_length>0){
    # ensures no effort during females_closed - we assume perfect compliance 
    hardshell_size[start_policy_month:end_policy_month] = largemale_cutoff
    
  }
  
  return(list(
    closed_seasons = closed_seasons,
    q_pot_males = q_pot_males,
    q_pot_females = q_pot_females,
    q_dredge = q_dredge,
    hardshell_size = hardshell_size,
    closed_seasons_no_policy=closed_seasons_no_policy
  ))
  
}


################################################################################
# Objective function to be optimized 
# Inputs 
# 1. E = effort
# Cost 
# Output is sum discounted NPV across all time 


objective_function = function(E){
  #source ("Variable_setup.R")
  #E =  scan(file = "results/res_actual.txt", what = numeric(), sep = " ") 
  #E =  scan(file = "results/res_1_exo_scenario0_04-Jul-2024 02.10.txt", what = numeric(), sep = " ")
  effort = effort_matrix(E=E)
  shared_env$mod = biological_model(N_0 = N_0, effort = effort)  # Store result in the environment
  mod = shared_env$mod
  mod$recruits
  har = mod[[5]]* lbs_conversion *model_ratio
  price = p(a=a, b=b, har=har, proportion_DE_CB=proportion_DE_CB, endo=endo)

  effort = effort *model_ratio
  pot_males = mod[[6]] * lbs_conversion*model_ratio
  pot_females = mod[[7]]  * lbs_conversion*model_ratio
  dredge = mod[[8]] * lbs_conversion*model_ratio

  profit_males = apply(price*pot_males,1,sum) - (effort[,1]*cost[,1])
  profit_females = apply(price*pot_females,1,sum) - (effort[,2]*cost[,2])
  profit_dredge = apply(price*dredge,1,sum) - (effort[,3]*cost[,3])
  
  profit_males[is.nan(profit_males)] <- 0
  profit_females[is.nan(profit_females)] <- 0
  profit_dredge[is.nan(profit_dredge)] <- 0

  profit = apply(cbind(profit_males, profit_females, profit_dredge),1,sum)  # sum by month (row=1)
  NPV = R*profit
  #print(sum(NPV))
  return(-sum(NPV)) #negative to maximize 
}


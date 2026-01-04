################################################################################
# Constraint function
# From nloptr documentation: 
# The constrain function: to evaluate (non-)linear inequality constraints that should hold in the solution. 
# constraints must be formulated as g(x)<= 0

inequality_constraints  = function (E){
 
  effort = effort_matrix(E=E)
  mod = shared_env$mod  # Access stored model output
  har = mod[[5]]* lbs_conversion *model_ratio
  price = p(a=a, b=b, har=har, proportion_DE_CB=proportion_DE_CB, endo=endo)
  effort = effort *model_ratio
  pot_males = mod[[6]] * lbs_conversion*model_ratio
  pot_females = mod[[7]]  * lbs_conversion*model_ratio
  dredge = mod[[8]] * lbs_conversion*model_ratio
  
  profit_males = apply(price*pot_males,1,sum) - (effort[,1]*cost[,1])
  profit_females = apply(price*pot_females,1,sum) - (effort[,2]*cost[,2])
  profit_dredge = apply(price*dredge,1,sum) - (effort[,3]*cost[,3])
  
  # Handle NaN values
  profit_males[is.nan(profit_males)] = 0
  profit_females[is.nan(profit_females)] = 0
  profit_dredge[is.nan(profit_dredge)] = 0
  
  ef_control = apply(effort[,1:2],1,sum)-Emax_pot
  return(c(-profit_males, -profit_females, -profit_males, ef_control))
}

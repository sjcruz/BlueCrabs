################################################################################
# Seleni Cruz 
# Jan, 2026
# Set up to run in MPI parallel processing in HPC environment 
# Dimensions: temp_change (warming intensity), policy, climate shock (timing of wamring) and policy length
# Output is a text file of E_star for each simulation ran 
################################################################################

require(dplyr)
require(tidyr)
require(nloptr)
require(data.table)

data_location =  paste0(getwd(), "/outputs/")
home = paste0(getwd(), "/")
inputs_dir = paste0(home, "inputs/")
source ("Variable_setup.R")

# Parse command-line arguments
args = commandArgs(trailingOnly = TRUE)
slurm_task_id = as.numeric(args[1])  # The SLURM task ID for batch processing

discountrate = 0.03 #also ran for 0 and 0.01

# Define temperature change and policy scenarios
temp_changes = 0:10  # Vector of possible temp_change values (11 values)
policies = c("males_size", "females_closed", "males_closed", "all_closed", "none")  # "none" is status quo

target_year= 3  #warming starts in year 3, allows cohorts and population to stabilize 
start_day = (target_year - 1) * 365 + 1
climate_shock_scen <- list(
  fall = start_day + 244:334,
  summer = start_day + 172:243,
  winter = 731:789,
  spring = start_day + 59:171,
  may_aug = start_day + 120:243,
  apr_aug = start_day + 90:243,
  mar_aug = start_day + 59:243,
  feb_aug = start_day + 31:243,
  jan_aug = start_day + 0:243,
  jun_sep = start_day + 151:273,
  jun_oct = start_day + 151:304,
  jun_nov = start_day + 151:334,
  jun_dec = start_day + 151:364,
  june = start_day + 151:180,
  july = start_day + 181:211,
  august = start_day + 212:242,
  jun_jul = start_day + 151:211,
  year = start_day:(start_day + 364)
  #all = start_day:(365 * 10)
)



num_temp = length(temp_changes)  # 11
num_policy = length(policies)  # 5
num_climate_shock = length(climate_shock_scen)  # 19

# Compute total number of jobs
total_jobs = num_temp * num_policy * num_climate_shock 

# Number of simulations per SLURM job
batch_size = 16  # We run 16 scenarios per SLURM job, can adjust depending on workers and availability 

# Compute which 16 jobs this SLURM task should handle
start_index = slurm_task_id * batch_size
end_index = min(start_index + batch_size - 1, total_jobs - 1)

# Compute job-specific indices
for (job_index in start_index:end_index) {
  # Convert linear index to (temp_change, policy, climate_shock) combination
  temp_index = ((job_index - 1) %% num_temp) + 1  # 1-based indexing
  policy_index = (((job_index - 1) %/% num_temp) %% num_policy) + 1
  climate_index = (((job_index - 1) %/% (num_temp * num_policy)) %% num_climate_shock) + 1
  
  # Extract values
  temp_change = temp_changes[temp_index]
  policy = policies[policy_index]
  climate_shock = climate_shock_scen[[climate_index]]  # Use double brackets
  climate_shock_name = names(climate_shock_scen)[climate_index]
  
  cat(paste0("SLURM Assigned Job ID: ", slurm_task_id, 
             " -> Running scenario: temp_change = ", temp_change, 
             ", policy = ", policy, 
             ", climate_shock = ", names(climate_shock_scen)[climate_index], "\n"))
}

# define period which policy is implemented 
begin_policy = max(climate_shock) #policy will be implemented the following month
start_policy_month = min(which(end_mon>=begin_policy))+1
policy_length_range = seq(from=(start_policy_month+5), to=120, by=6)

#ensured the whole time period is included 
if (max(policy_length_range)<120) {
  policy_length_range = c(policy_length_range, 120)
}


################################################################################
# Setup model parameters
################################################################################

daily_temperature = rep(read.csv(paste0(inputs_dir, "AvgDailyTempDEBay.csv"))$sst, 10)
endo=1          # prices are endogenous in all scenarios 
temp_dep =1     # running biological model with temperature dependence for natural mortality and growth 


# Runs optimized Effort (optim=1) versus actual effort (optim=0)
for(policy_length in policy_length_range){
  end_policy = end_mon[policy_length]
  
  daily_temp=daily_temperature #resets each run 
  daily_temp[climate_shock] <- daily_temp[climate_shock] + temp_change
  
  shared_env = new.env()
  
  E_0 = c(rep(13000, months), rep(6000, months), rep(8, months)) #initial guess 
  E_0[closed_seasons] = 0 #changes based on different policy scenarios, if season is closed the no fishing (perfect compliance); cuts down on processing time 
  E_0 = E_0[E_0 > 0]
  
  upper_lim = c(rep(Emax_pot, months), rep(fem_pot_max, months), rep(Emax_dredge, months))
  upper_lim[closed_seasons] = 0
  upper_lim = upper_lim[upper_lim > 0]
  
  opts_global = list(
    "algorithm" = "NLOPT_GN_ISRES",  # Global optimizer
    "xtol_rel" = 1e-2,   
    "ftol_rel" = 1e-2,  
    "maxeval" = 50000  
  )
  
  opts_local = list(
    "algorithm" = "NLOPT_LN_COBYLA",  # Supports inequality constraints!
    "xtol_rel" = 1e-4,   
    "ftol_rel" = 1e-4,  
    "maxeval" = 50000
  )
  
  
  for (sim in 1:100) { # run each simulation 100 times
    print(paste("Starting global search for simulation", sim))
    
    res_global = tryCatch({
      nloptr(
        x0 = E_0,
        eval_f = objective_function,
        eval_g_ineq = inequality_constraints,
        lb = rep(0, length(E_0)),
        ub = upper_lim,
        opts = opts_global
      )
    }, error = function(e) {
      print(paste("Global search failed:", e$message))
      return(NULL)
    })
    
    if (!is.null(res_global) && res_global$status > 0) {
      print("Starting local refinement...")
      res_local = tryCatch({
        nloptr(
          x0 = res_global$solution,
          eval_f = objective_function,
          eval_g_ineq = inequality_constraints,
          lb = rep(0, length(E_0)),
          ub = upper_lim,
          opts = opts_local
        )
      }, error = function(e) {
        print(paste("Local search failed:", e$message))
        return(NULL)
      })
    }
    
    final_solution = if (!is.null(res_local) && res_local$status > 0) res_local$solution else res_global$solution
    
    scenario_name = paste0(climate_shock_name,"_", policy,"_", temp_change, "_",policy_length)
    
    write(final_solution, paste0(home, "results/", scenario_name, format(Sys.time(), "%Y-%m-%d_%H-%M-%S"), ".txt"))
  }
  
  
  
}



################################################################################
# initial_population sets up initial population at time 0
# Inputs 
# 1. Total population (number of individuals)
# 2. Distribution is a table that includes 4 columns:
#     a. type : males (M), mature females (F), immature females (I) and sponge (S)
#     b. location: arithmetic mean (mm)
#     c. shape: arithmetic SD (mm)
#     d. prop: percent of the initial pop assigned to each type 
# 
# Output is a data.table of initial carapace width and sex
################################################################################

initial_population = function (cohort_size, male_ratio, mean_male_mm=mean_male_mm, sd_male_mm=sd_male_mm, 
                               mean_female_mm=mean_female_mm, sd_female_mm=sd_female_mm){
  
  all_males = round(male_ratio*cohort_size,2)
  all_females = cohort_size-all_males
  
  initial_pop = data.table(CW = c(rlnorm(all_males, meanlog = mean_male_mm, sdlog = sd_male_mm), 
                                  rlnorm(all_females, meanlog = mean_female_mm, sdlog = sd_female_mm)),
                           type = c(rep("M", all_males), rep("I", all_females)))
  
  initial_pop$sex = ifelse(initial_pop$type == "M", "males", "females")
  
  # Initial population set up 
  initial_pop[,CW:=ifelse(CW>200, 200, CW)]
  initial_pop[,IP:=time_between_molts(CW), by=1:nrow(initial_pop)]
  initial_pop[,degree_days_exposure := 0]
  initial_pop[,rho := peeler_thresh(CW, day=0), by=1:nrow(initial_pop)]
  initial_pop[,terminal_molt:=0]
  initial_pop[,maturity := ifelse(type=="M" & CW>=mature_male_cutoff, 1, 0)]
  initial_pop[,shell_status :=shell_status(molting=0, terminal_molt, degree_days_exposure, IP, rho), by=1:nrow(initial_pop)]
  initial_pop[,spawning_day := ifelse(terminal_molt==1, round(runif(1, min= start_spawn, max = end_spawn)), 0), by=1:nrow(initial_pop)]
  initial_pop[,max_sperm:= ifelse(type=="M", rlnorm(1, meanlog=mean_sperm, sdlog=sd_sperm), 0), by=1:nrow(initial_pop)]
  initial_pop[,max_sperm:= ifelse(max_sperm<700000000 || max_sperm > 6000000000, 1900000000, max_sperm), by=1:nrow(initial_pop)] 
  initial_pop[,male_sperm:= ifelse(maturity==1, max_sperm, 0), by=1:nrow(initial_pop)]
  initial_pop[,data:="model"]
  
  return(initial_pop)
  
}

################################################################################
# time_between_molts is a function which calculates time between molts (IP) for each crab
# IP are calculated at time 0 and at each molt 
# Unit in degree days celsius 
# Modeled to increase exponentially with CW 
# IP is modeled to be drawn from a shifted exponential density function
# 
# Inputs 
# 1. Crab Size (CW)
# 
# Output is a single value for IP
################################################################################
time_between_molts = function(CW){
  
  gamma = 69.70*(1.0149)^(CW)
  beta = (166.39*(1.0115)^(CW))-gamma
  z=runif(1,min=0,max=1) 
  return(gamma-beta*log(1-z))
}

################################################################################
# growth is a function that determines growth per molt 
# 
# Inputs 
# 1. Crab Size (CW)
# 2. Sex (male or female)
# 
# Output is a single value for growth per molt in  mm 
################################################################################
growth = function(CW, sex, temp_dep, daily_temp) {
  if(temp_dep==0){
    # mean_val = ifelse(sex == 2, 1.218 + (7.09 * 10^-4) * CW, 1.25)
    # sd_val = ifelse(sex == 2, 0.07, 0.06)
    # grow = rnorm(length(CW), mean = mean_val, sd = sd_val)
    size = CW/10
    grow = (1 + (0.142 - (0.000113*daily_temp) + (0.0287*size)))
  } else {
    size = CW/10
    grow = (1 + (0.142 - (0.000113*daily_temp) + (0.0287*size)))
  }
  
  return(grow)
}

################################################################################
# peeler_thresh is a function that determines peeler threshold 
# 
# Inputs 
# 1. Crab Size (CW)
# 3. day of year 
# 
# Output is a single value for rho
# rho is the predicted proportion of the IP that is obtained one week before molting
# this vaires during time of year 
################################################################################

peeler_thresh =function(CW, day=doy){
  #october (274-304)
  #rho is highest in oct because individuals carry with them relative high degree-day exposure cfuring winter
  
  rho <- ifelse(day > 120 & day <= 273,  #warm months of may 1 (121) -sept 30 (273) 
                0.60 + ((1.81 * (10^-3)) * CW),
                ifelse(day >= 274 & day <= 304, #Nov through April
                       0.8722 + ((2.94 * (10^-4)) * CW),
                       0.7355 + ((9.59 * (10^-4)) * CW)))
  
  return(rho)
  
}

################################################################################
# shell_status is a function that determines if crab shell status as peeler, soft or 
# hard shell 
# 
# Inputs 
# 1. Crab Size (CW)
# 2. Terminal molt: females who have terminally molted will not grow 
# 3. Degree days exposure to determine soft shell  status 
# 4. rho to determine peeler status  
# 
# Output is a string value for hard or soft shell, or peeler
################################################################################

shell_status = function(molting, terminal_molt, degree_days_exposure, IP, rho) {
  shell = rep(1, length(molting))
  shell[(molting==1)] = 2  # Crabs remain in soft shell for two days
  shell[(terminal_molt==1)] = 1  # Terminally molted crabs stay hard shelled
  shell[(degree_days_exposure > (IP * rho) & degree_days_exposure <= IP)] = 3  # Shell status is peeler
  shell[(degree_days_exposure> IP)] = 2 #crab becomes soft shell 
  
  return(shell)
}

################################################################################
# predict_maturity is a function that calculates the probability that female crabs will 
# mature, after which they stop growing  
#
# Inputs 
# 1. Crab Size (CW) of immature females
# 2. Days of the year they can mature 
# 3. Begin and end days for the mating season
#
# Output is a binary value for whether or not females mature in a given time 0
################################################################################

predict_maturity= function(CW, day, sex, molting, begin = begin_mature, end = end_mature, min_maturity_size = min_maturity_size, temp_dep=temp_dep, terminal_molt, tt=tt, daily_temp, size_mat_females) {
  random_number = runif(length(CW), min = 0, max = 1)
  pr_maturity = numeric(length(CW))  # Preallocate pr_maturity vector
  
  if(day >= begin & day <= end) {
    
    if (temp_dep == 0) {
      pr_maturity[CW >= min_maturity_size & sex == 2 & molting == 1] = 
        0.9994 / (1 + (CW[CW >= min_maturity_size & sex == 2 & molting == 1] / 117.9807)^-28.5056)
      
      mature =as.integer(pr_maturity > random_number)
      mature[which(terminal_molt == 1 | CW >= size_mat_females)] = 1
      
    } else {
      pr_maturity[CW >= min_maturity_size & sex == 2 & molting == 1] = 
        -15.116 + (0.094*CW[ CW >= min_maturity_size & sex == 2 & molting == 1]) + 
        (0.081*salinity) + (0.115*daily_temp)
      
      mature =as.integer(pr_maturity > random_number)
      mature[which(terminal_molt == 1 | CW >= size_mat_females)] = 1
      
    }
  }
  else {
    mature= numeric(length(CW))  
    mature[(which(terminal_molt == 1))] = 1
    mature[(which(CW >= size_mat_females & sex==2))] = 1
    
  }
  
  return(mature)
}

################################################################################
# findspawningday is a function that determines spawning day
# 
# Inputs 
# 1. Crab Size (CW)
# 2. spawning day - each mature female crab is ransombly assigned a spawning day
# 3. day of year - can only spawn from May 15 - Sept 15 
# 
# 
# Spawning potential equaled the sum of the total numbers of eggs predicted to be spawned by 
# females in the second year of the simulation
################################################################################

findspawningday =function(day, yr){
  
  if(yr == 1){
    #spawn day is based on day of calendar year, so crab matures in the first year it can spawn either the first year 
    #or second year based on the random draw
    spawn = round(runif(1, min = start_spawn, max = end_spawn))
    
  } else if(day>=start_spawn & day <=end_spawn & yr>1){
    #if it is year 2 during the spawning period 
    #spawning day can be that day until end of spawning season
    spawn = round(runif(1, min = day, max = end_spawn))
  } else if(day<start_spawn & yr>1){
    #if crab matures in year 2 before spawning season
    spawn = round(runif(1, min = start_spawn, max = end_spawn))
  } else {
    spawn = 0 
  }
  
  return (spawn)
  
}

################################################################################
# natural_mortality is a function that simulates probability of natural mortality
# 
# Inputs 
# 1. shell status = soft shell have twice the probability of natural mortality
# 3. M is yearly natural mortality (daily M is 1-e(^(-M/365)))
# 
# Output is a single binary value for mortality (2= dead, 1=alive)
################################################################################

natural_mortality = function (M_daily, shell_status, soft_mdaily, sex, tt=tt, CW, heat_stress=heat_stress, heat_stress_soft=heat_stress_soft, temp_dep=temp_dep){
  
  random_numbers = runif(length(CW), min = 0, max = 1)
  mort = rep(1, length(CW))  # Initialize mortality vector
  
  if(temp_dep == 0 ) {
    soft = shell_status == 2 
    mort[soft & random_numbers > soft_mdaily] = 2
    female_ = sex == 2
    mort[female_ & random_numbers > M_daily] = 2
    soft_or_female = shell_status == 2 | sex == 2
    other = !soft_or_female
    mort[other & random_numbers > M_daily] = 2
    
    
  } else {
    if(daily_temp[tt]<temp_threshold){ #overwintering months
      h_t = rep(1, length(CW))
      winter_days = 0
      
      for(dd in seq(tt:1)) {
        if(daily_temp[dd] < temp_threshold) {
          winter_days = winter_days +1
        } else {
          break
        }
      } 
      h_t = exp(-((1/0.45)*winter_days^((1/0.45)-1)*exp(-(1/0.45)*(3.59 + (0.1*daily_temp[tt]) + (0.02*salinity) +(0.03*CW)))))
      h_t = M_daily + h_t #overwintering mortality is at least daily natural mortality
      mort[random_numbers > h_t] = 2
      
    } else {
      if (daily_temp[tt]>=max_temp_threshold) {
        # heat stress mortality
        mort[random_numbers > heat_stress]  = 2
        soft = shell_status == 2 
        mort[soft & random_numbers > heat_stress_soft] = 2
        
      } else {
        soft = shell_status == 2 
        mort[soft & random_numbers > soft_mdaily] = 2
        female_ = sex == 2
        mort[female_ & random_numbers > M_daily] = 2
        soft_or_female = shell_status == 2 | sex == 2
        other = !soft_or_female
        mort[other & random_numbers > M_daily] = 2
      }
    }
  } #closes temp dependent loop 
  
  
  return(mort)
  
}

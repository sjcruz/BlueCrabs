################################################################################
# Script sets up all biological and economic variables for use in simulations
################################################################################
options(scipen=999)      #no scientific notation

shared_model = new.env()  # Create an environment to store shared data
home = paste0(getwd (),"/")
fun_dir = paste0(home, "functions/")
inputs_dir = paste0(home,"inputs/")

source (paste0(fun_dir, 'bio_econ_model.R'))
source (paste0(fun_dir, 'biological_functions.R'))
source (paste0(fun_dir, 'economic_policy_functions.R'))
source (paste0(fun_dir, "effort_matrix.R"))
source (paste0(fun_dir, 'obj_function.R'))
source (paste0(fun_dir, 'inequality_constraints.R'))
source (paste0(fun_dir, "recovery_difference_function.R"))


temp_dep =1   # temp dependent scenarios 
endo = 1      # endogenous prices 
years = 10               #Number of years  
days = round(365*years,0)      #Number of days 
months = 12*years              #total months 
N_actual = 10*10^6             # actual stock size (adults) from stock assessment 
first_cohort = 500
total_population = 7000        #this will get updated later as we add new cohorts to the model 
model_ratio = N_actual/1000
proportion_DE_CB = 0.06        #DE harvest is 6% of CB harvest
lbs_conversion = c(0.4, 0.3,0.3, 0.2,0.3)
lbs_conversion = do.call(rbind, replicate(months, lbs_conversion, simplify=FALSE))
max_female_size = 300 # setting up max variable for policy simulations later
Emax_pot = 20000      # max comes from landings data (mx = 16000)
fem_pot_max = 9000    # max number of pots in actual time series for female crabs 
Emax_dredge = 20      # also comes from landings data (max = 14)
Hmax_pot = 100000000  # basically no limit, but can cap for policy scenarios 
Hmax_dredge = 100000000 # basically no limit 
max_pot_catch_per_month = 780000    # average max monthly landings for pots 
max_dredge_catch_per_month = 170000 # average max monthly landings for pots 
dredge_cost = rep(500,12)            # of operating one dredge for a month
pot_cost =  rep(10,12)
fem_pot_cost = pot_cost             #cost of operating one pot for a month 
fem_pot_cost[6:12] = pot_cost[6:12] * 1.5 #increasing cost of female pots in the last half of the year 
cost = cbind(pot_cost, fem_pot_cost, dredge_cost)
cost = do.call(rbind, replicate(years, cost, simplify=FALSE))
a = 9715462     ## from Marty's Land Econ blue crab paper 
choke = c(6, 4, 3, 7, 3)
b = -0.8          ##from Marty's Land Econ blue crab paper 

ii = 3                          #Market categories being optimized (males, females pot and females dredge)
target_fishery = c("large_males", "females_pot", "females_dredge")
cat = 5                         #all market categories 
market_cat = c("large_males", "small_males", "females_pot", "peeler","females_dredge")
life_span = round(2 * 365)      #crabs only live for 2 years 
aug_last = 244 # last day ofaugust after which recuirtment happens
cohort_ends = c(life_span, aug_last + life_span, aug_last + life_span + (365 * (1:years)))
days_cal = rep(c(31,28,31,30,31,30,31,31,30,31,30,31), years)
end_mon = c(cumsum(days_cal))
start_mon = c(1,c(end_mon[1:length(end_mon)-1])+1)

################################################################################
# Biological variables 
sex_ratio = 0.5
temp=read_excel(paste0(inputs_dir, "AvgDailyTempDEBay.xlsx"))
salinity = 18 
mean_male_mm = 3.67             #using arithmetic mean CW for log normal dist from April trawl data 2000-2022 for upper bay
mean_female_mm = 3.61
sd_male_mm = 0.439
sd_female_mm = 0.425
recruit_males = 3.49            ##recruit size distribution in september from 2000-2022, crabs less than 60mm in length condisered new recruits 
recruit_females = 3.49
recruit_m_sd = 0.315
recruit_f_sd = 0.307
mature_male_cutoff = 107        ##Size cut off for mature males 
temp_threshold = 8.9            ##Temp in degrees C after which crabs do not grow 
max_temp_threshold = 33         ## Mortality rates increase after temps increase above this threshold 
M_daily = exp(-0.8/365)            ##daily probability of natural mortality 
soft_mdaily = exp(-(0.8*2)/365)        ##soft shells have 2X thfirst_cohorte natural mortality rate 
heat_stress = exp(-(0.8)/45)
heat_stress_soft = exp(-(0.8*2)/45)
start_spawn = 135               ##May 15 beginning of spawning season 
end_spawn = 258                 ##Sept 15 end of spawning season, used to evaluate spawning potential
begin_mature = 91               ##Females can mature April 1- October 1 
end_mature = 274                ##
min_maturity_size = 80         ##setting minimum maturity size for females based on trawl data 
size_mat_females = 120          ##female maturity =1 if they hit 120mm without maturing, this number comes from Rich's trawl survey analysis. at 120mm 96% of females are mature
softshell_size = 89             ##soft shell size limit (3.5in = 89 mm)
peeler_size = 76                ##peeler size limit (3in = 76mm)
sperm_thresh = 300000000        ##If males crab have less than threshold then they are not able to mate
largemale_cutoff = 140          ##Market category cutoff, using 140mm for DE crab that are smaller than CB
back_mean = 3012612613          ##Using back calculated numbers from Taylor sent May 22,2023 for male sperm distribution
back_sd = 1520905987            ##Using back calculated numbers from Taylor sent May 22,2023 for male sperm distribution
mean_sperm = log(back_mean^2/sqrt(back_sd^2 + back_mean^2))    ##Calculate arithmetic log mean for log normal 
sd_sperm = sqrt(log(1+(back_sd^2/back_mean^2)))                ##Numbers are very similar to Kendall et al 2021 (log mean = 21.49 and logsd = 0.56)
percent_viability = 0.8         #ratio of live:dead sperms- used as a 3rd reduction in sperm at mating, from Taylors work - ideally we have this at the monthly level to capture that variation
max_crab_size = 200             ##capping the size of crabs to the largest recorded crab in the trawl survey
days_bw_broods = 21             ##Number of days between brood production, she will remain sponge for this period of time (Dickinson et al. 2006)
recruit_alpha = 206.095         #alpha for recruitment ricker function (Stock assessment)
recruit_beta =-40.759           #beta for recruitment ricker function (Stock assessment)

################################################################################
#Economic variables for optimization
alpha = rep(1,ii)                ##effort elasticity 
beta = rep(1,ii)                 ##stock elasticity
pot_fem_proportion = rep(0.6,12)  # only 20% of total fishable females are available to the pot and trap fishery, this captures the movement of females to the lower bay
pot_fem_proportion[6:12] = 0.15
pot_fem_proportion[1:2] = 0.15
pot_fem_proportion = rep(pot_fem_proportion, years)

### Price data for now using average monthly prices for each category to capture seasonality 2015-2019
prices = fread(paste0(inputs_dir, "prices.csv"), drop = "V1")[,2:6]
prices = do.call(rbind, replicate(years, prices, simplify=FALSE))

price_ratios = prices 
price_ratios[,2] = price_ratios[,2]/price_ratios[,1]
price_ratios[,3] = price_ratios[,3]/price_ratios[,1]
price_ratios[,4] = price_ratios[,4]/price_ratios[,1]
price_ratios[,5] = 0.5
price_ratios[,1] = 1

################################################################################
# Set up  matrices to track variables 
N_0 = initial_population(cohort_size=first_cohort, male_ratio=sex_ratio, mean_male_mm=mean_male_mm, 
                         sd_male_mm=sd_male_mm, mean_female_mm=mean_female_mm, sd_female_mm=sd_female_mm)
shell = matrix(0,nrow=first_cohort, ncol=days+1)
CW = matrix(0, nrow=first_cohort, ncol=days+1)
growth_per_molt = rep(1, first_cohort)
degree_days_exposure =  rep(0, first_cohort)
IP = rep(0, first_cohort)
terminal_molt = rep(0, first_cohort)
sponge =  rep(0, first_cohort)
rho = rep(0, first_cohort)
survived = matrix(1,nrow=first_cohort, ncol=days+1)
molting =  matrix(0,nrow=first_cohort, ncol=days+1)
size = matrix(0, ncol = 5)

male1_fishery = matrix (0,nrow=total_population, ncol = months)
male2_fishery = matrix (0,nrow=total_population, ncol = months)
female_fishery = matrix (0,nrow=total_population, ncol = months)
peeler_fishery = matrix (0,nrow=total_population, ncol = months)

spawning_day =  matrix(0, nrow=total_population, ncol=days+1)
spawning_potential = matrix(0,nrow=total_population, ncol=days+1)
male_sperm = rep(0, total_population)
female_sperm = matrix(0,nrow=total_population, ncol=days+1)
female_sperm_at_mating = rep(0, total_population)
mated_males =  matrix(0,nrow=total_population, ncol=days+1)
broods =  matrix(0,nrow=total_population, ncol=days+1)
sponge =  matrix(0,nrow=total_population, ncol=days+1)
mated_female = rep(0, total_population)
mature_males =  rep(0, total_population)
spawning_day[,1] = N_0$spawning_day
max_sperm = N_0$max_sperm
mated_female = rep(0, total_population)

X = matrix(0, months, cat)
harvest = matrix (0,months,cat)  
pounds = matrix (0,months,cat) 
male_fish = matrix (0,months,cat)  
female_fish = matrix (0,months,cat) 
dred_fish  = matrix (0,months,cat) 
population =  matrix (0,months,11)  
################################################################################
# Biological Model 
shell[,1] = N_0$shell_status
CW[,1] = N_0$CW
degree_days_exposure = N_0$degree_days_exposure
IP = N_0$IP
terminal_molt=N_0$terminal_molt
#mature_males=N_0$maturity
rho=N_0$rho
survived[,1] = 1 #2 or 3 = died due to natural/fishing mort, 1 = survived
molting[,1] = 0 
sex=ifelse(N_0$sex=="males",1,2)
rec_by_cohort = matrix (0,nrow=years+1, ncol = 2)
rec_by_cohort[1,1] = first_cohort


day_and_year = function(day_of_year) {
  year = ceiling(day_of_year / 365)  # Calculate the year
  day_in_year = ifelse(day_of_year %% 365 == 0, 365, day_of_year %% 365)  # Calculate the day within the year
  
  return(list(year = year, day_in_year = day_in_year))
}





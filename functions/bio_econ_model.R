
################################################################################
# biological_model is a function that tracks all individual crabs through their
# life cycle using biological and economic functions
################################################################################

biological_model = function(N_0=N_0, effort=effort){
  ################################################################################
  ## time loop starts here 
  #effort = summary_E(scenario = "actual", type="mean")
  #profvis({
  pop = first_cohort
  
  for (tt in 2:(days)){
    #setting day of the year and year variable 
    doy = day_and_year(tt)$day_in_year
    year = day_and_year(tt)$year
    
    ################### ACCUMULATION OF DEGREE DAYS ##############################
    exposure = daily_temp[doy]-temp_threshold 
    exposure = ifelse(exposure<0, 0, exposure) # temp is less than threshold then no degree days are accumulated 
    exposure = ifelse(daily_temp[doy]>=max_temp_threshold, 0, exposure)  # max temp threshold 
    
    # Crab that are dead, terminally molted, or are already at max size do not accumulate degree days 
    condition_1 = survived[1:pop, tt-1] == 1 & sex[1:pop] == 1 & CW[1:pop, tt-1] < max_crab_size
    condition_2 = survived[1:pop, tt-1] == 1 & sex[1:pop] == 2 & terminal_molt[1:pop] == 0 & CW[1:pop, tt-1] < size_mat_females
    degree_days_exposure[condition_1] = degree_days_exposure[condition_1] + exposure
    degree_days_exposure[condition_2] = degree_days_exposure[condition_2] + exposure
    
    ################### SHELL STATUS AND MOLTING ##############################
    shell[1:pop, tt] = shell_status(molting=molting[1:pop,tt-1], terminal_molt=terminal_molt[1:pop], degree_days_exposure=degree_days_exposure[1:pop], 
                                    IP=IP[1:pop], rho=rho[1:pop])
    
    shell[(survived[1:pop,tt-1]!=1),tt] = 0 
    
    #if crabs become soft shell, and this is their first day in soft shell) and are smaller than the max size then they molt 
    molting[1:pop,tt] = ifelse(shell[1:pop, tt]==2 & molting[1:pop,tt-1]==0, 1,0) 
    #if they molt they grow, calculate new IP, rho and reset degree days
    
    growth_per_molt[1:pop] = growth(CW=CW[1:pop,tt-1], sex=sex[1:pop], daily_temp=daily_temp[tt], temp_dep = temp_dep)
    growth_per_molt[(molting[1:pop,tt]!=1)] =1 
    
    CW[1:pop,tt] =  CW[1:pop,tt-1]*growth_per_molt[1:pop] 
    
    CW[(sex[1:pop]==1 & CW[1:pop,tt]>=max_crab_size),tt] = max_crab_size #max size is capped for females and males 
    CW[(sex[1:pop]==2 & CW[1:pop,tt]>=size_mat_females),tt] = size_mat_females #
    
    idx_molting = which(molting[1:pop,tt]==1)
    
    IP[idx_molting] =  time_between_molts(CW=CW[idx_molting,tt])
    rho[idx_molting] =  peeler_thresh(CW[1:pop,tt], day=doy)
    degree_days_exposure[idx_molting] = 0
    
    ############################ MATURITY  ##############################
    #probability of females maturing depends on CW and time of year
    terminal_molt[1:pop] = predict_maturity(CW=CW[1:pop,tt], day=doy, sex=sex[1:pop], molting=molting[1:pop,tt], 
                                            begin=begin_mature, end=end_mature, min_maturity_size = min_maturity_size, 
                                            temp_dep=temp_dep, terminal_molt=terminal_molt[1:pop], daily_temp=daily_temp[tt], size_mat_females=size_mat_females)
    #after terminal molt they no longer accumulate degree days exposure and thus do not continue to grow 
    degree_days_exposure[(terminal_molt[1:pop]==1)]  =0
    
    #Males mature: size cutoff
    mature_males[which(sex[1:pop] == 1 & CW[1:pop,tt] >= mature_male_cutoff & survived[1:pop,tt-1]==1)] = 1
    #if they mature in this time step then they get assigned some max sperm
    male_sperm[which(mature_males[1:pop]==1 & male_sperm[1:pop]==0)] = max_sperm[which(mature_males[1:pop]==1 & male_sperm[1:pop]==0)]

    ##############  Rains et al RANDOM MATING ##################################
    #mating constrained to season when females terminally molt
    female_mates = which(terminal_molt[1:pop] == 1 & mated_female[1:pop] == 0 & survived[1:pop,tt-1]==1 & sex[1:pop]==2)
    male_mates = which(mature_males[1:pop]==1 & shell[1:pop, tt]== 1 & male_sperm[1:pop]> sperm_thresh & mated_males[1:pop,tt-1]==0 & survived[1:pop,tt-1]==1)

    if (length(male_mates) > 0 && length(female_mates) > 0) {
      for(fem in female_mates){ # loop through every female to find her a mate
        mate = sample(male_mates, 1) # randomly picks a male to mate with
        mated_males[mate,tt] = 1
        mated_female[fem]= 1
        female_sperm[fem,tt] = 0.5*0.5*male_sperm[mate]*percent_viability # she gets some sperm at mating
        female_sperm_at_mating[fem]= female_sperm[fem,tt]
        male_sperm[mate] = 0.5 * male_sperm[mate] # male sperm is reduced by 50%

        #Females become sponge females after 21 dyas of mating
        spawning_day[female_mates, tt] = ifelse( (doy + days_bw_broods) <= end_spawn & (doy + days_bw_broods) >= start_spawn,
                                                 doy + days_bw_broods,
                                                 start_spawn)
      }}

    non_mating_males = which(mated_males[1:pop,tt]==0 & sex[1:pop]==1 & male_sperm[1:pop]<max_sperm[1:pop])
    male_sperm[non_mating_males] = male_sperm[non_mating_males]*exp(0.057)

    #################### SPAWNING #############################################
    # Currently not size dependent, assuming an egg:sperm ratio of 1:4
    eligible_females = which(survived[1:pop,tt-1] == 1 & female_sperm[1:pop, tt] > 0 & spawning_day[1:pop, tt] == doy)
    #spawning_potential[ind,tt]= (-2.248 + (0.337*CW[ind,tt]))*100000
    spawning_potential[eligible_females,tt]= 3000000
    female_sperm[eligible_females, tt] = female_sperm[eligible_females, tt] - ( spawning_potential[eligible_females,tt] * 4)
    # if she has enough sperm then she makes one brood
    broods[eligible_females, tt] = ifelse(female_sperm[eligible_females, tt] >= 0, 1, 0)
    # if she has enough sperm and spawns she gets assigned another spawning day
    spawning_day[eligible_females & broods[eligible_females, tt] == 1, tt] =  ifelse( (doy + days_bw_broods) <= end_spawn & (doy + days_bw_broods) >= start_spawn,
                                                                                      doy + days_bw_broods,
                                                                                      start_spawn)

    spawning_day[eligible_females & broods[eligible_females, tt] == 0, tt] = 0

    #is she spawns she remains a sponge crab until next spawning day, ensures she cannot be fished
    sponge[eligible_females, doy:(doy + days_bw_broods)] = ifelse(broods[eligible_females, tt] == 1, 1, 0)
    
    #################### NATURAL MORTALITY ######################################
    survived[1:pop,tt] = ifelse(survived[1:pop,tt-1]==1, 
                                natural_mortality(M_daily=M_daily, shell_status=shell[1:pop,tt], soft_mdaily=soft_mdaily, sex=sex[1:pop], tt=tt, CW=CW[1:pop, tt],
                                                  heat_stress=heat_stress, heat_stress_soft=heat_stress_soft, temp_dep=temp_dep), survived[1:pop,tt-1])
    
    
    CW[which(survived[1:pop, tt] == 2), tt] = 0
    
    #################### FISHING MORTALITY ######################################
    
    if(tt%in%end_mon){ #end mon is the last calendar day of each month for 
      month = which(tt==end_mon)
      # crabs vulnerable to fishing in each market category
      
      male1_fishery[1:pop,month] = ifelse(survived[1:pop,tt]==1 & CW[1:pop,tt]>=largemale_cutoff & shell[1:pop,tt]==1 & sex[1:pop]==1, 1, 0)
      male2_fishery[1:pop,month] = ifelse(survived[1:pop,tt]==1 & CW[1:pop,tt]>= hardshell_size[month] &  CW[1:pop,tt]<largemale_cutoff & shell[1:pop,tt]==1 & sex[1:pop]==1, 1, 0)
      female_fishery[1:pop,month] = ifelse(survived[1:pop,tt]==1 & shell[1:pop,tt]==1 & sex[1:pop]==2 & terminal_molt[1:pop]==1 & sponge[1:pop]==0 & CW[1:pop,tt]<= max_female_size,1,0) #mature hard shell females no size limit
      peeler_fishery[1:pop,month] = ifelse((survived[1:pop,tt]==1 & shell[1:pop,tt] == 3 & CW[1:pop,tt]>=peeler_size) | # crabs in peeler stage can be fished 
                                             #(survived[1:pop,tt]==1 & shell[1:pop,tt]==2 & sex[1:pop]==2 & terminal_molt[1:pop]==1) | #females right after terminal molt can be fished, soft shell state
                                             (survived[1:pop,tt]==1 & shell[1:pop,tt]==2 & CW[1:pop,tt]>= softshell_size), 1, 0) #crabs in soft shell stage can be harvested
      
      total_females = sum(female_fishery[1:pop,month])
      X[month, ] = c(sum(male1_fishery[1:pop,month]), sum(male2_fishery[1:pop,month]), total_females*pot_fem_proportion[month], sum(peeler_fishery[1:pop,month]),total_females*(1-pot_fem_proportion[month])) 
      
      population[month, ] =c(sum(survived[1:pop,tt]==1 & sex[1:pop]==1), # crabs that are alive
                             sum(survived[1:pop,tt]==2 & sex[1:pop]==1), #natural mortality
                             sum(survived[1:pop,tt]==1 & CW[1:pop,tt]>=largemale_cutoff & shell[1:pop,tt]==1 & sex[1:pop]==1), #large males, hard shell 
                             sum(survived[1:pop,tt]==1 & CW[1:pop,tt]<largemale_cutoff & shell[1:pop,tt]==1 & sex[1:pop]==1), #small males, hard shell
                             sum(survived[1:pop,tt]==1 & shell[1:pop,tt]>1 & sex[1:pop]==1), #peeler and soft shell males 
                             
                             sum(survived[1:pop,tt]==1 & sex[1:pop]==2), #female crabs that are alive
                             sum(survived[1:pop,tt]==2 & sex[1:pop]==2), #female natural mortality
                             sum(survived[1:pop,tt]==1 & shell[1:pop,tt]==1 & sex[1:pop]==2 & terminal_molt[1:pop]==1), #mature females
                             sum(survived[1:pop,tt]==1 & shell[1:pop,tt]==1 & sex[1:pop]==2 & terminal_molt[1:pop]==0), #immature females 
                             sum(survived[1:pop,tt]==1 & sex[1:pop]==2 & sponge[1:pop]==1), #sponge females
                             sum(survived[1:pop,tt]==1 & shell[1:pop,tt]>1 & sex[1:pop]==1)) #soft shell/peeler females
      
      ####### Imperfect selectivity 
      #calculate harvest and bycatch of pot male fishery 
      male_fish[month,] = c(q_pot_males[month, 1]*(effort[month,1]^alpha[1])*(X[month,1]^beta[1]), #direct harvest of large males 
                            q_pot_males[month, 2]*(effort[month,1]^alpha[1])*(X[month,2]^beta[1]), #bycatch of small males 
                            q_pot_males[month, 3]*(effort[month,1]^alpha[1])*(X[month,3]^beta[1]), #bycatch of females
                            q_pot_males[month, 4]*(effort[month,1]^alpha[1])*(X[month,4]^beta[1]), #bycatch of peelers
                            q_pot_males[month, 5]*(effort[month,1]^alpha[1])*(X[month,5]^beta[1])) #bycatch of  female dredge = 0 
      
      #calculate harvest and bycatch of pot female fishery 
      female_fish[month,] = c(q_pot_females[month, 1]*(effort[month,2]^alpha[2])*(X[month,1]^beta[2]), #bycatch large males 
                              q_pot_females[month, 2]*(effort[month,2]^alpha[2])*(X[month,2]^beta[2]), #bycatch of small males 
                              q_pot_females[month, 3]*(effort[month,2]^alpha[2])*(X[month,3]^beta[2]), #direct harvest of females
                              q_pot_females[month, 4]*(effort[month,2]^alpha[2])*(X[month,4]^beta[2]), #bycatch of peelers 
                              q_pot_females[month, 5]*(effort[month,2]^alpha[2])*(X[month,5]^beta[2])) #bycatch of  female dredge = 0 
      
      #calculate harvest and bycatch of dredge female fishery 
      dred_fish[month,] = c(q_dredge[month, 1]*(effort[month,3]^alpha[3])*(X[month,1]^beta[3]), #bycatch large males 
                            q_dredge[month, 2]*(effort[month,3]^alpha[3])*(X[month,2]^beta[3]), #bycatch of small males 
                            q_dredge[month, 3]*(effort[month,3]^alpha[3])*(X[month,3]^beta[3]), #bycatch of females pots = 0 
                            q_dredge[month, 4]*(effort[month,3]^alpha[3])*(X[month,4]^beta[3]), #bycatch of peelers = 0
                            q_dredge[month, 5]*(effort[month,3]^alpha[3])*(X[month,5]^beta[3])) #direct harvest of females using dredge 
      
      #sums each market category to get total harvest 
      harvest[month, ] = as.numeric(male_fish[month,]) + as.numeric(female_fish[month,]) + as.numeric(dred_fish[month,])
      #harvest[month, ] = ifelse(harvest[month, ] > X[month, ], X[month, ], harvest[month, ])
      
      harvest[month, ] = ifelse(harvest[month, ]<0, 0, harvest[month, ])
      harvest[month,] = ifelse(harvest[month,]>X[month,], X[month,], harvest[month,])
      #Randomly assign fishing mortality based on Number of crabs harvested 
      HH = round(harvest, 0)
      
      
      total_females_harvested = harvest[month,3] + harvest[month,5] 
      #if harvest is more than X then just fish all
      rand_males1 = ifelse(HH[month, 1]>=X[month,1], which(male1_fishery[,month]==1), sample(x=which(male1_fishery[,month]==1), HH[month,1], replace = F))
      rand_males2 = ifelse(HH[month, 2]>=X[month,2], which(male2_fishery[,month]==1), sample(x=which(male2_fishery[,month]==1), HH[month,2], replace = F))
      rand_females = ifelse(HH[month,3]>=X[month,3], which(female_fishery[,month]==1), sample(x=which(female_fishery[,month]==1), HH[month,3], replace = F))
      rand_peeler = ifelse(HH[month, 4]>=X[month,4], which(peeler_fishery[,month]==1), sample(x=which(peeler_fishery[,month]==1), HH[month,4], replace = F))
      dredgefished = ifelse(HH[month,5]>=X[month,5], which(!female_fishery[,month]==1%in%rand_females), sample(x=which(!female_fishery[,month]==1%in%rand_females), HH[month,5], replace = F))
      
      survived[c(rand_males1, rand_males2, rand_females, rand_peeler),tt]= 3
      survived[dredgefished,tt] = 4 
      
      pounds[month, ] = c(sum((0.000631 * (CW[rand_males1,tt]^2.47))), 
                          sum((0.000631 * (CW[rand_males2,tt]^2.47))),
                          sum((0.000631 * (CW[rand_females,tt]^2.47))),
                          sum((0.000631 * (CW[rand_peeler,tt]^2.47))),
                          sum((0.000631 * (CW[dredgefished,tt]^2.47))))
      
      CW[c(rand_males1, rand_males2, rand_females, rand_peeler, dredgefished),tt]= 0
      
      size_dis = cbind(month = rep(month, pop), CW = CW[1:pop,tt], survived = survived[1:pop,tt], sex = sex[1:pop], shell = shell[1:pop,tt])
      size = rbind(size, size_dis)
      
    } #close fishing at the end of the month 
    
    
    ################  ends fishing loop 
    
    ############## RECRUITEMNT #################################################
    # Cohorts recruit on last day of July, to begin simulation on Aug 1st
    for(coh in 1:10) {
      if(tt == c(aug_last+(365*(coh-1)))){
        # Density dependent recruitment - Ricker function (Rich -Stock assessment)
        ssb_cw = CW[(survived[1:pop,tt] == 1 & terminal_molt[1:pop]==1), tt]
        mat_females = length(ssb_cw)
        #print(mat_females)
        ssb_kg = sum(0.000631 * (ssb_cw^2.47))/1000000
        rec_per_spawners = recruit_alpha*ssb_kg*exp(recruit_beta*ssb_kg)
        #print(rec_per_spawners)
        cohort = round(((mat_females* rec_per_spawners))/2)*2
        #print(cohort)
        #cohort = round(((sum(spawning_potential[1:pop,1:aug_last])*per)/100)/2)*2
        
        if(cohort>0){
          recruits = initial_population(cohort_size=cohort, male_ratio=sex_ratio, mean_male_mm=recruit_males, sd_male_mm=recruit_m_sd, 
                                        mean_female_mm=recruit_females, sd_female_mm=recruit_f_sd)
          
          shell = rbind(shell, matrix(0, nrow = cohort, ncol = days + 1))
          CW = rbind(CW, matrix(0, nrow = cohort, ncol = days + 1))
          survived = rbind(survived, matrix(0, nrow = cohort, ncol = days + 1))
          molting = rbind(molting, matrix(0, nrow = cohort, ncol = days + 1))
          
          degree_days_exposure = c(degree_days_exposure, rep(0, cohort))
          IP = c(IP, rep(0, cohort))
          rho = c(rho, rep(0, cohort))
          terminal_molt = c(terminal_molt, rep(0, cohort))
          sponge = c(sponge, rep(0, cohort))
          sex = c(sex, rep(NA, cohort))
          
          idx = (pop+1):(pop+cohort)
          
          shell[idx,tt] = recruits$shell_status
          CW[idx,tt]  = recruits$CW
          degree_days_exposure[idx] = recruits$degree_days_exposure
          IP[(pop+1):(pop+cohort)] = recruits$IP
          terminal_molt[idx]=recruits$terminal_molt
          rho[idx]=recruits$rho
          survived[idx,tt] = 1 #2 or 3 = died due to natural/fishing mort, 1 = survived
          molting[idx,tt] = 0 
          sex[idx]= ifelse(recruits$sex=="males",1,2)
          max_sperm[(pop+1):(pop+cohort)] = recruits$max_sperm
          male_sperm[(pop+1):(pop+cohort)] = recruits$male_sperm
          mature_males[(pop+1):(pop+cohort)]=recruits$maturity
          mated_female[(pop+1):(pop+cohort)] = rep(0, cohort)
          spawning_day[(pop+1):(pop+cohort),tt] = recruits$spawning_day
          
        } else {
          next
        }
        pop = pop + cohort
        rec_by_cohort[coh+1,1] =  cohort
        rec_by_cohort[coh+1,2] =  mat_females
      }
    }
    
    ############## Cohorts only live for two years
    for (ccc in seq_along(cohort_ends)) {
      if (tt == cohort_ends[ccc] && rec_by_cohort[ccc, 1] > 0) {
        if (ccc == 1) {
          idx = 1:first_cohort
        } else {
          start = sum(rec_by_cohort[1:(ccc - 1), 1]) + 1
          end   = sum(rec_by_cohort[1:ccc, 1])
          idx = start:end
        }
        
        CW[idx, tt] = 0
        survived[idx, tt] = ifelse(survived[idx, tt] == 1, 2, survived[idx, tt])
      }
    }
    
    ############## 
    
  } # ends daily loop 
  
  # })
  return(list(fishable = X, recruits = rec_by_cohort, stock = population, pounds = pounds, harvest = harvest, male_pot = male_fish, female_pot = female_fish, dred_fish = dred_fish, 
              size = size ))
}





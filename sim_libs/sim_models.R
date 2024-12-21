library(dplyr) ## Needed for group by hh_size
log_lik<-function(simulated,real){
  if(length(simulated)>1){
    result_<-rep(0,length(simulated))
    result_[simulated>0]<- log(simulated[simulated>0])*real[simulated>0]-simulated[simulated>0]
      return(result_)
  }else{
    if(simulated>0){result_<-  log(simulated)*real- simulated}else{result_<-0}
    return(result_)
  }
}

minusLS<-function(simulated,real){
  result_<- -1* (simulated-real)**2
  return(result_)
}
# 
# target<-2.5
# generate_gauss<-function(this_mean){
#   dum<-mean(rnorm(100, this_mean, 0.2))
#   return(list("Score"=minusLS(target,dum)))
# }
# 
# 
# optObj <- bayesOpt(
#   #FUN=SEIR_simulate_two_networks_forBOptim_Binloglik_Ncases_by_HH_averageN
#   FUN=generate_gauss,
#   initPoints = 100
#   , bounds = list(this_mean=c(0,5))
#   , iters.n = 100
# )


SEIR_simulate_var_comm_trans <-function(HH_network,beta1, beta2){
  ### Read the comm incidence
  df_comm_transm<-read.csv("../epi_fitting/community_transmission.csv")
  sigmaz<-1./1.
  gamma<-1./(7)
  ## but there is the MINIMUM of 1 day of infectiousness
  df_epi_info<-read.csv(file='../epi_fitting/epidemic_info/epidemic_info_net_onset_yest.csv')  
  df_epi_info$cdc_primary_onsetdt<-as.Date(df_epi_info$cdc_primary_onsetdt)  ## Define as date
  df_epi_info$delay_isolation<-df_epi_info$delay_isolation#+ceiling(1./sigmaz) ## add sigma to time
  ## (transition I->R done on infected at t-1)
  
  N_timesteps<-max(df_epi_info$max_t_simulation)  ## NB it is 42
  ### This is the maximum of variable "max_t_simulation"
  ## which is defined as 28 days (follow up) + time from cdc_primary_onset (time of symptom onset for index) to doe_dt (enrollment)
  ## So 28 + doe_dt - cdc_primary_onset
  #### NB: later, only infections happening before max_t_simulation are considered
  nrVertices <- network.size(HH_network)
  vertexNames <- HH_network %v% "vertex.names"  ## Needed to check isolation status
  infected <- rep(0, nrVertices)
  susceptible <- rep(1, nrVertices)
  ids_to_infec<-which(network::get.vertex.attribute(x = HH_network,"index")==1)
  infected[ids_to_infec] <- 1		
  susceptible[ids_to_infec] <- 0		
  infectedMat <- infected
  exposed <- rep(0, nrVertices)
  incidenceMat <- infected
  susceptibleMat<-susceptible
  # Individual attributes
  vertex_id <- 1:nrVertices
  timeSteps <- 1
  while(timeSteps<N_timesteps){ 
    hhfoi <- rep(0, nrVertices)
    # Expose each susceptible with a certain probability
    incidence <- rep(0, nrVertices)
    exposed_incidence <- rep(0, nrVertices)
    infectedPrev <- infected	# indicator of infection in the previous time step
    exposedPrev  <- exposed   # indicator of exposure in the previous time step
    susceptiblePrev<-susceptible
    # Transition from S---> E
    for(i in vertex_id[susceptible == 1]){
      vid<-vertexNames[i]
      t_isolation_this_vid<-df_epi_info$delay_isolation[which(df_epi_info$memid==vid)]
      #print(paste0("t=", timeSteps, "i= ",i," exposed= ", exposed[i]," susceptible= ",susceptible[i]))
      
      hh_contacts <- get.neighborhood(HH_network, i)
      nrInfHH <- sum(infectedPrev[vertex_id %in% hh_contacts ], na.rm = TRUE)
      ### Match household to timing and site
      absolute_time<-df_epi_info$cdc_primary_onsetdt[df_epi_info$memid==vid]+timeSteps
      site_of_id<-df_epi_info$site[df_epi_info$memid==vid]
      date_to_next_wednesday<-lubridate::ceiling_date(absolute_time,unit="weeks",week_start=3)
      ### debug: lubridate::wday(lubridate::ceiling_date(absolute_time,unit="weeks",week_start=3),label=TRUE)
      ### Get the relevant factor
      selection<-df_comm_transm$location==site_of_id& df_comm_transm$date==date_to_next_wednesday
      commModulationFactor<-df_comm_transm$ratio_weekly_cases[selection]
      #      print(paste0("For participant ",vid, " after ",timeSteps, " timesteps the date is ", absolute_time, " which is ceiled to ", date_to_next_wednesday, " for ", commModulationFactor ))
      infProb <- 1-(1-beta1)^nrInfHH*(1-beta2*commModulationFactor)
      #comfoi[i] <- beta2
      SuccessExposure<-runif(1) < infProb
      if(SuccessExposure){
        exposed[i]<-1
        susceptible[i]<-0
        #   print(paste0("At time=",timeSteps , " id= ",i, "is removed from susceptibles"))
      }
      exposed_incidence[i]<-exposed[i]
    }
    susceptibleMat<-cbind(susceptibleMat, susceptible)
    # Transition from E---> I (TBD on exposed at previous time step)
    for(i in vertex_id[exposedPrev == 1]){
      SuccessInfection<-runif(1) < sigmaz
      if(SuccessInfection){
        infected[i]<-1
        exposed[i]<-0
        incidence[i] <- 1
      }
    }
    # Transition from I---> R (TBD on infected at previous time step)
    for(i in vertex_id[infectedPrev == 1]){
      SuccessRecover<-runif(1) < gamma
      if(SuccessRecover){
        #infected[i]<-NA ## NAs mark who has eventually been infected
        infected[i]<-0 ## NAs mark who has eventually been infected
      }
    }
    # Append current status to matrix and update time variable
    infectedMat <- cbind(infectedMat, infected)
    incidenceMat <- cbind(incidenceMat, incidence)
    timeSteps <- timeSteps + 1    
  }
  
  ## Ignoring infections happening after the max_time in the study
  infectedMat_censored<-infectedMat
  for(irow in 1:nrow(infectedMat_censored)){
    censoring_time_this_part<-df_epi_info$max_t_simulation[df_epi_info$memid==vertexNames[irow]]
    infectedMat_censored[irow,censoring_time_this_part:ncol(infectedMat_censored)]<-0
  }
  ever_infected<-base::rowSums(infectedMat_censored,na.rm = TRUE)
  ever_infected[ever_infected>0]<-1
  df_final_infection<-data.frame( "memid"=network::get.vertex.attribute(HH_network,"vertex.names"),
                                  "index"=network::get.vertex.attribute(HH_network,"index"),
                                  "ever_infected"=ever_infected)
  return(df_final_infection) 
}

get_int_vertex_attribute_list_network<-function(x,attrname){
  Nvertices<-sum(unlist(lapply(x, network.size)))
  vec_attr<-NULL
  for(i in 1:length(x)){
    dum_attr<-network::get.vertex.attribute(x[[i]],attrname)
    vec_attr<-append(vec_attr,dum_attr) 
  }
  return(as.integer(vec_attr))
  
}


###############################################################################
#### Minimal working example                           ####
###############################################################################
#### Simulation only functions 
SEIR_simulate <- function(HH_network,beta1, beta2){
  sigmaz<-1./1.
  gamma<-1./(7)    ### Average duration_pos is 9.796 (computed from swabs)
  ## but there is the MINIMUM of 1 day of infectiousness
  df_epi_info<-read.csv(file='../epi_fitting/epidemic_info/epidemic_info_net_onset_yest.csv')  
  ## (transition I->R done on infected at t-1)
  N_timesteps<-max(df_epi_info$max_t_simulation)  ## NB it is 42
  ## (transition I->R done on infected at t-1)
  nrVertices <- network.size(HH_network)
  infected <- rep(0, nrVertices)
  susceptible <- rep(1, nrVertices)
  ids_to_infec<-which(network::get.vertex.attribute(x = HH_network,"index")==1)
  infected[ids_to_infec] <- 1		
  susceptible[ids_to_infec] <- 0		
  infectedMat <- infected
  exposed <- rep(0, nrVertices)
  incidenceMat <- infected
  # hhfoiMat <- rep(0, nrVertices)
  # comfoiMat <- rep(0, nrVertices)
  susceptibleMat<-susceptible
  # exposedMat <- exposed
  # exposed_incidenceMat <- exposed
  # Individual attributes
  hh_id <- HH_network %v% "hh_id"
  vertexNames <- HH_network %v% "vertex.names"
  vertex_id <- 1:nrVertices
  timeSteps <- 1
  while(timeSteps<N_timesteps){ 
    # Repeat loop for 29 time steps
    hhfoi <- rep(0, nrVertices)
    #comfoi<- rep(0, nrVertices)      
    # Expose each susceptible with a certain probability
    incidence <- rep(0, nrVertices)
    exposed_incidence <- rep(0, nrVertices)
    infectedPrev <- infected	# indicator of infection in the previous time step
    exposedPrev  <- exposed   # indicator of exposure in the previous time step
    susceptiblePrev<-susceptible
    # Transition from S---> E
    for(i in vertex_id[susceptible == 1]){
      #print(paste0("t=", timeSteps, "i= ",i," exposed= ", exposed[i]," susceptible= ",susceptible[i]))
      hh_contacts <- get.neighborhood(HH_network, i)
      nrInfHH <- sum(infectedPrev[vertex_id %in% hh_contacts ], na.rm = TRUE)
      infProb <- 1-(1-beta1)^nrInfHH*(1-beta2)
      #comfoi[i] <- beta2
      SuccessExposure<-runif(1) < infProb
      if(SuccessExposure){
        exposed[i]<-1
        susceptible[i]<-0
        #   print(paste0("At time=",timeSteps , " id= ",i, "is removed from susceptibles"))
      }
      exposed_incidence[i]<-exposed[i]
    }
    susceptibleMat<-cbind(susceptibleMat, susceptible)
    # Transition from E---> I (TBD on exposed at previous time step)
    for(i in vertex_id[exposedPrev == 1]){
      SuccessInfection<-runif(1) < sigmaz
      if(SuccessInfection){
        infected[i]<-1
        exposed[i]<-0
        incidence[i] <- 1
      }
    }
    # Transition from I---> R (TBD on infected at previous time step)
    for(i in vertex_id[infectedPrev == 1]){
      SuccessRecover<-runif(1) < gamma
      if(SuccessRecover){
        #infected[i]<-NA ## NAs mark who has eventually been infected
        infected[i]<-0 ## NAs mark who has eventually been infected
      }
    }
    # Append current status to matrix and update time variable
    infectedMat <- cbind(infectedMat, infected)
    incidenceMat <- cbind(incidenceMat, incidence)
    # exposedMat <- cbind(exposedMat, exposed)
    # exposed_incidenceMat <- cbind(exposed_incidenceMat, exposed_incidence)
    # hhfoiMat <- cbind(hhfoiMat, hhfoi)
    # comfoiMat <- cbind(comfoiMat, comfoi)
    deltaSusceptible<-susceptible-susceptiblePrev
    if(any(deltaSusceptible>0)){
      print(paste0("List of ids that at time ", timeSteps, " have become susceptible again"))
      print(which(deltaSusceptible<0))
      
    }
    infected_noNAs<-infected
    infected_noNAs[is.na(infected_noNAs)]<-0
    dum_check<-susceptible+infected_noNAs+exposed
    if(any(dum_check>1)){
      print(paste0("List of ids that at time ", timeSteps, " belong to multiple status"))
      print(which(dum_check>1))
      
    }
    
    timeSteps <- timeSteps + 1    
  }
  
  infectedMat_censored<-infectedMat
  for(irow in 1:nrow(infectedMat_censored)){
    censoring_time_this_part<-df_epi_info$max_t_simulation[df_epi_info$memid==vertexNames[irow]]
    infectedMat_censored[irow,censoring_time_this_part:ncol(infectedMat_censored)]<-0
  }
  ever_infected<-base::rowSums(infectedMat_censored,na.rm = TRUE)
  ever_infected[ever_infected>0]<-1
  df_final_infection<-data.frame( "memid"=network::get.vertex.attribute(HH_network,"vertex.names"),
                                  "index"=network::get.vertex.attribute(HH_network,"index"),
                                  "ever_infected"=ever_infected)
  return(df_final_infection) 
}
SEIR_simulate_two_networks <-function(HH_network_no_isolation,HH_network_w_isolation,beta1, beta2){
  sigmaz<-1./1.
  gamma<-1./(7)    
  ## but there is the MINIMUM of 1 day of infectiousness
  df_epi_info<-read.csv(file='../epi_fitting/epidemic_info/epidemic_info_net_onset_yest.csv')  
  df_epi_info$cdc_primary_onsetdt<-as.Date(df_epi_info$cdc_primary_onsetdt)  ## Define as date
  df_epi_info$delay_isolation<-df_epi_info$delay_isolation+2
  ## (transition I->R done on infected at t-1)
  N_timesteps<-max(df_epi_info$max_t_simulation)  ## NB it is 42
  ## (transition I->R done on infected at t-1)
  nrVertices <- network.size(HH_network_no_isolation)
  vertexNames <- HH_network_no_isolation %v% "vertex.names"  ## Needed to check isolation status
  infected <- rep(0, nrVertices)
  susceptible <- rep(1, nrVertices)
  ids_to_infec<-which(network::get.vertex.attribute(x = HH_network_no_isolation,"index")==1)
  infected[ids_to_infec] <- 1		
  susceptible[ids_to_infec] <- 0		
  infectedMat <- infected
  exposed <- rep(0, nrVertices)
  incidenceMat <- infected
  susceptibleMat<-susceptible
  # Individual attributes
  vertex_id <- 1:nrVertices
  timeSteps <- 1
  while(timeSteps<N_timesteps){ 
    hhfoi <- rep(0, nrVertices)
    # Expose each susceptible with a certain probability
    incidence <- rep(0, nrVertices)
    exposed_incidence <- rep(0, nrVertices)
    infectedPrev <- infected	# indicator of infection in the previous time step
    exposedPrev  <- exposed   # indicator of exposure in the previous time step
    susceptiblePrev<-susceptible
    # Transition from S---> E
    for(i in vertex_id[susceptible == 1]){
      vid<-vertexNames[i]
      t_isolation_this_vid<-df_epi_info$delay_isolation[which(df_epi_info$memid==vid)]
      #print(paste0("t=", timeSteps, "i= ",i," exposed= ", exposed[i]," susceptible= ",susceptible[i]))
      ### Choose network to consider depending on time of isolations
      if(timeSteps>t_isolation_this_vid){hh_contacts <- get.neighborhood(HH_network_w_isolation, i)}else{
        hh_contacts <- get.neighborhood(HH_network_no_isolation, i)
      }
      
      nrInfHH <- sum(infectedPrev[vertex_id %in% hh_contacts ], na.rm = TRUE)
      infProb <- 1-(1-beta1)^nrInfHH*(1-beta2)
      #comfoi[i] <- beta2
      SuccessExposure<-runif(1) < infProb
      if(SuccessExposure){
        exposed[i]<-1
        susceptible[i]<-0
        #   print(paste0("At time=",timeSteps , " id= ",i, "is removed from susceptibles"))
      }
      exposed_incidence[i]<-exposed[i]
    }
    susceptibleMat<-cbind(susceptibleMat, susceptible)
    # Transition from E---> I (TBD on exposed at previous time step)
    for(i in vertex_id[exposedPrev == 1]){
      SuccessInfection<-runif(1) < sigmaz
      if(SuccessInfection){
        infected[i]<-1
        exposed[i]<-0
        incidence[i] <- 1
      }
    }
    # Transition from I---> R (TBD on infected at previous time step)
    for(i in vertex_id[infectedPrev == 1]){
      SuccessRecover<-runif(1) < gamma
      if(SuccessRecover){
        #infected[i]<-NA ## NAs mark who has eventually been infected
        infected[i]<-0 ## NAs mark who has eventually been infected
      }
    }
    # Append current status to matrix and update time variable
    infectedMat <- cbind(infectedMat, infected)
    incidenceMat <- cbind(incidenceMat, incidence)
    timeSteps <- timeSteps + 1    
  }
  
  infectedMat_censored<-infectedMat
  for(irow in 1:nrow(infectedMat_censored)){
    censoring_time_this_part<-df_epi_info$max_t_simulation[df_epi_info$memid==vertexNames[irow]]
    infectedMat_censored[irow,censoring_time_this_part:ncol(infectedMat_censored)]<-0
  }
  
  
  ever_infected<-base::rowSums(infectedMat_censored,na.rm = TRUE)
  ever_infected[ever_infected>0]<-1
  df_final_infection<-data.frame( "memid"=network::get.vertex.attribute(HH_network_no_isolation,"vertex.names"),
                                  "index"=network::get.vertex.attribute(HH_network_no_isolation,"index"),
                                  "ever_infected"=ever_infected)
  return(df_final_infection) 
}
# Variable community transmisison (based on county data)
SEIR_simulate_two_networks_var_comm_trans <-function(HH_network_no_isolation,HH_network_w_isolation,beta1, beta2){
  ### Read the comm incidence
  df_comm_transm<-read.csv("../epi_fitting/community_transmission.csv")
  sigmaz<-1./1.
  gamma<-1./(7)    ### Average duration_pos 7
  ## but there is the MINIMUM of 1 day of infectiousness
  df_epi_info<-read.csv(file='../epi_fitting/epidemic_info/epidemic_info_net_onset_yest.csv')  
  df_epi_info$cdc_primary_onsetdt<-as.Date(df_epi_info$cdc_primary_onsetdt)  ## Define as date
  df_epi_info$delay_isolation<-df_epi_info$delay_isolation+2
  ## (transition I->R done on infected at t-1)
  
  N_timesteps<-max(df_epi_info$max_t_simulation)  ## NB it is 42
  ### This is the maximum of variable "max_t_simulation"
  ## which is defined as 28 days (follow up) + time from cdc_primary_onset (time of symptom onset for index) to doe_dt (enrollment)
  ## So 28 + doe_dt - cdc_primary_onset
  #### NB: later, only infections happening before max_t_simulation are considered
  nrVertices <- network.size(HH_network_no_isolation)
  vertexNames <- HH_network_no_isolation %v% "vertex.names"  ## Needed to check isolation status
  infected <- rep(0, nrVertices)
  susceptible <- rep(1, nrVertices)
  ids_to_infec<-which(network::get.vertex.attribute(x = HH_network_no_isolation,"index")==1)
  infected[ids_to_infec] <- 1		
  susceptible[ids_to_infec] <- 0		
  infectedMat <- infected
  exposed <- rep(0, nrVertices)
  incidenceMat <- infected
  susceptibleMat<-susceptible
  # Individual attributes
  vertex_id <- 1:nrVertices
  timeSteps <- 1
  while(timeSteps<N_timesteps){ 
    hhfoi <- rep(0, nrVertices)
    # Expose each susceptible with a certain probability
    incidence <- rep(0, nrVertices)
    exposed_incidence <- rep(0, nrVertices)
    infectedPrev <- infected	# indicator of infection in the previous time step
    exposedPrev  <- exposed   # indicator of exposure in the previous time step
    susceptiblePrev<-susceptible
    # Transition from S---> E
    for(i in vertex_id[susceptible == 1]){
      vid<-vertexNames[i]
      t_isolation_this_vid<-df_epi_info$delay_isolation[which(df_epi_info$memid==vid)]
      #print(paste0("t=", timeSteps, "i= ",i," exposed= ", exposed[i]," susceptible= ",susceptible[i]))
      ### Choose network to consider depending on time of isolations
      if(timeSteps>t_isolation_this_vid){hh_contacts <- get.neighborhood(HH_network_w_isolation, i)}else{
        hh_contacts <- get.neighborhood(HH_network_no_isolation, i)
      }
      
      nrInfHH <- sum(infectedPrev[vertex_id %in% hh_contacts ], na.rm = TRUE)
      
      
      ### Match household to timing and site
      absolute_time<-df_epi_info$cdc_primary_onsetdt[df_epi_info$memid==vid]+timeSteps
      site_of_id<-df_epi_info$site[df_epi_info$memid==vid]
      date_to_next_wednesday<-lubridate::ceiling_date(absolute_time,unit="weeks",week_start=3)
      ### debug: lubridate::wday(lubridate::ceiling_date(absolute_time,unit="weeks",week_start=3),label=TRUE)
      ### Get the relevant factor
      selection<-df_comm_transm$location==site_of_id& df_comm_transm$date==date_to_next_wednesday
      commModulationFactor<-df_comm_transm$ratio_weekly_cases[selection]
      #      print(paste0("For participant ",vid, " after ",timeSteps, " timesteps the date is ", absolute_time, " which is ceiled to ", date_to_next_wednesday, " for ", commModulationFactor ))
      infProb <- 1-(1-beta1)^nrInfHH*(1-beta2*commModulationFactor)
      #comfoi[i] <- beta2
      SuccessExposure<-runif(1) < infProb
      if(SuccessExposure){
        exposed[i]<-1
        susceptible[i]<-0
        #   print(paste0("At time=",timeSteps , " id= ",i, "is removed from susceptibles"))
      }
      exposed_incidence[i]<-exposed[i]
    }
    susceptibleMat<-cbind(susceptibleMat, susceptible)
    # Transition from E---> I (TBD on exposed at previous time step)
    for(i in vertex_id[exposedPrev == 1]){
      SuccessInfection<-runif(1) < sigmaz
      if(SuccessInfection){
        infected[i]<-1
        exposed[i]<-0
        incidence[i] <- 1
      }
    }
    # Transition from I---> R (TBD on infected at previous time step)
    for(i in vertex_id[infectedPrev == 1]){
      SuccessRecover<-runif(1) < gamma
      if(SuccessRecover){
        #infected[i]<-NA ## NAs mark who has eventually been infected
        infected[i]<-0 ## NAs mark who has eventually been infected
      }
    }
    # Append current status to matrix and update time variable
    infectedMat <- cbind(infectedMat, infected)
    incidenceMat <- cbind(incidenceMat, incidence)
    timeSteps <- timeSteps + 1    
  }
  
  ## Ignoring infections happening after the max_time in the study
  infectedMat_censored<-infectedMat
  for(irow in 1:nrow(infectedMat_censored)){
    censoring_time_this_part<-df_epi_info$max_t_simulation[df_epi_info$memid==vertexNames[irow]]
    infectedMat_censored[irow,censoring_time_this_part:ncol(infectedMat_censored)]<-0
  }
  ever_infected<-base::rowSums(infectedMat_censored,na.rm = TRUE)
  ever_infected[ever_infected>0]<-1
  df_final_infection<-data.frame( "memid"=network::get.vertex.attribute(HH_network_no_isolation,"vertex.names"),
                                  "index"=network::get.vertex.attribute(HH_network_no_isolation,"index"),
                                  "ever_infected"=ever_infected)
  return(df_final_infection) 
}
# Variable community tranmission (based on state data)
SEIR_simulate_two_networks_var_comm_trans_state <-function(HH_network_no_isolation,HH_network_w_isolation,beta1, beta2){
  ### Read the comm incidence
  df_comm_transm<-read.csv("../epi_fitting/community_transmission_state.csv")
  sigmaz<-1./1.
  gamma<-1./(7)    ### Average duration_pos 7
  ## but there is the MINIMUM of 1 day of infectiousness
  df_epi_info<-read.csv(file='../epi_fitting/epidemic_info/epidemic_info_net_onset_yest.csv')  
  df_epi_info$cdc_primary_onsetdt<-as.Date(df_epi_info$cdc_primary_onsetdt)  ## Define as date
  df_epi_info$delay_isolation<-df_epi_info$delay_isolation+2
  ## (transition I->R done on infected at t-1)
  
  N_timesteps<-max(df_epi_info$max_t_simulation)  ## NB it is 42
  ### This is the maximum of variable "max_t_simulation"
  ## which is defined as 28 days (follow up) + time from cdc_primary_onset (time of symptom onset for index) to doe_dt (enrollment)
  ## So 28 + doe_dt - cdc_primary_onset
  #### NB: later, only infections happening before max_t_simulation are considered
  nrVertices <- network.size(HH_network_no_isolation)
  vertexNames <- HH_network_no_isolation %v% "vertex.names"  ## Needed to check isolation status
  infected <- rep(0, nrVertices)
  susceptible <- rep(1, nrVertices)
  ids_to_infec<-which(network::get.vertex.attribute(x = HH_network_no_isolation,"index")==1)
  infected[ids_to_infec] <- 1		
  susceptible[ids_to_infec] <- 0		
  infectedMat <- infected
  exposed <- rep(0, nrVertices)
  incidenceMat <- infected
  susceptibleMat<-susceptible
  # Individual attributes
  vertex_id <- 1:nrVertices
  timeSteps <- 1
  while(timeSteps<N_timesteps){ 
    hhfoi <- rep(0, nrVertices)
    # Expose each susceptible with a certain probability
    incidence <- rep(0, nrVertices)
    exposed_incidence <- rep(0, nrVertices)
    infectedPrev <- infected	# indicator of infection in the previous time step
    exposedPrev  <- exposed   # indicator of exposure in the previous time step
    susceptiblePrev<-susceptible
    # Transition from S---> E
    for(i in vertex_id[susceptible == 1]){
      vid<-vertexNames[i]
      t_isolation_this_vid<-df_epi_info$delay_isolation[which(df_epi_info$memid==vid)]
      #print(paste0("t=", timeSteps, "i= ",i," exposed= ", exposed[i]," susceptible= ",susceptible[i]))
      ### Choose network to consider depending on time of isolations
      if(timeSteps>t_isolation_this_vid){hh_contacts <- get.neighborhood(HH_network_w_isolation, i)}else{
        hh_contacts <- get.neighborhood(HH_network_no_isolation, i)
      }
      
      nrInfHH <- sum(infectedPrev[vertex_id %in% hh_contacts ], na.rm = TRUE)
      
      
      ### Match household to timing and site
      absolute_time<-df_epi_info$cdc_primary_onsetdt[df_epi_info$memid==vid]+timeSteps
      site_of_id<-df_epi_info$site[df_epi_info$memid==vid]
      date_to_next_wednesday<-lubridate::ceiling_date(absolute_time,unit="weeks",week_start=3)
      ### debug: lubridate::wday(lubridate::ceiling_date(absolute_time,unit="weeks",week_start=3),label=TRUE)
      ### Get the relevant factor
      selection<-df_comm_transm$location==site_of_id& df_comm_transm$date==date_to_next_wednesday
      commModulationFactor<-df_comm_transm$ratio_weekly_cases[selection]
      #      print(paste0("For participant ",vid, " after ",timeSteps, " timesteps the date is ", absolute_time, " which is ceiled to ", date_to_next_wednesday, " for ", commModulationFactor ))
      infProb <- 1-(1-beta1)^nrInfHH*(1-beta2*commModulationFactor)
      #comfoi[i] <- beta2
      SuccessExposure<-runif(1) < infProb
      if(SuccessExposure){
        exposed[i]<-1
        susceptible[i]<-0
        #   print(paste0("At time=",timeSteps , " id= ",i, "is removed from susceptibles"))
      }
      exposed_incidence[i]<-exposed[i]
    }
    susceptibleMat<-cbind(susceptibleMat, susceptible)
    # Transition from E---> I (TBD on exposed at previous time step)
    for(i in vertex_id[exposedPrev == 1]){
      SuccessInfection<-runif(1) < sigmaz
      if(SuccessInfection){
        infected[i]<-1
        exposed[i]<-0
        incidence[i] <- 1
      }
    }
    # Transition from I---> R (TBD on infected at previous time step)
    for(i in vertex_id[infectedPrev == 1]){
      SuccessRecover<-runif(1) < gamma
      if(SuccessRecover){
        #infected[i]<-NA ## NAs mark who has eventually been infected
        infected[i]<-0 ## NAs mark who has eventually been infected
      }
    }
    # Append current status to matrix and update time variable
    infectedMat <- cbind(infectedMat, infected)
    incidenceMat <- cbind(incidenceMat, incidence)
    timeSteps <- timeSteps + 1    
  }
  
  ## Ignoring infections happening after the max_time in the study
  infectedMat_censored<-infectedMat
  for(irow in 1:nrow(infectedMat_censored)){
    censoring_time_this_part<-df_epi_info$max_t_simulation[df_epi_info$memid==vertexNames[irow]]
    infectedMat_censored[irow,censoring_time_this_part:ncol(infectedMat_censored)]<-0
  }
  ever_infected<-base::rowSums(infectedMat_censored,na.rm = TRUE)
  ever_infected[ever_infected>0]<-1
  df_final_infection<-data.frame( "memid"=network::get.vertex.attribute(HH_network_no_isolation,"vertex.names"),
                                  "index"=network::get.vertex.attribute(HH_network_no_isolation,"index"),
                                  "ever_infected"=ever_infected)
  return(df_final_infection) 
}


# Variable infectivity
# Variable community transmisison
SEIR_simulate_two_networks_var_infectivity <-function(HH_network_no_isolation,HH_network_w_isolation,beta1, beta2){
  sigmaz<-1./1
  gammaz<-1./7    ### Average duration of infectivity is 7
  ### Read the comm incidence
  df_comm_transm<-read.csv("../epi_fitting/community_transmission.csv")
  ## but there is the MINIMUM of 1 day of infectiousness
  df_epi_info<-read.csv(file='../epi_fitting/epidemic_info/epidemic_info_net_onset_yest.csv')  
  df_epi_info$delay_isolation<-df_epi_info$delay_isolation + 2 ## Infectious from 2 days before
  
  df_infectiousness<-read.csv("../epi_fitting/infectiousness_over_time.csv")
  dum_df<-data.frame(time=(12:100), infectuosness=rep(0,89))
  df_infectiousness<-rbind(df_infectiousness,dum_df)  ## File defined up to 12
  ## (transition I->R done on infected at t-1)
  N_timesteps<-max(df_epi_info$max_t_simulation)  ## NB it is 42
  ### This is the maximum of variable "max_t_simulation"
  ## which is defined as 28 days (follow up) + time from cdc_primary_onset (time of symptom onset for index) to doe_dt (enrollment)
  ## So 28 + doe_dt - cdc_primary_onset
  #### NB: later, only infections happening before max_t_simulation are considered
  nrVertices <- network.size(HH_network_no_isolation)
  vertexNames <- HH_network_no_isolation %v% "vertex.names"  ## Needed to check isolation status
  infected <- rep(0, nrVertices)
  time_of_infection<-rep(NA, nrVertices)
  time_of_recovery<-rep(NA, nrVertices)
  susceptible <- rep(1, nrVertices)
  ids_to_infec<-which(network::get.vertex.attribute(x = HH_network_no_isolation,"index")==1)
  infected[ids_to_infec] <- 1
  time_of_infection[ids_to_infec]<-0
  susceptible[ids_to_infec] <- 0		
  infectedMat <- infected
  exposed <- rep(0, nrVertices)
  incidenceMat <- infected
  susceptibleMat<-susceptible
  # Individual attributes
  vertex_id <- 1:nrVertices
  timeSteps <- 1
  while(timeSteps<N_timesteps){ 
    hhfoi <- rep(0, nrVertices)
    # Expose each susceptible with a certain probability
    incidence <- rep(0, nrVertices)
    exposed_incidence <- rep(0, nrVertices)
    infectedPrev <- infected	# indicator of infection in the previous time step
    exposedPrev  <- exposed   # indicator of exposure in the previous time step
    susceptiblePrev<-susceptible
    # Transition from S---> E
    for(i in vertex_id[susceptible == 1]){
      vid<-vertexNames[i]
      t_isolation_this_vid<-df_epi_info$delay_isolation[which(df_epi_info$memid==vid)]
      #print(paste0("t=", timeSteps, "i= ",i," exposed= ", exposed[i]," susceptible= ",susceptible[i]))
      ### Choose network to consider depending on time of isolations
      if(timeSteps>t_isolation_this_vid){hh_contacts <- get.neighborhood(HH_network_w_isolation, i)}else{
        hh_contacts <- get.neighborhood(HH_network_no_isolation, i)
      }
      #nrInfHH <- sum(infectedPrev[vertex_id %in% hh_contacts ], na.rm = TRUE)
      ## Loop on hh_contacts and sum infection parameters
      prod_escape_infection_hh<-1
      if(length(hh_contacts)>0){
        for(i_hh_member in 1:length(hh_contacts)){
          this_contact_id<-hh_contacts[i_hh_member]
          if(infectedPrev[this_contact_id]==1){
            time_from_infection= timeSteps-time_of_infection[this_contact_id]
            
            factor<-df_infectiousness$infectuosness[time_from_infection]
            contribution<-(1-beta1*factor)
            prod_escape_infection_hh<-prod_escape_infection_hh*contribution
            #print(paste0("contact with ",this_contact_id,"t_from_infection= ",time_from_infection," infection prob= ",infProb))
          }
        }
      }
      ## For each element 1-beta1*factor
      infProb <- 1-prod_escape_infection_hh*(1-beta2)
      #print(paste0("Escape infection= ",prod_escape_infection_hh," infection prob= ",infProb))
      #comfoi[i] <- beta2
      SuccessExposure<-runif(1) < infProb
      if(SuccessExposure){
        exposed[i]<-1
        susceptible[i]<-0
        #   print(paste0("At time=",timeSteps , " id= ",i, "is removed from susceptibles"))
      }
      exposed_incidence[i]<-exposed[i]
    }
    susceptibleMat<-cbind(susceptibleMat, susceptible)
    # Transition from E---> I (TBD on exposed at previous time step)
    for(i in vertex_id[exposedPrev == 1]){
      SuccessInfection<-runif(1) < sigmaz
      if(SuccessInfection){
        infected[i]<-1
        exposed[i]<-0
        incidence[i] <- 1
        time_of_infection[i]<-timeSteps
        
      }
    }
    # Transition from I---> R (TBD on infected at previous time step)
    for(i in vertex_id[infectedPrev == 1]){
      SuccessRecover<-runif(1) < gammaz
      if(SuccessRecover){
        #infected[i]<-NA ## NAs mark who has eventually been infected
        infected[i]<-0 
        time_of_recovery[i]<-timeSteps
      }
    }
    # Append current status to matrix and update time variable
    infectedMat <- cbind(infectedMat, infected)
    incidenceMat <- cbind(incidenceMat, incidence)
    timeSteps <- timeSteps + 1    
  }
  
  ## Ignoring infections happening after the max_time in the study
  infectedMat_censored<-infectedMat
  for(irow in 1:nrow(infectedMat_censored)){
    censoring_time_this_part<-df_epi_info$max_t_simulation[df_epi_info$memid==vertexNames[irow]]
    infectedMat_censored[irow,censoring_time_this_part:ncol(infectedMat_censored)]<-0
  }
  ever_infected<-base::rowSums(infectedMat_censored,na.rm = TRUE)
  ever_infected[ever_infected>0]<-1
  df_final_infection<-data.frame( "memid"=network::get.vertex.attribute(HH_network_no_isolation,"vertex.names"),
                                  "index"=network::get.vertex.attribute(HH_network_no_isolation,"index"),
                                  "ever_infected"=ever_infected)
  return(df_final_infection) 
}


#### Simulation for GOF
SEIR_simulate_two_networks_return_incidence<-function(beta1,beta2){
  sigmaz<-1./3
  gamma<-1./(9.796-1)    ### Average duration_pos is 9.796 (computed from swabs)
  ## but there is the MINIMUM of 1 day of infectiousness
  df_epi_info<-read.csv(file='../epi_fitting/epidemic_info/epidemic_info_net_onset_yest.csv')  
  df_epi_info$cdc_primary_onsetdt<-as.Date(df_epi_info$cdc_primary_onsetdt)  ## Define as date
  df_epi_info$delay_isolation<-df_epi_info$delay_isolation
  ## (transition I->R done on infected at t-1)
  N_timesteps<-max(df_epi_info$max_t_simulation)  ## NB it is 42
  ## (transition I->R done on infected at t-1)
  nrVertices <- network.size(HH_network_no_isolation)
  vertexNames <- HH_network_no_isolation %v% "vertex.names"  ## Needed to check isolation status
  infected <- rep(0, nrVertices)
  susceptible <- rep(1, nrVertices)
  ids_to_infec<-which(network::get.vertex.attribute(x = HH_network_no_isolation,"index")==1)
  infected[ids_to_infec] <- 1		
  susceptible[ids_to_infec] <- 0		
  infectedMat <- infected
  exposed <- rep(0, nrVertices)
  incidenceMat <- infected
  susceptibleMat<-susceptible
  # Individual attributes
  vertex_id <- 1:nrVertices
  timeSteps <- 1
  while(timeSteps<N_timesteps){ 
    hhfoi <- rep(0, nrVertices)
    # Expose each susceptible with a certain probability
    incidence <- rep(0, nrVertices)
    exposed_incidence <- rep(0, nrVertices)
    infectedPrev <- infected	# indicator of infection in the previous time step
    exposedPrev  <- exposed   # indicator of exposure in the previous time step
    susceptiblePrev<-susceptible
    # Transition from S---> E
    for(i in vertex_id[susceptible == 1]){
      vid<-vertexNames[i]
      t_isolation_this_vid<-df_epi_info$delay_isolation[which(df_epi_info$memid==vid)]
      #print(paste0("t=", timeSteps, "i= ",i," exposed= ", exposed[i]," susceptible= ",susceptible[i]))
      ### Choose network to consider depending on time of isolations
      if(timeSteps>t_isolation_this_vid){hh_contacts <- get.neighborhood(HH_network_w_isolation, i)}else{
        hh_contacts <- get.neighborhood(HH_network_no_isolation, i)
      }
      
      nrInfHH <- sum(infectedPrev[vertex_id %in% hh_contacts ], na.rm = TRUE)
      infProb <- 1-(1-beta1)^nrInfHH*(1-beta2)
      #comfoi[i] <- beta2
      SuccessExposure<-runif(1) < infProb
      if(SuccessExposure){
        exposed[i]<-1
        susceptible[i]<-0
        #   print(paste0("At time=",timeSteps , " id= ",i, "is removed from susceptibles"))
      }
      exposed_incidence[i]<-exposed[i]
    }
    susceptibleMat<-cbind(susceptibleMat, susceptible)
    # Transition from E---> I (TBD on exposed at previous time step)
    for(i in vertex_id[exposedPrev == 1]){
      SuccessInfection<-runif(1) < sigmaz
      if(SuccessInfection){
        infected[i]<-1
        exposed[i]<-0
        incidence[i] <- 1
      }
    }
    # Transition from I---> R (TBD on infected at previous time step)
    for(i in vertex_id[infectedPrev == 1]){
      SuccessRecover<-runif(1) < gamma
      if(SuccessRecover){
        #infected[i]<-NA ## NAs mark who has eventually been infected
        infected[i]<-0 ## NAs mark who has eventually been infected
      }
    }
    # Append current status to matrix and update time variable
    infectedMat <- cbind(infectedMat, infected)
    incidenceMat <- cbind(incidenceMat, incidence)
    timeSteps <- timeSteps + 1    
  }
  infectedMat_censored<-infectedMat
  for(irow in 1:nrow(infectedMat_censored)){
    censoring_time_this_part<-df_epi_info$max_t_simulation[df_epi_info$memid==vertexNames[irow]]
    infectedMat_censored[irow,censoring_time_this_part:ncol(infectedMat_censored)]<-0
  }
  
  ever_infected<-base::rowSums(infectedMat_censored,na.rm = TRUE)
  incidenceVector<-base::colSums(incidenceMat,na.rm = TRUE)
  ## Removing t=0 (incidence is defined as N indexes)
  incidenceVector<-incidenceVector[2:length(incidenceVector)]
  
  
  return(incidenceVector)
}
SEIR_simulate_two_networks_var_comm_trans_return_incidence<-function(beta1, beta2){
  ### Read the comm incidence
  df_comm_transm<-read.csv("../epi_fitting/community_transmission.csv")
  sigmaz<-1./3
  gamma<-1./(9.796-1)    ### Average duration_pos is 9.796 (computed from swabs)
  ## but there is the MINIMUM of 1 day of infectiousness
  df_epi_info<-read.csv(file='../epi_fitting/epidemic_info/epidemic_info_net_onset_yest.csv')  
  df_epi_info$cdc_primary_onsetdt<-as.Date(df_epi_info$cdc_primary_onsetdt)  ## Define as date
  df_epi_info$delay_isolation<-df_epi_info$delay_isolation#+ceiling(1./sigmaz) ## add sigma to time
  ## (transition I->R done on infected at t-1)
  
  N_timesteps<-max(df_epi_info$max_t_simulation)  ## NB it is 42
  ### This is the maximum of variable "max_t_simulation"
  ## which is defined as 28 days (follow up) + time from cdc_primary_onset (time of symptom onset for index) to doe_dt (enrollment)
  ## So 28 + doe_dt - cdc_primary_onset
  #### NB: later, only infections happening before max_t_simulation are considered
  nrVertices <- network.size(HH_network_no_isolation)
  vertexNames <- HH_network_no_isolation %v% "vertex.names"  ## Needed to check isolation status
  infected <- rep(0, nrVertices)
  susceptible <- rep(1, nrVertices)
  ids_to_infec<-which(network::get.vertex.attribute(x = HH_network_no_isolation,"index")==1)
  infected[ids_to_infec] <- 1		
  susceptible[ids_to_infec] <- 0		
  infectedMat <- infected
  exposed <- rep(0, nrVertices)
  incidenceMat <- infected
  susceptibleMat<-susceptible
  # Individual attributes
  vertex_id <- 1:nrVertices
  timeSteps <- 1
  while(timeSteps<N_timesteps){ 
    hhfoi <- rep(0, nrVertices)
    # Expose each susceptible with a certain probability
    incidence <- rep(0, nrVertices)
    exposed_incidence <- rep(0, nrVertices)
    infectedPrev <- infected	# indicator of infection in the previous time step
    exposedPrev  <- exposed   # indicator of exposure in the previous time step
    susceptiblePrev<-susceptible
    # Transition from S---> E
    for(i in vertex_id[susceptible == 1]){
      vid<-vertexNames[i]
      t_isolation_this_vid<-df_epi_info$delay_isolation[which(df_epi_info$memid==vid)]
      #print(paste0("t=", timeSteps, "i= ",i," exposed= ", exposed[i]," susceptible= ",susceptible[i]))
      ### Choose network to consider depending on time of isolations
      if(timeSteps>t_isolation_this_vid){hh_contacts <- get.neighborhood(HH_network_w_isolation, i)}else{
        hh_contacts <- get.neighborhood(HH_network_no_isolation, i)
      }
      
      nrInfHH <- sum(infectedPrev[vertex_id %in% hh_contacts ], na.rm = TRUE)
      
      
      ### Match household to timing and site
      absolute_time<-df_epi_info$cdc_primary_onsetdt[df_epi_info$memid==vid]+timeSteps
      site_of_id<-df_epi_info$site[df_epi_info$memid==vid]
      date_to_next_wednesday<-lubridate::ceiling_date(absolute_time,unit="weeks",week_start=3)
      ### debug: lubridate::wday(lubridate::ceiling_date(absolute_time,unit="weeks",week_start=3),label=TRUE)
      ### Get the relevant factor
      selection<-df_comm_transm$location==site_of_id& df_comm_transm$date==date_to_next_wednesday
      commModulationFactor<-df_comm_transm$ratio_weekly_cases[selection]
      #      print(paste0("For participant ",vid, " after ",timeSteps, " timesteps the date is ", absolute_time, " which is ceiled to ", date_to_next_wednesday, " for ", commModulationFactor ))
      infProb <- 1-(1-beta1)^nrInfHH*(1-beta2*commModulationFactor)
      #comfoi[i] <- beta2
      SuccessExposure<-runif(1) < infProb
      if(SuccessExposure){
        exposed[i]<-1
        susceptible[i]<-0
        #   print(paste0("At time=",timeSteps , " id= ",i, "is removed from susceptibles"))
      }
      exposed_incidence[i]<-exposed[i]
    }
    susceptibleMat<-cbind(susceptibleMat, susceptible)
    # Transition from E---> I (TBD on exposed at previous time step)
    for(i in vertex_id[exposedPrev == 1]){
      SuccessInfection<-runif(1) < sigmaz
      if(SuccessInfection){
        infected[i]<-1
        exposed[i]<-0
        incidence[i] <- 1
      }
    }
    # Transition from I---> R (TBD on infected at previous time step)
    for(i in vertex_id[infectedPrev == 1]){
      SuccessRecover<-runif(1) < gamma
      if(SuccessRecover){
        #infected[i]<-NA ## NAs mark who has eventually been infected
        infected[i]<-0 ## NAs mark who has eventually been infected
      }
    }
    # Append current status to matrix and update time variable
    infectedMat <- cbind(infectedMat, infected)
    incidenceMat <- cbind(incidenceMat, incidence)
    timeSteps <- timeSteps + 1    
  }
  
  ## Ignoring infections happening after the max_time in the study
  infectedMat_censored<-infectedMat
  for(irow in 1:nrow(infectedMat_censored)){
    censoring_time_this_part<-df_epi_info$max_t_simulation[df_epi_info$memid==vertexNames[irow]]
    infectedMat_censored[irow,censoring_time_this_part:ncol(infectedMat_censored)]<-0
  }
  
  ever_infected<-base::rowSums(infectedMat_censored,na.rm = TRUE)
  incidenceVector<-base::colSums(incidenceMat,na.rm = TRUE)
  ## Removing t=0 (incidence is defined as N indexes)
  incidenceVector<-incidenceVector[2:length(incidenceVector)]
  
  
  return(incidenceVector)
}
#### Wrappers for fit


SEIR_simulate_two_networks_forBOptim_Binloglik_Ncases_by_HH<-function(beta1,beta2){
  sigmaz<-1./1.  ## 3 days to develop symptoms, but infectious from 2 days before symtpoms
  gammaz<-1./7    ### Average duration of infectivity is 7
  ## but there is the MINIMUM of 1 day of infectiousness
  df_epi_info<-read.csv(file='../epi_fitting/epidemic_info/epidemic_info_net_onset_yest.csv')  
  df_epi_info$delay_isolation<-df_epi_info$delay_isolation + 2 ## Infectious from 2 days before
  N_timesteps<-max(df_epi_info$max_t_simulation)  ## NB it is 42
  ## (transition I->R done on infected at t-1)
  nrVertices <- network.size(HH_network_no_isolation)
  vertexNames <- HH_network_no_isolation %v% "vertex.names"  ## Needed to check isolation status
  infected <- rep(0, nrVertices)
  susceptible <- rep(1, nrVertices)
  ids_to_infec<-which(network::get.vertex.attribute(x = HH_network_no_isolation,"index")==1)
  infected[ids_to_infec] <- 1		
  susceptible[ids_to_infec] <- 0		
  infectedMat <- infected
  exposed <- rep(0, nrVertices)
  incidenceMat <- infected
  susceptibleMat<-susceptible
  # Individual attributes
  vertex_id <- 1:nrVertices
  timeSteps <- 1
  while(timeSteps<N_timesteps){ 
    hhfoi <- rep(0, nrVertices)
    # Expose each susceptible with a certain probability
    incidence <- rep(0, nrVertices)
    exposed_incidence <- rep(0, nrVertices)
    infectedPrev <- infected	# indicator of infection in the previous time step
    exposedPrev  <- exposed   # indicator of exposure in the previous time step
    susceptiblePrev<-susceptible
    # Transition from S---> E
    for(i in vertex_id[susceptible == 1]){
      vid<-vertexNames[i]
      t_isolation_this_vid<-df_epi_info$delay_isolation[which(df_epi_info$memid==vid)]
      #print(paste0("t=", timeSteps, "i= ",i," exposed= ", exposed[i]," susceptible= ",susceptible[i]))
      ### Choose network to consider depending on time of isolations
      if(timeSteps>t_isolation_this_vid){hh_contacts <- get.neighborhood(HH_network_w_isolation, i)}else{
        hh_contacts <- get.neighborhood(HH_network_no_isolation, i)
      }
      nrInfHH <- sum(infectedPrev[vertex_id %in% hh_contacts ], na.rm = TRUE)
      infProb <- 1-(1-beta1)^nrInfHH*(1-beta2)
      #comfoi[i] <- beta2
      SuccessExposure<-runif(1) < infProb
      if(SuccessExposure){
        exposed[i]<-1
        susceptible[i]<-0
        #   print(paste0("At time=",timeSteps , " id= ",i, "is removed from susceptibles"))
      }
      exposed_incidence[i]<-exposed[i]
    }
    susceptibleMat<-cbind(susceptibleMat, susceptible)
    # Transition from E---> I (TBD on exposed at previous time step)
    for(i in vertex_id[exposedPrev == 1]){
      SuccessInfection<-runif(1) < sigmaz
      if(SuccessInfection){
        infected[i]<-1
        exposed[i]<-0
        incidence[i] <- 1
      }
    }
    # Transition from I---> R (TBD on infected at previous time step)
    for(i in vertex_id[infectedPrev == 1]){
      SuccessRecover<-runif(1) < gammaz
      if(SuccessRecover){
        #infected[i]<-NA ## NAs mark who has eventually been infected
        infected[i]<-0 ## NAs mark who has eventually been infected
      }
    }
    # Append current status to matrix and update time variable
    infectedMat <- cbind(infectedMat, infected)
    incidenceMat <- cbind(incidenceMat, incidence)
    timeSteps <- timeSteps + 1    
  }
  
  infectedMat_censored<-infectedMat
  for(irow in 1:nrow(infectedMat_censored)){
    censoring_time_this_part<-df_epi_info$max_t_simulation[df_epi_info$memid==vertexNames[irow]]
    infectedMat_censored[irow,censoring_time_this_part:ncol(infectedMat_censored)]<-0
  }
  
  
  ## Compute whether a participant has been infected
  ever_infected<-base::rowSums(infectedMat_censored,na.rm = TRUE)
  ever_infected[ever_infected>0]<-1
  ## Dropping non-index
  ever_infected<-ever_infected[which(network::get.vertex.attribute(x = HH_network_no_isolation,"index")==0)]
  df_ever_infected<-data.frame(memid=vertexNames[which(network::get.vertex.attribute(x = HH_network_no_isolation,"index")==0)],ever_infected=ever_infected)
  
  
  ### Load information about hh size of the individual
  subject_info<-read.csv("../epi_simulation/subject_data_for_analysis_of_infection.csv")
  ## merge i lmemid, seleziona su household size and generate final df
  df_ever_infected<-merge(df_ever_infected,subject_info[,c("memid","hhmem_num")])
  df_ever_infected_by_hh<-df_ever_infected %>% group_by(hhmem_num) %>% summarize(tot_infected=sum(ever_infected))
  
  
  ### Load infection data
  df_data_to_fit<-read.csv("../epi_fitting/time_infection_matrix_net_onset_yest.csv",header=TRUE)
  colnames(df_data_to_fit)[1]<-"memid"
  dum<-as.matrix(df_data_to_fit[,2:dim(df_data_to_fit)[2]]) ## Dropping the Ids name column
  df_data_to_fit$ever_infected<-base::rowSums(dum,na.rm = TRUE) ## Changing to vector of ever infected
  df_data_to_fit$ever_infected[df_data_to_fit$ever_infected>0]<-1
  ## Dropping non-index
  df_data_to_fit<-df_data_to_fit %>% filter(memid %in% vertexNames[which(network::get.vertex.attribute(x = HH_network_no_isolation,"index")==0)])
  df_data_to_fit<-df_data_to_fit[,c("memid","ever_infected")]
  df_data_to_fit<-merge(df_data_to_fit,subject_info[,c("memid","hhmem_num")])
  ## Generate final data
  data_to_fit_by_hh<-df_data_to_fit %>% group_by(hhmem_num) %>% summarize(tot_infected=sum(ever_infected))
  
  
  ## Compute likelihood as Binomial log-likelihood for each hh_size
  error<- sum(log_lik(df_ever_infected_by_hh$tot_infected,data_to_fit_by_hh$tot_infected))
  return(list("Score"=error)) 
}

SEIR_simulate_two_networks_forBOptim_Binloglik_Ncases_by_HH_no2<-function(beta1,beta2){
  sigmaz<-1./1.  ## 3 days to develop symptoms, but infectious from 2 days before symtpoms
  gammaz<-1./7    ### Average duration of infectivity is 7
  ## but there is the MINIMUM of 1 day of infectiousness
  df_epi_info<-read.csv(file='../epi_fitting/epidemic_info/epidemic_info_net_onset_yest.csv')  
  df_epi_info$delay_isolation<-df_epi_info$delay_isolation + 2 ## Infectious from 2 days before
  N_timesteps<-max(df_epi_info$max_t_simulation)  ## NB it is 42
  ## (transition I->R done on infected at t-1)
  nrVertices <- network.size(HH_network_no_isolation)
  vertexNames <- HH_network_no_isolation %v% "vertex.names"  ## Needed to check isolation status
  infected <- rep(0, nrVertices)
  susceptible <- rep(1, nrVertices)
  ids_to_infec<-which(network::get.vertex.attribute(x = HH_network_no_isolation,"index")==1)
  infected[ids_to_infec] <- 1		
  susceptible[ids_to_infec] <- 0		
  infectedMat <- infected
  exposed <- rep(0, nrVertices)
  incidenceMat <- infected
  susceptibleMat<-susceptible
  # Individual attributes
  vertex_id <- 1:nrVertices
  timeSteps <- 1
  while(timeSteps<N_timesteps){ 
    hhfoi <- rep(0, nrVertices)
    # Expose each susceptible with a certain probability
    incidence <- rep(0, nrVertices)
    exposed_incidence <- rep(0, nrVertices)
    infectedPrev <- infected	# indicator of infection in the previous time step
    exposedPrev  <- exposed   # indicator of exposure in the previous time step
    susceptiblePrev<-susceptible
    # Transition from S---> E
    for(i in vertex_id[susceptible == 1]){
      vid<-vertexNames[i]
      t_isolation_this_vid<-df_epi_info$delay_isolation[which(df_epi_info$memid==vid)]
      #print(paste0("t=", timeSteps, "i= ",i," exposed= ", exposed[i]," susceptible= ",susceptible[i]))
      ### Choose network to consider depending on time of isolations
      if(timeSteps>t_isolation_this_vid){hh_contacts <- get.neighborhood(HH_network_w_isolation, i)}else{
        hh_contacts <- get.neighborhood(HH_network_no_isolation, i)
      }
      nrInfHH <- sum(infectedPrev[vertex_id %in% hh_contacts ], na.rm = TRUE)
      infProb <- 1-(1-beta1)^nrInfHH*(1-beta2)
      #comfoi[i] <- beta2
      SuccessExposure<-runif(1) < infProb
      if(SuccessExposure){
        exposed[i]<-1
        susceptible[i]<-0
        #   print(paste0("At time=",timeSteps , " id= ",i, "is removed from susceptibles"))
      }
      exposed_incidence[i]<-exposed[i]
    }
    susceptibleMat<-cbind(susceptibleMat, susceptible)
    # Transition from E---> I (TBD on exposed at previous time step)
    for(i in vertex_id[exposedPrev == 1]){
      SuccessInfection<-runif(1) < sigmaz
      if(SuccessInfection){
        infected[i]<-1
        exposed[i]<-0
        incidence[i] <- 1
      }
    }
    # Transition from I---> R (TBD on infected at previous time step)
    for(i in vertex_id[infectedPrev == 1]){
      SuccessRecover<-runif(1) < gammaz
      if(SuccessRecover){
        #infected[i]<-NA ## NAs mark who has eventually been infected
        infected[i]<-0 ## NAs mark who has eventually been infected
      }
    }
    # Append current status to matrix and update time variable
    infectedMat <- cbind(infectedMat, infected)
    incidenceMat <- cbind(incidenceMat, incidence)
    timeSteps <- timeSteps + 1    
  }
  
  infectedMat_censored<-infectedMat
  for(irow in 1:nrow(infectedMat_censored)){
    censoring_time_this_part<-df_epi_info$max_t_simulation[df_epi_info$memid==vertexNames[irow]]
    infectedMat_censored[irow,censoring_time_this_part:ncol(infectedMat_censored)]<-0
  }
  
  
  ## Compute whether a participant has been infected
  ever_infected<-base::rowSums(infectedMat_censored,na.rm = TRUE)
  ever_infected[ever_infected>0]<-1
  ## Dropping non-index
  ever_infected<-ever_infected[which(network::get.vertex.attribute(x = HH_network_no_isolation,"index")==0)]
  df_ever_infected<-data.frame(memid=vertexNames[which(network::get.vertex.attribute(x = HH_network_no_isolation,"index")==0)],ever_infected=ever_infected)
  
  
  ### Load information about hh size of the individual
  subject_info<-read.csv("../epi_simulation/subject_data_for_analysis_of_infection.csv")
  ## merge i lmemid, seleziona su household size and generate final df
  df_ever_infected<-merge(df_ever_infected,subject_info[,c("memid","hhmem_num")])
  df_ever_infected_by_hh<-df_ever_infected %>% group_by(hhmem_num) %>% summarize(tot_infected=sum(ever_infected))
  
  
  ### Load infection data
  df_data_to_fit<-read.csv("../epi_fitting/time_infection_matrix_net_onset_yest.csv",header=TRUE)
  colnames(df_data_to_fit)[1]<-"memid"
  dum<-as.matrix(df_data_to_fit[,2:dim(df_data_to_fit)[2]]) ## Dropping the Ids name column
  df_data_to_fit$ever_infected<-base::rowSums(dum,na.rm = TRUE) ## Changing to vector of ever infected
  df_data_to_fit$ever_infected[df_data_to_fit$ever_infected>0]<-1
  ## Dropping non-index
  df_data_to_fit<-df_data_to_fit %>% filter(memid %in% vertexNames[which(network::get.vertex.attribute(x = HH_network_no_isolation,"index")==0)])
  df_data_to_fit<-df_data_to_fit[,c("memid","ever_infected")]
  df_data_to_fit<-merge(df_data_to_fit,subject_info[,c("memid","hhmem_num")])
  ## Generate final data
  data_to_fit_by_hh<-df_data_to_fit %>% group_by(hhmem_num) %>% summarize(tot_infected=sum(ever_infected))
  
  
  ## Compute likelihood as Binomial log-likelihood for each hh_size EXCLUDING 2
  error<- sum(log_lik(df_ever_infected_by_hh$tot_infected,data_to_fit_by_hh$tot_infected)[-1])
  return(list("Score"=error)) 
}

SEIR_simulate_two_networks_forSimpleOptim_Binloglik_Ncases_by_HH<-function(full_list_params){
  beta1<-full_list_params[1] 
  beta2<-full_list_params[2] 
  sigmaz<-1./1.  ## 3 days to develop symptoms, but infectious from 2 days before symtpoms
  gammaz<-1./7    ### Average duration of infectivity is 7
  ## but there is the MINIMUM of 1 day of infectiousness
  df_epi_info<-read.csv(file='../epi_fitting/epidemic_info/epidemic_info_net_onset_yest.csv')  
  df_epi_info$delay_isolation<-df_epi_info$delay_isolation + 2 ## Infectious from 2 days before
  N_timesteps<-max(df_epi_info$max_t_simulation)  ## NB it is 42
  ## (transition I->R done on infected at t-1)
  nrVertices <- network.size(HH_network_no_isolation)
  vertexNames <- HH_network_no_isolation %v% "vertex.names"  ## Needed to check isolation status
  infected <- rep(0, nrVertices)
  susceptible <- rep(1, nrVertices)
  ids_to_infec<-which(network::get.vertex.attribute(x = HH_network_no_isolation,"index")==1)
  infected[ids_to_infec] <- 1		
  susceptible[ids_to_infec] <- 0		
  infectedMat <- infected
  exposed <- rep(0, nrVertices)
  incidenceMat <- infected
  susceptibleMat<-susceptible
  # Individual attributes
  vertex_id <- 1:nrVertices
  timeSteps <- 1
  while(timeSteps<N_timesteps){ 
    hhfoi <- rep(0, nrVertices)
    # Expose each susceptible with a certain probability
    incidence <- rep(0, nrVertices)
    exposed_incidence <- rep(0, nrVertices)
    infectedPrev <- infected	# indicator of infection in the previous time step
    exposedPrev  <- exposed   # indicator of exposure in the previous time step
    susceptiblePrev<-susceptible
    # Transition from S---> E
    for(i in vertex_id[susceptible == 1]){
      vid<-vertexNames[i]
      t_isolation_this_vid<-df_epi_info$delay_isolation[which(df_epi_info$memid==vid)]
      #print(paste0("t=", timeSteps, "i= ",i," exposed= ", exposed[i]," susceptible= ",susceptible[i]))
      ### Choose network to consider depending on time of isolations
      if(timeSteps>t_isolation_this_vid){hh_contacts <- get.neighborhood(HH_network_w_isolation, i)}else{
        hh_contacts <- get.neighborhood(HH_network_no_isolation, i)
      }
      nrInfHH <- sum(infectedPrev[vertex_id %in% hh_contacts ], na.rm = TRUE)
      infProb <- 1-(1-beta1)^nrInfHH*(1-beta2)
      #comfoi[i] <- beta2
      SuccessExposure<-runif(1) < infProb
      if(SuccessExposure){
        exposed[i]<-1
        susceptible[i]<-0
        #   print(paste0("At time=",timeSteps , " id= ",i, "is removed from susceptibles"))
      }
      exposed_incidence[i]<-exposed[i]
    }
    susceptibleMat<-cbind(susceptibleMat, susceptible)
    # Transition from E---> I (TBD on exposed at previous time step)
    for(i in vertex_id[exposedPrev == 1]){
      SuccessInfection<-runif(1) < sigmaz
      if(SuccessInfection){
        infected[i]<-1
        exposed[i]<-0
        incidence[i] <- 1
      }
    }
    # Transition from I---> R (TBD on infected at previous time step)
    for(i in vertex_id[infectedPrev == 1]){
      SuccessRecover<-runif(1) < gammaz
      if(SuccessRecover){
        #infected[i]<-NA ## NAs mark who has eventually been infected
        infected[i]<-0 ## NAs mark who has eventually been infected
      }
    }
    # Append current status to matrix and update time variable
    infectedMat <- cbind(infectedMat, infected)
    incidenceMat <- cbind(incidenceMat, incidence)
    timeSteps <- timeSteps + 1    
  }
  
  infectedMat_censored<-infectedMat
  for(irow in 1:nrow(infectedMat_censored)){
    censoring_time_this_part<-df_epi_info$max_t_simulation[df_epi_info$memid==vertexNames[irow]]
    infectedMat_censored[irow,censoring_time_this_part:ncol(infectedMat_censored)]<-0
  }
  
  
  ## Compute whether a participant has been infected
  ever_infected<-base::rowSums(infectedMat_censored,na.rm = TRUE)
  ever_infected[ever_infected>0]<-1
  ## Dropping non-index
  ever_infected<-ever_infected[which(network::get.vertex.attribute(x = HH_network_no_isolation,"index")==0)]
  df_ever_infected<-data.frame(memid=vertexNames[which(network::get.vertex.attribute(x = HH_network_no_isolation,"index")==0)],ever_infected=ever_infected)
  
  
  ### Load information about hh size of the individual
  subject_info<-read.csv("../epi_simulation/subject_data_for_analysis_of_infection.csv")
  ## merge i lmemid, seleziona su household size and generate final df
  df_ever_infected<-merge(df_ever_infected,subject_info[,c("memid","hhmem_num")])
  df_ever_infected_by_hh<-df_ever_infected %>% group_by(hhmem_num) %>% summarize(tot_infected=sum(ever_infected))
  
  
  ### Load infection data
  df_data_to_fit<-read.csv("../epi_fitting/time_infection_matrix_net_onset_yest.csv",header=TRUE)
  colnames(df_data_to_fit)[1]<-"memid"
  dum<-as.matrix(df_data_to_fit[,2:dim(df_data_to_fit)[2]]) ## Dropping the Ids name column
  df_data_to_fit$ever_infected<-base::rowSums(dum,na.rm = TRUE) ## Changing to vector of ever infected
  df_data_to_fit$ever_infected[df_data_to_fit$ever_infected>0]<-1
  ## Dropping non-index
  df_data_to_fit<-df_data_to_fit %>% filter(memid %in% vertexNames[which(network::get.vertex.attribute(x = HH_network_no_isolation,"index")==0)])
  df_data_to_fit<-df_data_to_fit[,c("memid","ever_infected")]
  df_data_to_fit<-merge(df_data_to_fit,subject_info[,c("memid","hhmem_num")])
  ## Generate final data
  data_to_fit_by_hh<-df_data_to_fit %>% group_by(hhmem_num) %>% summarize(tot_infected=sum(ever_infected))
  
  
  ## Compute likelihood as Binomial log-likelihood for each hh_size
  error<- sum(log_lik(df_ever_infected_by_hh$tot_infected,data_to_fit_by_hh$tot_infected))
  return(list("Score"=-error))  ## Simple optim MINIMIZES by default
}
SEIR_simulate_two_networks_forSimpleOptim_Binloglik_Ncases_by_HH_averageNfunction<-function(full_list_params){
  Nsims<-100
  sum<-0
  counter<-0
  for(i_sim in 1:Nsims){
    dum<-unlist(SEIR_simulate_two_networks_forSimpleOptim_Binloglik_Ncases_by_HH(full_list_params))
    if(!is.na(dum)){
      sum<-sum+dum
        counter<-counter+1
    }
    
  }
  avg<-sum/counter
  return(list("Score"=as.numeric(avg)))
}




SEIR_simulate_two_networks_var_comm_trans_forBOptim_loglik_Ncases_by_HH_averageN<-function(beta1,beta2){
  Nsims<-50
  sum<-0
  for(i_sim in 1:Nsims){
    sum<-sum+unlist(SEIR_simulate_two_networks_var_comm_trans_forBOptim_loglik_Ncases_by_HH(beta1,beta2))
  }
  avg<-sum/Nsims
  return(list("Score"=avg))
}
SEIR_simulate_two_networks_var_comm_trans_forBOptim_loglik_Ncases_by_HH<-function(beta1,beta2){
  sigmaz<-1./1
  gammaz<-1./7    ### Average duration of infectivity is 7
  ### Read the comm incidence
  df_comm_transm<-read.csv("../epi_fitting/community_transmission.csv")
  ## but there is the MINIMUM of 1 day of infectiousness
  df_epi_info<-read.csv(file='../epi_fitting/epidemic_info/epidemic_info_net_onset_yest.csv')  
  df_epi_info$delay_isolation<-df_epi_info$delay_isolation + 2 ## Infectious from 2 days before
  ## (transition I->R done on infected at t-1)
  N_timesteps<-max(df_epi_info$max_t_simulation)  ## NB it is 42
  ### This is the maximum of variable "max_t_simulation"
  ## which is defined as 28 days (follow up) + time from cdc_primary_onset (time of symptom onset for index) to doe_dt (enrollment)
  ## So 28 + doe_dt - cdc_primary_onset
  #### NB: later, only infections happening before max_t_simulation are considered
  nrVertices <- network.size(HH_network_no_isolation)
  vertexNames <- HH_network_no_isolation %v% "vertex.names"  ## Needed to check isolation status
  infected <- rep(0, nrVertices)
  susceptible <- rep(1, nrVertices)
  ids_to_infec<-which(network::get.vertex.attribute(x = HH_network_no_isolation,"index")==1)
  infected[ids_to_infec] <- 1		
  susceptible[ids_to_infec] <- 0		
  infectedMat <- infected
  exposed <- rep(0, nrVertices)
  incidenceMat <- infected
  susceptibleMat<-susceptible
  # Individual attributes
  vertex_id <- 1:nrVertices
  timeSteps <- 1
  while(timeSteps<N_timesteps){ 
    hhfoi <- rep(0, nrVertices)
    # Expose each susceptible with a certain probability
    incidence <- rep(0, nrVertices)
    exposed_incidence <- rep(0, nrVertices)
    infectedPrev <- infected	# indicator of infection in the previous time step
    exposedPrev  <- exposed   # indicator of exposure in the previous time step
    susceptiblePrev<-susceptible
    # Transition from S---> E
    for(i in vertex_id[susceptible == 1]){
      vid<-vertexNames[i]
      t_isolation_this_vid<-df_epi_info$delay_isolation[which(df_epi_info$memid==vid)]
      #print(paste0("t=", timeSteps, "i= ",i," exposed= ", exposed[i]," susceptible= ",susceptible[i]))
      ### Choose network to consider depending on time of isolations
      if(timeSteps>t_isolation_this_vid){hh_contacts <- get.neighborhood(HH_network_w_isolation, i)}else{
        hh_contacts <- get.neighborhood(HH_network_no_isolation, i)
      }
      nrInfHH <- sum(infectedPrev[vertex_id %in% hh_contacts ], na.rm = TRUE)
      ### Match household to timing and site
      absolute_time<-as.Date(df_epi_info$cdc_primary_onsetdt[df_epi_info$memid==vid])+timeSteps
      site_of_id<-df_epi_info$site[df_epi_info$memid==vid]
      date_to_next_wednesday<-lubridate::ceiling_date(absolute_time,unit="weeks",week_start=3)
      ### debug: lubridate::wday(lubridate::ceiling_date(absolute_time,unit="weeks",week_start=3),label=TRUE)
      ### Get the relevant factor
      selection<-df_comm_transm$location==site_of_id& df_comm_transm$date==date_to_next_wednesday
      commModulationFactor<-df_comm_transm$ratio_weekly_cases[selection]
      #      print(paste0("For participant ",vid, " after ",timeSteps, " timesteps the date is ", absolute_time, " which is ceiled to ", date_to_next_wednesday, " for ", commModulationFactor ))
      infProb <- 1-(1-beta1)^nrInfHH*(1-beta2*commModulationFactor)
      #comfoi[i] <- beta2
      SuccessExposure<-runif(1) < infProb
      if(SuccessExposure){
        exposed[i]<-1
        susceptible[i]<-0
        #   print(paste0("At time=",timeSteps , " id= ",i, "is removed from susceptibles"))
      }
      exposed_incidence[i]<-exposed[i]
    }
    susceptibleMat<-cbind(susceptibleMat, susceptible)
    # Transition from E---> I (TBD on exposed at previous time step)
    for(i in vertex_id[exposedPrev == 1]){
      SuccessInfection<-runif(1) < sigmaz
      if(SuccessInfection){
        infected[i]<-1
        exposed[i]<-0
        incidence[i] <- 1
      }
    }
    # Transition from I---> R (TBD on infected at previous time step)
    for(i in vertex_id[infectedPrev == 1]){
      SuccessRecover<-runif(1) < gammaz
      if(SuccessRecover){
        #infected[i]<-NA ## NAs mark who has eventually been infected
        infected[i]<-0 ## NAs mark who has eventually been infected
      }
    }
    # Append current status to matrix and update time variable
    infectedMat <- cbind(infectedMat, infected)
    incidenceMat <- cbind(incidenceMat, incidence)
    timeSteps <- timeSteps + 1    
  }
  
  ## Ignoring infections happening after the max_time in the study
  infectedMat_censored<-infectedMat
  for(irow in 1:nrow(infectedMat_censored)){
    censoring_time_this_part<-df_epi_info$max_t_simulation[df_epi_info$memid==vertexNames[irow]]
    infectedMat_censored[irow,censoring_time_this_part:ncol(infectedMat_censored)]<-0
  }
  
  ## Compute whether a participant has been infected
  ever_infected<-base::rowSums(infectedMat_censored,na.rm = TRUE)
  ever_infected[ever_infected>0]<-1
  ## Dropping non-index
  ever_infected<-ever_infected[which(network::get.vertex.attribute(x = HH_network_no_isolation,"index")==0)]
  df_ever_infected<-data.frame(memid=vertexNames[which(network::get.vertex.attribute(x = HH_network_no_isolation,"index")==0)],ever_infected=ever_infected)
  
  
  ### Load information about hh size of the individual
  subject_info<-read.csv("../epi_simulation/subject_data_for_analysis_of_infection.csv")
  ## merge i lmemid, seleziona su household size and generate final df
  df_ever_infected<-merge(df_ever_infected,subject_info[,c("memid","hhmem_num")])
  df_ever_infected_by_hh<-df_ever_infected %>% group_by(hhmem_num) %>% summarize(tot_infected=sum(ever_infected))
  
  
  ### Load infection data
  df_data_to_fit<-read.csv("../epi_fitting/time_infection_matrix_net_onset_yest.csv",header=TRUE)
  colnames(df_data_to_fit)[1]<-"memid"
  dum<-as.matrix(df_data_to_fit[,2:dim(df_data_to_fit)[2]]) ## Dropping the Ids name column
  df_data_to_fit$ever_infected<-base::rowSums(dum,na.rm = TRUE) ## Changing to vector of ever infected
  df_data_to_fit$ever_infected[df_data_to_fit$ever_infected>0]<-1
  ## Dropping non-index
  df_data_to_fit<-df_data_to_fit %>% filter(memid %in% vertexNames[which(network::get.vertex.attribute(x = HH_network_no_isolation,"index")==0)])
  df_data_to_fit<-df_data_to_fit[,c("memid","ever_infected")]
  df_data_to_fit<-merge(df_data_to_fit,subject_info[,c("memid","hhmem_num")])
  ## Generate final data
  data_to_fit_by_hh<-df_data_to_fit %>% group_by(hhmem_num) %>% summarize(tot_infected=sum(ever_infected))
  ## Compute likelihood as Binomial log-likelihood for each hh_size
  error<- sum(log_lik(df_ever_infected_by_hh$tot_infected,data_to_fit_by_hh$tot_infected))
  return(list("Score"=error)) 
  
}


SEIR_simulate_two_networks_var_comm_trans_state_forBOptim_loglik_Ncases_by_HH_averageN<-function(beta1,beta2){
  Nsims<-10
  sum<-0
  for(i_sim in 1:Nsims){
    sum<-sum+unlist(SEIR_simulate_two_networks_var_comm_trans_state_forBOptim_loglik_Ncases_by_HH(beta1,beta2))
  }
  avg<-sum/Nsims
  return(list("Score"=avg))
}
SEIR_simulate_two_networks_var_comm_trans_state_forBOptim_loglik_Ncases_by_HH<-function(beta1,beta2){
  sigmaz<-1./1
  gammaz<-1./7    ### Average duration of infectivity is 7
  ### Read the comm incidence
  df_comm_transm<-read.csv("../epi_fitting/community_transmission_state.csv")
  ## but there is the MINIMUM of 1 day of infectiousness
  df_epi_info<-read.csv(file='../epi_fitting/epidemic_info/epidemic_info_net_onset_yest.csv')  
  df_epi_info$delay_isolation<-df_epi_info$delay_isolation + 2 ## Infectious from 2 days before
  ## (transition I->R done on infected at t-1)
  N_timesteps<-max(df_epi_info$max_t_simulation)  ## NB it is 42
  ### This is the maximum of variable "max_t_simulation"
  ## which is defined as 28 days (follow up) + time from cdc_primary_onset (time of symptom onset for index) to doe_dt (enrollment)
  ## So 28 + doe_dt - cdc_primary_onset
  #### NB: later, only infections happening before max_t_simulation are considered
  nrVertices <- network.size(HH_network_no_isolation)
  vertexNames <- HH_network_no_isolation %v% "vertex.names"  ## Needed to check isolation status
  infected <- rep(0, nrVertices)
  susceptible <- rep(1, nrVertices)
  ids_to_infec<-which(network::get.vertex.attribute(x = HH_network_no_isolation,"index")==1)
  infected[ids_to_infec] <- 1		
  susceptible[ids_to_infec] <- 0		
  infectedMat <- infected
  exposed <- rep(0, nrVertices)
  incidenceMat <- infected
  susceptibleMat<-susceptible
  # Individual attributes
  vertex_id <- 1:nrVertices
  timeSteps <- 1
  while(timeSteps<N_timesteps){ 
    hhfoi <- rep(0, nrVertices)
    # Expose each susceptible with a certain probability
    incidence <- rep(0, nrVertices)
    exposed_incidence <- rep(0, nrVertices)
    infectedPrev <- infected	# indicator of infection in the previous time step
    exposedPrev  <- exposed   # indicator of exposure in the previous time step
    susceptiblePrev<-susceptible
    # Transition from S---> E
    for(i in vertex_id[susceptible == 1]){
      vid<-vertexNames[i]
      t_isolation_this_vid<-df_epi_info$delay_isolation[which(df_epi_info$memid==vid)]
      #print(paste0("t=", timeSteps, "i= ",i," exposed= ", exposed[i]," susceptible= ",susceptible[i]))
      ### Choose network to consider depending on time of isolations
      if(timeSteps>t_isolation_this_vid){hh_contacts <- get.neighborhood(HH_network_w_isolation, i)}else{
        hh_contacts <- get.neighborhood(HH_network_no_isolation, i)
      }
      nrInfHH <- sum(infectedPrev[vertex_id %in% hh_contacts ], na.rm = TRUE)
      ### Match household to timing and site
      absolute_time<-as.Date(df_epi_info$cdc_primary_onsetdt[df_epi_info$memid==vid])+timeSteps
      site_of_id<-df_epi_info$site[df_epi_info$memid==vid]
      date_to_next_wednesday<-lubridate::ceiling_date(absolute_time,unit="weeks",week_start=3)
      ### debug: lubridate::wday(lubridate::ceiling_date(absolute_time,unit="weeks",week_start=3),label=TRUE)
      ### Get the relevant factor
      selection<-df_comm_transm$location==site_of_id& df_comm_transm$date==date_to_next_wednesday
      commModulationFactor<-df_comm_transm$ratio_weekly_cases[selection]
      #      print(paste0("For participant ",vid, " after ",timeSteps, " timesteps the date is ", absolute_time, " which is ceiled to ", date_to_next_wednesday, " for ", commModulationFactor ))
      infProb <- 1-(1-beta1)^nrInfHH*(1-beta2*commModulationFactor)
      #comfoi[i] <- beta2
      SuccessExposure<-runif(1) < infProb
      if(SuccessExposure){
        exposed[i]<-1
        susceptible[i]<-0
        #   print(paste0("At time=",timeSteps , " id= ",i, "is removed from susceptibles"))
      }
      exposed_incidence[i]<-exposed[i]
    }
    susceptibleMat<-cbind(susceptibleMat, susceptible)
    # Transition from E---> I (TBD on exposed at previous time step)
    for(i in vertex_id[exposedPrev == 1]){
      SuccessInfection<-runif(1) < sigmaz
      if(SuccessInfection){
        infected[i]<-1
        exposed[i]<-0
        incidence[i] <- 1
      }
    }
    # Transition from I---> R (TBD on infected at previous time step)
    for(i in vertex_id[infectedPrev == 1]){
      SuccessRecover<-runif(1) < gammaz
      if(SuccessRecover){
        #infected[i]<-NA ## NAs mark who has eventually been infected
        infected[i]<-0 ## NAs mark who has eventually been infected
      }
    }
    # Append current status to matrix and update time variable
    infectedMat <- cbind(infectedMat, infected)
    incidenceMat <- cbind(incidenceMat, incidence)
    timeSteps <- timeSteps + 1    
  }
  
  ## Ignoring infections happening after the max_time in the study
  infectedMat_censored<-infectedMat
  for(irow in 1:nrow(infectedMat_censored)){
    censoring_time_this_part<-df_epi_info$max_t_simulation[df_epi_info$memid==vertexNames[irow]]
    infectedMat_censored[irow,censoring_time_this_part:ncol(infectedMat_censored)]<-0
  }
  
  ## Compute whether a participant has been infected
  ever_infected<-base::rowSums(infectedMat_censored,na.rm = TRUE)
  ever_infected[ever_infected>0]<-1
  ## Dropping non-index
  ever_infected<-ever_infected[which(network::get.vertex.attribute(x = HH_network_no_isolation,"index")==0)]
  df_ever_infected<-data.frame(memid=vertexNames[which(network::get.vertex.attribute(x = HH_network_no_isolation,"index")==0)],ever_infected=ever_infected)
  
  
  ### Load information about hh size of the individual
  subject_info<-read.csv("../epi_simulation/subject_data_for_analysis_of_infection.csv")
  ## merge i lmemid, seleziona su household size and generate final df
  df_ever_infected<-merge(df_ever_infected,subject_info[,c("memid","hhmem_num")])
  df_ever_infected_by_hh<-df_ever_infected %>% group_by(hhmem_num) %>% summarize(tot_infected=sum(ever_infected))
  
  
  ### Load infection data
  df_data_to_fit<-read.csv("../epi_fitting/time_infection_matrix_net_onset_yest.csv",header=TRUE)
  colnames(df_data_to_fit)[1]<-"memid"
  dum<-as.matrix(df_data_to_fit[,2:dim(df_data_to_fit)[2]]) ## Dropping the Ids name column
  df_data_to_fit$ever_infected<-base::rowSums(dum,na.rm = TRUE) ## Changing to vector of ever infected
  df_data_to_fit$ever_infected[df_data_to_fit$ever_infected>0]<-1
  ## Dropping non-index
  df_data_to_fit<-df_data_to_fit %>% filter(memid %in% vertexNames[which(network::get.vertex.attribute(x = HH_network_no_isolation,"index")==0)])
  df_data_to_fit<-df_data_to_fit[,c("memid","ever_infected")]
  df_data_to_fit<-merge(df_data_to_fit,subject_info[,c("memid","hhmem_num")])
  ## Generate final data
  data_to_fit_by_hh<-df_data_to_fit %>% group_by(hhmem_num) %>% summarize(tot_infected=sum(ever_infected))
  ## Compute likelihood as Binomial log-likelihood for each hh_size
  error<- sum(log_lik(df_ever_infected_by_hh$tot_infected,data_to_fit_by_hh$tot_infected))
  return(list("Score"=error)) 
  
}

SEIR_simulate_two_networks_forBOptim_Binloglik_Ncases_by_HH<-function(beta1,beta2){
  sigmaz<-1./1.  ## 3 days to develop symptoms, but infectious from 2 days before symtpoms
  gammaz<-1./7    ### Average duration of infectivity is 7
  ## but there is the MINIMUM of 1 day of infectiousness
  df_epi_info<-read.csv(file='../epi_fitting/epidemic_info/epidemic_info_net_onset_yest.csv')  
  df_epi_info$delay_isolation<-df_epi_info$delay_isolation + 2 ## Infectious from 2 days before
  N_timesteps<-max(df_epi_info$max_t_simulation)  ## NB it is 42
  ## (transition I->R done on infected at t-1)
  nrVertices <- network.size(HH_network_no_isolation)
  vertexNames <- HH_network_no_isolation %v% "vertex.names"  ## Needed to check isolation status
  infected <- rep(0, nrVertices)
  susceptible <- rep(1, nrVertices)
  ids_to_infec<-which(network::get.vertex.attribute(x = HH_network_no_isolation,"index")==1)
  infected[ids_to_infec] <- 1		
  susceptible[ids_to_infec] <- 0		
  infectedMat <- infected
  exposed <- rep(0, nrVertices)
  incidenceMat <- infected
  susceptibleMat<-susceptible
  # Individual attributes
  vertex_id <- 1:nrVertices
  timeSteps <- 1
  while(timeSteps<N_timesteps){ 
    hhfoi <- rep(0, nrVertices)
    # Expose each susceptible with a certain probability
    incidence <- rep(0, nrVertices)
    exposed_incidence <- rep(0, nrVertices)
    infectedPrev <- infected	# indicator of infection in the previous time step
    exposedPrev  <- exposed   # indicator of exposure in the previous time step
    susceptiblePrev<-susceptible
    # Transition from S---> E
    for(i in vertex_id[susceptible == 1]){
      vid<-vertexNames[i]
      t_isolation_this_vid<-df_epi_info$delay_isolation[which(df_epi_info$memid==vid)]
      #print(paste0("t=", timeSteps, "i= ",i," exposed= ", exposed[i]," susceptible= ",susceptible[i]))
      ### Choose network to consider depending on time of isolations
      if(timeSteps>t_isolation_this_vid){hh_contacts <- get.neighborhood(HH_network_w_isolation, i)}else{
        hh_contacts <- get.neighborhood(HH_network_no_isolation, i)
      }
      nrInfHH <- sum(infectedPrev[vertex_id %in% hh_contacts ], na.rm = TRUE)
      infProb <- 1-(1-beta1)^nrInfHH*(1-beta2)
      #comfoi[i] <- beta2
      SuccessExposure<-runif(1) < infProb
      if(SuccessExposure){
        exposed[i]<-1
        susceptible[i]<-0
        #   print(paste0("At time=",timeSteps , " id= ",i, "is removed from susceptibles"))
      }
      exposed_incidence[i]<-exposed[i]
    }
    susceptibleMat<-cbind(susceptibleMat, susceptible)
    # Transition from E---> I (TBD on exposed at previous time step)
    for(i in vertex_id[exposedPrev == 1]){
      SuccessInfection<-runif(1) < sigmaz
      if(SuccessInfection){
        infected[i]<-1
        exposed[i]<-0
        incidence[i] <- 1
      }
    }
    # Transition from I---> R (TBD on infected at previous time step)
    for(i in vertex_id[infectedPrev == 1]){
      SuccessRecover<-runif(1) < gammaz
      if(SuccessRecover){
        #infected[i]<-NA ## NAs mark who has eventually been infected
        infected[i]<-0 ## NAs mark who has eventually been infected
      }
    }
    # Append current status to matrix and update time variable
    infectedMat <- cbind(infectedMat, infected)
    incidenceMat <- cbind(incidenceMat, incidence)
    timeSteps <- timeSteps + 1    
  }
  
  infectedMat_censored<-infectedMat
  for(irow in 1:nrow(infectedMat_censored)){
    censoring_time_this_part<-df_epi_info$max_t_simulation[df_epi_info$memid==vertexNames[irow]]
    infectedMat_censored[irow,censoring_time_this_part:ncol(infectedMat_censored)]<-0
  }
  
  
  ## Compute whether a participant has been infected
  ever_infected<-base::rowSums(infectedMat_censored,na.rm = TRUE)
  ever_infected[ever_infected>0]<-1
  ## Dropping non-index
  ever_infected<-ever_infected[which(network::get.vertex.attribute(x = HH_network_no_isolation,"index")==0)]
  df_ever_infected<-data.frame(memid=vertexNames[which(network::get.vertex.attribute(x = HH_network_no_isolation,"index")==0)],ever_infected=ever_infected)
  
  
  ### Load information about hh size of the individual
  subject_info<-read.csv("../epi_simulation/subject_data_for_analysis_of_infection.csv")
  ## merge i lmemid, seleziona su household size and generate final df
  df_ever_infected<-merge(df_ever_infected,subject_info[,c("memid","hhmem_num")])
  df_ever_infected_by_hh<-df_ever_infected %>% group_by(hhmem_num) %>% summarize(tot_infected=sum(ever_infected))
  
  
  ### Load infection data
  df_data_to_fit<-read.csv("../epi_fitting/time_infection_matrix_net_onset_yest.csv",header=TRUE)
  colnames(df_data_to_fit)[1]<-"memid"
  dum<-as.matrix(df_data_to_fit[,2:dim(df_data_to_fit)[2]]) ## Dropping the Ids name column
  df_data_to_fit$ever_infected<-base::rowSums(dum,na.rm = TRUE) ## Changing to vector of ever infected
  df_data_to_fit$ever_infected[df_data_to_fit$ever_infected>0]<-1
  ## Dropping non-index
  df_data_to_fit<-df_data_to_fit %>% filter(memid %in% vertexNames[which(network::get.vertex.attribute(x = HH_network_no_isolation,"index")==0)])
  df_data_to_fit<-df_data_to_fit[,c("memid","ever_infected")]
  df_data_to_fit<-merge(df_data_to_fit,subject_info[,c("memid","hhmem_num")])
  ## Generate final data
  data_to_fit_by_hh<-df_data_to_fit %>% group_by(hhmem_num) %>% summarize(tot_infected=sum(ever_infected))
  
  
  ## Compute likelihood as Binomial log-likelihood for each hh_size
  error<- sum(log_lik(df_ever_infected_by_hh$tot_infected,data_to_fit_by_hh$tot_infected))
  return(list("Score"=error)) 
}

SEIR_simulate_two_networks_forBOptim_Binloglik_Ncases_by_HH_averageN20<-function(beta1,beta2){
  Nsims<-20
  sum<-0
  for(i_sim in 1:Nsims){
    sum<-sum+unlist(SEIR_simulate_two_networks_forBOptim_Binloglik_Ncases_by_HH(beta1,beta2))
  }
  avg<-sum/Nsims
  return(list("Score"=as.numeric(avg)))
}

SEIR_simulate_two_networks_forBOptim_Binloglik_Ncases_by_HH_averageN50<-function(beta1,beta2){
  Nsims<-50
  sum<-0
  for(i_sim in 1:Nsims){
    sum<-sum+unlist(SEIR_simulate_two_networks_forBOptim_Binloglik_Ncases_by_HH(beta1,beta2))
  }
  avg<-sum/Nsims
  return(list("Score"=as.numeric(avg)))
}

SEIR_simulate_two_networks_forBOptim_Binloglik_Ncases_by_HH_averageN200<-function(beta1,beta2){
  Nsims<-200
  sum<-0
  for(i_sim in 1:Nsims){
    sum<-sum+unlist(SEIR_simulate_two_networks_forBOptim_Binloglik_Ncases_by_HH(beta1,beta2))
  }
  avg<-sum/Nsims
  return(list("Score"=as.numeric(avg)))
}


SEIR_simulate_two_networks_forBOptim_Binloglik_Ncases_by_HH_no2_averageN50<-function(beta1,beta2){
  Nsims<-50
  sum<-0
  for(i_sim in 1:Nsims){
    sum<-sum+unlist(SEIR_simulate_two_networks_forBOptim_Binloglik_Ncases_by_HH_no2(beta1,beta2))
  }
  avg<-sum/Nsims
  return(list("Score"=as.numeric(avg)))
}


SEIR_simulate_two_networks_forBOptim_LS_Ncases_by_HH<-function(beta1,beta2){
  sigmaz<-1./1.  ## 3 days to develop symptoms, but infectious from 2 days before symtpoms
  gammaz<-1./7    ### Average duration of infectivity is 7
  ## but there is the MINIMUM of 1 day of infectiousness
  df_epi_info<-read.csv(file='../epi_fitting/epidemic_info/epidemic_info_net_onset_yest.csv')  
  df_epi_info$delay_isolation<-df_epi_info$delay_isolation + 2 ## Infectious from 2 days before
  N_timesteps<-max(df_epi_info$max_t_simulation)  ## NB it is 42
  ## (transition I->R done on infected at t-1)
  nrVertices <- network.size(HH_network_no_isolation)
  vertexNames <- HH_network_no_isolation %v% "vertex.names"  ## Needed to check isolation status
  infected <- rep(0, nrVertices)
  susceptible <- rep(1, nrVertices)
  ids_to_infec<-which(network::get.vertex.attribute(x = HH_network_no_isolation,"index")==1)
  infected[ids_to_infec] <- 1		
  susceptible[ids_to_infec] <- 0		
  infectedMat <- infected
  exposed <- rep(0, nrVertices)
  incidenceMat <- infected
  susceptibleMat<-susceptible
  # Individual attributes
  vertex_id <- 1:nrVertices
  timeSteps <- 1
  while(timeSteps<N_timesteps){ 
    hhfoi <- rep(0, nrVertices)
    # Expose each susceptible with a certain probability
    incidence <- rep(0, nrVertices)
    exposed_incidence <- rep(0, nrVertices)
    infectedPrev <- infected	# indicator of infection in the previous time step
    exposedPrev  <- exposed   # indicator of exposure in the previous time step
    susceptiblePrev<-susceptible
    # Transition from S---> E
    for(i in vertex_id[susceptible == 1]){
      vid<-vertexNames[i]
      t_isolation_this_vid<-df_epi_info$delay_isolation[which(df_epi_info$memid==vid)]
      #print(paste0("t=", timeSteps, "i= ",i," exposed= ", exposed[i]," susceptible= ",susceptible[i]))
      ### Choose network to consider depending on time of isolations
      if(timeSteps>t_isolation_this_vid){hh_contacts <- get.neighborhood(HH_network_w_isolation, i)}else{
        hh_contacts <- get.neighborhood(HH_network_no_isolation, i)
      }
      nrInfHH <- sum(infectedPrev[vertex_id %in% hh_contacts ], na.rm = TRUE)
      infProb <- 1-(1-beta1)^nrInfHH*(1-beta2)
      #comfoi[i] <- beta2
      SuccessExposure<-runif(1) < infProb
      if(SuccessExposure){
        exposed[i]<-1
        susceptible[i]<-0
        #   print(paste0("At time=",timeSteps , " id= ",i, "is removed from susceptibles"))
      }
      exposed_incidence[i]<-exposed[i]
    }
    susceptibleMat<-cbind(susceptibleMat, susceptible)
    # Transition from E---> I (TBD on exposed at previous time step)
    for(i in vertex_id[exposedPrev == 1]){
      SuccessInfection<-runif(1) < sigmaz
      if(SuccessInfection){
        infected[i]<-1
        exposed[i]<-0
        incidence[i] <- 1
      }
    }
    # Transition from I---> R (TBD on infected at previous time step)
    for(i in vertex_id[infectedPrev == 1]){
      SuccessRecover<-runif(1) < gammaz
      if(SuccessRecover){
        #infected[i]<-NA ## NAs mark who has eventually been infected
        infected[i]<-0 ## NAs mark who has eventually been infected
      }
    }
    # Append current status to matrix and update time variable
    infectedMat <- cbind(infectedMat, infected)
    incidenceMat <- cbind(incidenceMat, incidence)
    timeSteps <- timeSteps + 1    
  }
  
  infectedMat_censored<-infectedMat
  for(irow in 1:nrow(infectedMat_censored)){
    censoring_time_this_part<-df_epi_info$max_t_simulation[df_epi_info$memid==vertexNames[irow]]
    infectedMat_censored[irow,censoring_time_this_part:ncol(infectedMat_censored)]<-0
  }
  
  
  ## Compute whether a participant has been infected
  ever_infected<-base::rowSums(infectedMat_censored,na.rm = TRUE)
  ever_infected[ever_infected>0]<-1
  ## Dropping non-index
  ever_infected<-ever_infected[which(network::get.vertex.attribute(x = HH_network_no_isolation,"index")==0)]
  df_ever_infected<-data.frame(memid=vertexNames[which(network::get.vertex.attribute(x = HH_network_no_isolation,"index")==0)],ever_infected=ever_infected)
  
  
  ### Load information about hh size of the individual
  subject_info<-read.csv("../epi_simulation/subject_data_for_analysis_of_infection.csv")
  ## merge i lmemid, seleziona su household size and generate final df
  df_ever_infected<-merge(df_ever_infected,subject_info[,c("memid","hhmem_num")])
  df_ever_infected_by_hh<-df_ever_infected %>% group_by(hhmem_num) %>% summarize(tot_infected=sum(ever_infected))
  
  
  ### Load infection data
  df_data_to_fit<-read.csv("../epi_fitting/time_infection_matrix_net_onset_yest.csv",header=TRUE)
  colnames(df_data_to_fit)[1]<-"memid"
  dum<-as.matrix(df_data_to_fit[,2:dim(df_data_to_fit)[2]]) ## Dropping the Ids name column
  df_data_to_fit$ever_infected<-base::rowSums(dum,na.rm = TRUE) ## Changing to vector of ever infected
  df_data_to_fit$ever_infected[df_data_to_fit$ever_infected>0]<-1
  ## Dropping non-index
  df_data_to_fit<-df_data_to_fit %>% filter(memid %in% vertexNames[which(network::get.vertex.attribute(x = HH_network_no_isolation,"index")==0)])
  df_data_to_fit<-df_data_to_fit[,c("memid","ever_infected")]
  df_data_to_fit<-merge(df_data_to_fit,subject_info[,c("memid","hhmem_num")])
  ## Generate final data
  data_to_fit_by_hh<-df_data_to_fit %>% group_by(hhmem_num) %>% summarize(tot_infected=sum(ever_infected))
  
  
  ## Compute likelihood as Binomial log-likelihood for each hh_size
  error<- sum(minusLS(df_ever_infected_by_hh$tot_infected,data_to_fit_by_hh$tot_infected))
  return(list("Score"=error)) 
}
SEIR_simulate_two_networks_forBOptim_LS_Ncases_by_HH_averageN<-function(beta1,beta2){
  Nsims<-20
  sum<-0
  for(i_sim in 1:Nsims){
    sum<-sum+unlist(SEIR_simulate_two_networks_forBOptim_LS_Ncases_by_HH(beta1,beta2))
  }
  avg<-sum/Nsims
  return(list("Score"=avg))
}




SEIR_simulate_two_networks_var_infectivity_forBOptim_loglik_Ncases_by_HH_averageN<-function(beta1,beta2){
  Nsims<-50
  sum<-0
  for(i_sim in 1:Nsims){
    sum<-sum+unlist(SEIR_simulate_two_networks_var_infectivity_forBOptim_loglik_Ncases_by_HH(beta1,beta2))
  }
  avg<-sum/Nsims
  return(list("Score"=avg))
}




SEIR_simulate_two_networks_var_infectivity_forBOptim_loglik_Ncases_by_HH<-function(beta1,beta2){
  sigmaz<-1./1
  gammaz<-1./7    ### Average duration of infectivity is 7
  ### Read the comm incidence
  df_comm_transm<-read.csv("../epi_fitting/community_transmission.csv")
  ## but there is the MINIMUM of 1 day of infectiousness
  df_epi_info<-read.csv(file='../epi_fitting/epidemic_info/epidemic_info_net_onset_yest.csv')  
  df_epi_info$delay_isolation<-df_epi_info$delay_isolation + 2 ## Infectious from 2 days before
  
  df_infectiousness<-read.csv("../epi_fitting/infectiousness_over_time.csv")
  dum_df<-data.frame(time=(12:100), infectuosness=rep(0,89))
  df_infectiousness<-rbind(df_infectiousness,dum_df)  ## File defined up to 12
  ## (transition I->R done on infected at t-1)
  N_timesteps<-max(df_epi_info$max_t_simulation)  ## NB it is 42
  ### This is the maximum of variable "max_t_simulation"
  ## which is defined as 28 days (follow up) + time from cdc_primary_onset (time of symptom onset for index) to doe_dt (enrollment)
  ## So 28 + doe_dt - cdc_primary_onset
  #### NB: later, only infections happening before max_t_simulation are considered
  nrVertices <- network.size(HH_network_no_isolation)
  vertexNames <- HH_network_no_isolation %v% "vertex.names"  ## Needed to check isolation status
  infected <- rep(0, nrVertices)
  time_of_infection<-rep(NA, nrVertices)
  time_of_recovery<-rep(NA, nrVertices)
  susceptible <- rep(1, nrVertices)
  ids_to_infec<-which(network::get.vertex.attribute(x = HH_network_no_isolation,"index")==1)
  infected[ids_to_infec] <- 1
  time_of_infection[ids_to_infec]<-0
  susceptible[ids_to_infec] <- 0		
  infectedMat <- infected
  exposed <- rep(0, nrVertices)
  incidenceMat <- infected
  susceptibleMat<-susceptible
  # Individual attributes
  vertex_id <- 1:nrVertices
  timeSteps <- 1
  while(timeSteps<N_timesteps){ 
    hhfoi <- rep(0, nrVertices)
    # Expose each susceptible with a certain probability
    incidence <- rep(0, nrVertices)
    exposed_incidence <- rep(0, nrVertices)
    infectedPrev <- infected	# indicator of infection in the previous time step
    exposedPrev  <- exposed   # indicator of exposure in the previous time step
    susceptiblePrev<-susceptible
    # Transition from S---> E
    for(i in vertex_id[susceptible == 1]){
      vid<-vertexNames[i]
      t_isolation_this_vid<-df_epi_info$delay_isolation[which(df_epi_info$memid==vid)]
      #print(paste0("t=", timeSteps, "i= ",i," exposed= ", exposed[i]," susceptible= ",susceptible[i]))
      ### Choose network to consider depending on time of isolations
      if(timeSteps>t_isolation_this_vid){hh_contacts <- get.neighborhood(HH_network_w_isolation, i)}else{
        hh_contacts <- get.neighborhood(HH_network_no_isolation, i)
      }
      #nrInfHH <- sum(infectedPrev[vertex_id %in% hh_contacts ], na.rm = TRUE)
      ## Loop on hh_contacts and sum infection parameters
      prod_escape_infection_hh<-1
      if(length(hh_contacts)>0){
      for(i_hh_member in 1:length(hh_contacts)){
        this_contact_id<-hh_contacts[i_hh_member]
        if(infectedPrev[this_contact_id]==1){
          time_from_infection= timeSteps-time_of_infection[this_contact_id]
          
          factor<-df_infectiousness$infectuosness[time_from_infection]
          contribution<-(1-beta1*factor)
          prod_escape_infection_hh<-prod_escape_infection_hh*contribution
          #print(paste0("contact with ",this_contact_id,"t_from_infection= ",time_from_infection," infection prob= ",infProb))
        }
      }
      }
      ## For each element 1-beta1*factor
      infProb <- 1-prod_escape_infection_hh*(1-beta2)
      #print(paste0("Escape infection= ",prod_escape_infection_hh," infection prob= ",infProb))
      #comfoi[i] <- beta2
      SuccessExposure<-runif(1) < infProb
      if(SuccessExposure){
        exposed[i]<-1
        susceptible[i]<-0
        #   print(paste0("At time=",timeSteps , " id= ",i, "is removed from susceptibles"))
      }
      exposed_incidence[i]<-exposed[i]
    }
    susceptibleMat<-cbind(susceptibleMat, susceptible)
    # Transition from E---> I (TBD on exposed at previous time step)
    for(i in vertex_id[exposedPrev == 1]){
      SuccessInfection<-runif(1) < sigmaz
      if(SuccessInfection){
        infected[i]<-1
        exposed[i]<-0
        incidence[i] <- 1
        time_of_infection[i]<-timeSteps
        
      }
    }
    # Transition from I---> R (TBD on infected at previous time step)
    for(i in vertex_id[infectedPrev == 1]){
      SuccessRecover<-runif(1) < gammaz
      if(SuccessRecover){
        #infected[i]<-NA ## NAs mark who has eventually been infected
        infected[i]<-0 
        time_of_recovery[i]<-timeSteps
      }
    }
    # Append current status to matrix and update time variable
    infectedMat <- cbind(infectedMat, infected)
    incidenceMat <- cbind(incidenceMat, incidence)
    timeSteps <- timeSteps + 1    
  }
  
  ## Ignoring infections happening after the max_time in the study
  infectedMat_censored<-infectedMat
  for(irow in 1:nrow(infectedMat_censored)){
    censoring_time_this_part<-df_epi_info$max_t_simulation[df_epi_info$memid==vertexNames[irow]]
    infectedMat_censored[irow,censoring_time_this_part:ncol(infectedMat_censored)]<-0
  }
  
  ## Compute whether a participant has been infected
  ever_infected<-base::rowSums(infectedMat_censored,na.rm = TRUE)
  ever_infected[ever_infected>0]<-1
  ## Dropping non-index
  ever_infected<-ever_infected[which(network::get.vertex.attribute(x = HH_network_no_isolation,"index")==0)]
  df_ever_infected<-data.frame(memid=vertexNames[which(network::get.vertex.attribute(x = HH_network_no_isolation,"index")==0)],ever_infected=ever_infected)
  
  
  ### Load information about hh size of the individual
  subject_info<-read.csv("../epi_simulation/subject_data_for_analysis_of_infection.csv")
  ## merge i lmemid, seleziona su household size and generate final df
  df_ever_infected<-merge(df_ever_infected,subject_info[,c("memid","hhmem_num")])
  df_ever_infected_by_hh<-df_ever_infected %>% group_by(hhmem_num) %>% summarize(tot_infected=sum(ever_infected))
  
  
  ### Load infection data
  df_data_to_fit<-read.csv("../epi_fitting/time_infection_matrix_net_onset_yest.csv",header=TRUE)
  colnames(df_data_to_fit)[1]<-"memid"
  dum<-as.matrix(df_data_to_fit[,2:dim(df_data_to_fit)[2]]) ## Dropping the Ids name column
  df_data_to_fit$ever_infected<-base::rowSums(dum,na.rm = TRUE) ## Changing to vector of ever infected
  df_data_to_fit$ever_infected[df_data_to_fit$ever_infected>0]<-1
  ## Dropping non-index
  df_data_to_fit<-df_data_to_fit %>% filter(memid %in% vertexNames[which(network::get.vertex.attribute(x = HH_network_no_isolation,"index")==0)])
  df_data_to_fit<-df_data_to_fit[,c("memid","ever_infected")]
  df_data_to_fit<-merge(df_data_to_fit,subject_info[,c("memid","hhmem_num")])
  ## Generate final data
  data_to_fit_by_hh<-df_data_to_fit %>% group_by(hhmem_num) %>% summarize(tot_infected=sum(ever_infected))
  ## Compute likelihood as Binomial log-likelihood for each hh_size
  error<- sum(log_lik(df_ever_infected_by_hh$tot_infected,data_to_fit_by_hh$tot_infected))
  return(list("Score"=error)) 
  
}

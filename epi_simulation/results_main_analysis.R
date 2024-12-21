# Rcode household contact networks
# Last update by Pietro on 02/06/2023
# Last run by Pietro on 29/11/2024
  library(igraph)
  library(ergm)
  library(ergm.multi)
  library(intergraph)
  library(sjPlot)
  library(ggpubr)
  library(sjlabelled)
  library(sjmisc)
  library(dplyr)
  library(modelsummary)
  library(statnet)
  library(ParBayesianOptimization)
  library(tidyverse)
  
  rm(list = ls())
  ########################################################
  `%notin%` <- function(a,b) ! a %in% b  
  ### Automatically set working directory
  if(require(rstudioapi) && isAvailable()){
    current_path <- getActiveDocumentContext()$path 
    setwd(dirname(current_path ))
  }
  
  
  source("../sim_libs/sim_models.R")

  fit_tag<-"main_analysis"
  optObj<-readRDS(paste0("../epi_fitting/best_params_final_switch_isolation_",fit_tag,".RDS"))
  onset_network<-readRDS("../network_obj_generation/net_objects/graph_full_dataset_onset_phys_in_common_yest_onlypart_symmetric.RDS")
  yest_network<-readRDS("../network_obj_generation/net_objects/graph_full_dataset_yest_phys_in_common_onset_onlypart_symmetric.RDS")
    
  SIMULATE_NORM_ISO<-TRUE # Set to FALSE if simualation already done #(contacts: onset/yest)
  SIMULATE_PRE_ISO<-TRUE  # Set to FALSE if simualation already done #(contacts: yest/yest)
  SIMULATE_NO_ISO<-TRUE   # Set to FALSE if simualation already done #(contacts: onset/onset)
  N_simulations<-1000
  Net_per_batch<-10 ## this should be the same as in the network simulation ("simulate_networks.R")
  N_batches<-floor(N_simulations/Net_per_batch)
  best_params<-getBestPars(optObj)
  
  ## Load best fit values from profile likelihood rather than from optObj
  best_params$beta1<-0.1617169
  best_params$beta2<-0.001051652
  
  
  set.seed(24112022)
  if(SIMULATE_NORM_ISO){
    if(exists("normal_isolation_inf_df")){rm(normal_isolation_inf_df)}
    batch_counter<-1
    for(batch_counter in 1:N_batches){
      onset_networks<-readRDS(paste0("../network_fitting/simulated_networks/onset_batch_", batch_counter,".RDS"))
      yest_networks<-readRDS(paste0("../network_fitting/simulated_networks/yest_batch_", batch_counter,".RDS"))
      print(paste0("batch ",batch_counter," of ",N_batches ))
      for(counter in 1:Net_per_batch){
        print(paste0("    ",counter))
        dum_list<-SEIR_simulate_two_networks(onset_networks[[counter]],
                                             yest_networks[[counter]],
                                             beta1=best_params[[1]],
                                             beta2=best_params[[2]]
        )
        dum_list$n<-(batch_counter-1)*Net_per_batch+counter
        if(!exists("normal_isolation_inf_df")){normal_isolation_inf_df<-dum_list}else{normal_isolation_inf_df<-rbind(normal_isolation_inf_df,dum_list)}
      }
      rm(onset_networks)
      rm(yest_networks)
    }
    write.csv(normal_isolation_inf_df,paste0("simulation_normal_isolation_infection_by_memid_",fit_tag ,".csv"),row.names = FALSE)
  }
  set.seed(24112022)
  if(SIMULATE_PRE_ISO){
    if(exists("preventive_isolation_df")){rm(preventive_isolation_df)}
    batch_counter<-1
    for(batch_counter in 1:N_batches){
      yest_networks<-readRDS(paste0("../network_fitting/simulated_networks/yest_batch_", batch_counter,".RDS"))
      print(paste0("batch ",batch_counter," of ",N_batches ))
      counter<-1
      for(counter in 1:Net_per_batch){
        print(paste0("    ",counter))
        dum_list<-SEIR_simulate_two_networks(yest_networks[[counter]],
                                   yest_networks[[counter]],
                                   beta1=best_params[[1]],
                                   beta2=best_params[[2]])
        dum_list$n<-(batch_counter-1)*Net_per_batch+counter
      if(!exists("preventive_isolation_df")){preventive_isolation_df<-dum_list}else{preventive_isolation_df<-rbind(preventive_isolation_df,dum_list)}
      }
      rm(yest_networks)
    }
    write.csv(preventive_isolation_df,paste0("simulation_preventive_isolation_infection_by_memid_",fit_tag ,".csv"),row.names = FALSE)
  }
  set.seed(24112022)
  if(SIMULATE_NO_ISO){
    if(exists("no_isolation_df")){rm(no_isolation_df)}
    batch_counter<-1
    for(batch_counter in 1:N_batches){
      onset_networks<-readRDS(paste0("../network_fitting/simulated_networks/onset_batch_", batch_counter,".RDS"))
      print(paste0("batch ",batch_counter," of ",N_batches ))
      counter<-1
      for(counter in 1:Net_per_batch){
        print(paste0("    ",counter))
        dum_list<-SEIR_simulate_two_networks(onset_networks[[counter]],
                                             onset_networks[[counter]],
                                             beta1=best_params[[1]],
                                             beta2=best_params[[2]])
        dum_list$n<-(batch_counter-1)*Net_per_batch+counter
        if(!exists("no_isolation_df")){no_isolation_df<-dum_list}else{no_isolation_df<-rbind(no_isolation_df,dum_list)}
      }
      rm(onset_networks)
    }
    write.csv(no_isolation_df,paste0("simulation_no_isolation_infection_by_memid_", fit_tag,".csv"),row.names = FALSE)
  }


subject_info<-read.csv("../epi_simulation/subject_data_for_analysis_of_infection.csv")

data_to_fit<-read.csv("../epi_fitting/time_infection_matrix_net_onset_yest.csv",header=TRUE)
memids<-data_to_fit$X
data_to_fit<-as.matrix(data_to_fit[,2:dim(data_to_fit)[2]]) ## Dropping the Ids name column
data_to_fit<-base::rowSums(data_to_fit,na.rm = TRUE) ## Changing to vector of ever infected
data_to_fit[data_to_fit>0]<-1
df_ever_infected<-data.frame(memid=memids,ever_infected=data_to_fit)
df_ever_infected<-merge(subject_info[,c("memid","age","hhmem_num","index")],df_ever_infected,by="memid")

normal_isolation_inf_df<-read.csv(paste0("simulation_normal_isolation_infection_by_memid_",fit_tag,".csv"))
preventive_isolation_df<-read.csv(paste0("simulation_preventive_isolation_infection_by_memid_",fit_tag,".csv"))
no_isolation_df<-read.csv(paste0("simulation_no_isolation_infection_by_memid_", fit_tag,".csv"))

#colnames(subject_info)

normal_isolation_inf_df<-merge(subject_info[,c("memid","age","hhmem_num")],normal_isolation_inf_df,by="memid")
preventive_isolation_df<-merge(subject_info[,c("memid","age","hhmem_num")],preventive_isolation_df,by="memid")
no_isolation_df<-merge(subject_info[,c("memid","age","hhmem_num")],no_isolation_df,by="memid")
      
      
### Filtering indexes out      
normal_isolation_inf_df<-normal_isolation_inf_df %>% filter(index==0)
preventive_isolation_df<-preventive_isolation_df %>% filter(index==0)
no_isolation_df<-no_isolation_df %>% filter(index==0)
df_ever_infected<-df_ever_infected %>% filter(index==0)




normal_isolation_inf_df$hhmem_num_cat<-normal_isolation_inf_df$hhmem_num
normal_isolation_inf_df$hhmem_num_cat[normal_isolation_inf_df$hhmem_num_cat>4]<-"5+"
normal_isolation_inf_df$hhmem_num_cat<-paste0("HH size ",normal_isolation_inf_df$hhmem_num_cat)

preventive_isolation_df$hhmem_num_cat<-preventive_isolation_df$hhmem_num
preventive_isolation_df$hhmem_num_cat[preventive_isolation_df$hhmem_num_cat>4]<-"5+"
preventive_isolation_df$hhmem_num_cat<-paste0("HH size ",preventive_isolation_df$hhmem_num_cat)

no_isolation_df$hhmem_num_cat<-no_isolation_df$hhmem_num
no_isolation_df$hhmem_num_cat[no_isolation_df$hhmem_num_cat>4]<-"5+"
no_isolation_df$hhmem_num_cat<-paste0("HH size ",no_isolation_df$hhmem_num_cat)

df_ever_infected$hhmem_num_cat<-df_ever_infected$hhmem_num
df_ever_infected$hhmem_num_cat[df_ever_infected$hhmem_num_cat>4]<-"5+"
df_ever_infected$hhmem_num_cat<-paste0("HH size ",df_ever_infected$hhmem_num_cat)

### Gropuping by hh_size_cat
N_infected_norm_iso_by_hh_cat<-normal_isolation_inf_df %>% group_by(n,hhmem_num_cat) %>% 
  summarise(Sum_infected = sum(ever_infected),Perc_infected=sum(ever_infected)/length(ever_infected),Potential_infected=length(ever_infected))
N_infected_norm_iso_by_hh_cat$net<-"Isolation after test"
N_infected_pre_iso_by_hh_cat<-preventive_isolation_df %>% group_by(n,hhmem_num_cat) %>% 
  summarise(Sum_infected = sum(ever_infected),Perc_infected=sum(ever_infected)/length(ever_infected),Potential_infected=length(ever_infected))
N_infected_pre_iso_by_hh_cat$net<-"Isolation after symptoms"
N_infected_no_iso_by_hh_cat<-no_isolation_df %>% group_by(n,hhmem_num_cat) %>% 
  summarise(Sum_infected = sum(ever_infected),Perc_infected=sum(ever_infected)/length(ever_infected),Potential_infected=length(ever_infected))
N_infected_no_iso_by_hh_cat$net<-"No isolation"
### Rbidnding
N_infected_by_hh_cat<-rbind(N_infected_norm_iso_by_hh_cat,N_infected_pre_iso_by_hh_cat,N_infected_no_iso_by_hh_cat)
### Adding data from FluTes
data_N_infected_by_hh_cat<-df_ever_infected %>% group_by(hhmem_num_cat) %>% 
  summarise(Sum_infected = sum(ever_infected),Perc_infected=sum(ever_infected)/length(ever_infected),Potential_infected=length(ever_infected))
N_infected_by_hh_cat$net<-factor(x=N_infected_by_hh_cat$net,levels=c("No isolation","Isolation after test","Isolation after symptoms"))
N_infected_by_hh_cat%>% group_by(net,hhmem_num_cat) %>% summarize(mean=mean(Sum_infected),sd=sd(Sum_infected))


### stat test for hh size 3
a<-N_infected_by_hh_cat %>% filter(hhmem_num_cat=="HH size 3",net=="Isolation after test")
b<-N_infected_by_hh_cat %>% filter(hhmem_num_cat=="HH size 3",net=="Isolation after symptoms")
ks.test(a$Sum_infected,b$Sum_infected)
t.test(a$Sum_infected,b$Sum_infected)
print(1-mean(b$Sum_infected)/mean(a$Sum_infected))
# ### stat test for hh size 4
a<-N_infected_by_hh_cat %>% filter(hhmem_num_cat=="HH size 4",net=="Isolation after test")
b<-N_infected_by_hh_cat %>% filter(hhmem_num_cat=="HH size 4",net=="Isolation after symptoms")
ks.test(a$Sum_infected,b$Sum_infected)
t.test(a$Sum_infected,b$Sum_infected)
print(1-mean(b$Sum_infected)/mean(a$Sum_infected))
# 
# ### stat test for hh size 5
a<-N_infected_by_hh_cat %>% filter(hhmem_num_cat=="HH size 5+",net=="Isolation after test")
b<-N_infected_by_hh_cat %>% filter(hhmem_num_cat=="HH size 5+",net=="Isolation after symptoms")
# ks.test(a$Sum_infected,b$Sum_infected)
t.test(a$Sum_infected,b$Sum_infected)
print(1-mean(b$Sum_infected)/mean(a$Sum_infected))

ggplot(data=N_infected_by_hh_cat)+geom_boxplot(outlier.size = -1,mapping = aes(y=Perc_infected,x=net))+ylab("Secondary infections")+xlab("Scenario")+
  facet_wrap(~hhmem_num_cat)+theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))#+geom_point(data=data_N_infected_by_hh_cat,mapping=aes(x="Isolation after test",y=Perc_infected),color="red")

ggsave(paste0("perc_infections_by_hh_cat_",fit_tag,".png"),width = 10,height=6)
write.csv(N_infected_by_hh_cat,paste0("perc_infections_by_hh_cat_",fit_tag,".csv") ,row.names = FALSE)


### Gropuping by age_cat
normal_isolation_inf_df$age_cat<-cut(normal_isolation_inf_df$age,breaks=c(-1,17,34,64,100),labels=c("Age: 0-17","Age: 18-34","Age: 35-64","Age: 65+"))
preventive_isolation_df$age_cat<-cut(preventive_isolation_df$age,breaks=c(-1,17,34,64,100),labels=c("Age: 0-17","Age: 18-34","Age: 35-64","Age: 65+"))
no_isolation_df$age_cat<-cut(no_isolation_df$age,breaks=c(-1,17,34,64,100),labels=c("Age: 0-17","Age: 18-34","Age: 35-64","Age: 65+"))


N_infected_norm_iso_by_agecat<-normal_isolation_inf_df %>% group_by(n,age_cat) %>% 
  summarise(Sum_infected = sum(ever_infected),Perc_infected=sum(ever_infected)/length(ever_infected),Potential_infected=length(ever_infected))
N_infected_norm_iso_by_agecat$net<-"Isolation after test"

N_infected_pre_iso_by_agecat<-preventive_isolation_df %>% group_by(n,age_cat) %>% 
  summarise(Sum_infected = sum(ever_infected),Perc_infected=sum(ever_infected)/length(ever_infected),Potential_infected=length(ever_infected))
N_infected_pre_iso_by_agecat$net<-"Isolation after symptoms"

N_infected_no_isolation_by_agecat<-no_isolation_df %>% group_by(n,age_cat) %>% 
  summarise(Sum_infected = sum(ever_infected),Perc_infected=sum(ever_infected)/length(ever_infected),Potential_infected=length(ever_infected))
N_infected_no_isolation_by_agecat$net<-"No isolation"

N_infected_by_agecat<-rbind(N_infected_no_isolation_by_agecat,N_infected_norm_iso_by_agecat,N_infected_pre_iso_by_agecat)

N_infected_by_agecat$net<-factor(x=N_infected_by_agecat$net,levels=c("No isolation","Isolation after test","Isolation after symptoms"))


### Adding data from FluTes
df_ever_infected$age_cat<-cut(df_ever_infected$age,breaks=c(-1,17,34,64,100),labels=c("Age: 0-17","Age: 18-34","Age: 35-64","Age: 65+"))
data_N_infected_by_agecat<-df_ever_infected %>% filter(index==0) %>% group_by(age_cat) %>% 
  summarise(Sum_infected = sum(ever_infected),Perc_infected=sum(ever_infected)/length(ever_infected),Potential_infected=length(ever_infected))

ggplot(data=N_infected_by_agecat)+geom_boxplot(outlier.size = -1,mapping = aes(y=Perc_infected,x=net))+ylab("Secondary infections")+xlab("Scenario")+
  facet_wrap(~age_cat,scales = "free_y")+theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
ggsave(paste0("perc_infections_by_age_",fit_tag,".png"),width = 10,height=6)

write.csv(N_infected_by_agecat,paste0("perc_infections_by_age_",fit_tag,".csv") ,row.names = FALSE)





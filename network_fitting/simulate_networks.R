# Rcode household contact networks
# Last update: v20211123
# Last updat by Pietro on 18/12/2023
rm(list = ls())
########################################################
`%notin%` <- function(a,b) ! a %in% b  
### Automatically set working directory
if(require(rstudioapi) && isAvailable()){
  current_path <- getActiveDocumentContext()$path 
  setwd(dirname(current_path ))
}

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
source("functions_recursive.R") ## Needed for simulation to define "AC_levels"
N_sims<-10       ## How many networks to simulate
Net_per_batch<-10 ## How many to save together
N_batches<-ceiling(N_sims/Net_per_batch)
dum<-read.csv("../epi_fitting/epidemic_info/epidemic_info_net_onset_yest.csv")
memids_in_common_yest_onset<-dum$memid
##########################################################
### Simulate N_sims onset networks and selecting only hhids in common with yest
##########################################################
model_onset_full<-readRDS("model_results/best_models/model_full_data_net_onset_edges_quadratic_2stars_quadratic_triangles_linear.RDS")
set.seed(24112022) ## Set the seed
this_sim_network_onset<-simulate(model_onset_full,nsim=1)
vids_to_keep<-which(get.vertex.attribute(this_sim_network_onset,"vertex.names") %in%  memids_in_common_yest_onset)

for(batch_counter in 1:N_batches){
  sim_network_onset<-list()
  print(paste0("batch ",batch_counter," of ",N_batches ))
  for(counter in 1:Net_per_batch){
    print(paste0("    Sim in batch:",counter," of ",Net_per_batch ))  
  this_sim_network_onset<-simulate(model_onset_full,nsim=1)
  ## Select only memids that are in common between onset and yest
  this_sim_network_onset<-get.inducedSubgraph(this_sim_network_onset, vids_to_keep)
  sim_network_onset[[counter]]<-this_sim_network_onset
  }
  saveRDS(sim_network_onset,paste0("simulated_networks/onset_batch_", batch_counter,".RDS"))
}

##########################################################
### Simulate N_sims yest networks and selecting only hhids in common with onset
##########################################################
model_yest_full<-readRDS("model_results/model_full_data_net_yest_edges_quadratic_2stars_quadratic_triangles_linearwith_interaction.RDS")
set.seed(24112022) ## Re-set the seed
this_sim_network_yest<-simulate(model_yest_full,nsim=1)
vids_to_keep<-which(get.vertex.attribute(this_sim_network_yest,"vertex.names") %in%  memids_in_common_yest_onset)
for(batch_counter in 1:N_batches){
  sim_network_yest<-list()
  print(paste0("batch ",batch_counter," of ",N_batches ))
  for(counter in 1:Net_per_batch){
    print(paste0("    Sim in batch:",counter," of ",Net_per_batch ))  
    this_sim_network_yest<-simulate(model_yest_full,nsim=1)
    ## Select only memids that are in common between onset and yest
    this_sim_network_yest<-get.inducedSubgraph(this_sim_network_yest, vids_to_keep)
    sim_network_yest[[counter]]<-this_sim_network_yest
  }
  saveRDS(sim_network_yest,paste0("simulated_networks/yest_batch_", batch_counter,".RDS"))
}

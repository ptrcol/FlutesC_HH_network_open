# Rcode household contact networks
## This file generates the HH network object starting from data
# Last update: 02/06/2023

###################################
rm(list = ls())
########################################################
`%notin%` <- function(a,b) ! a %in% b  
### Automatically set working directory
if(require(rstudioapi) && isAvailable()){
  current_path <- getActiveDocumentContext()$path 
  setwd(dirname(current_path ))
}
## Load libraries
library(ggplot2)
library(ergm)
library(ergm.multi)
library(statnet.common)
library(dplyr)


## Functions
`%notin%` <- function(a,b) ! a %in% b  

source("general_functions.R")
source("plot_functions.R")
load("../full_data_set/FluTESC_shortversion_26_01_2023.RData")

## Compute the time between enrollment and onset of symptoms according to CDC
out$delay_isolation<-difftime(out$doe_dt,out$cdc_primary_onsetdt,units = "days") #Yuwei suggestion from call 20/01/2023
## Removing individuals with negative value
hhid_to_remove<-out$hhid[which(out$delay_isolation<0)] ## Selecting the corresponding HH
out<-out[-which(out$hhid %in% hhid_to_remove),]        ## Dropping the hh


out<-out %>%
  group_by(hhid) %>%
  add_count(name = "hhmem_num")

## Load  old data (these is needed to have the name of the columns of the
## "separated" files, which will then be used to generate from the merged data 
## file the same structure as before
baselinedd_old<-read.csv("../Data_Belgium_202201/baselinediary.csv")
baselineiiq_old<-read.csv("../Data_Belgium_202201/baselineiiq.csv")
dailydiary_old<-read.csv("../Data_Belgium_202201/dailydiary.csv")
iiq_old<-read.csv("../Data_Belgium_202201/iiq.csv")
lab_combined_old<-read.csv("../Data_Belgium_202201/lab_combined.csv")
subject_data_old<-read.csv("../Data_Belgium_202201/subject_data_nophi.csv")

## Storing the colnames for each individual file
baselinedd_colnames<-colnames(baselinedd_old)
baselineiiq_colnames<-colnames(baselineiiq_old)
dailydiary_colnames<-colnames(dailydiary_old)
iiq_colnames<-colnames(iiq_old)
lab_combined_data_colnames<-colnames(lab_combined_old)
subject_data_colnames<-colnames(subject_data_old)
### Adding site colname
subject_data_colnames<-c(subject_data_colnames,"site")

## "divide" the merged data file into the separated ones (to "recycle" the old
## way of generating the network objects
baselinedd<-out[,baselinedd_colnames[baselinedd_colnames %in% colnames(out)]]
baselineiiq<-out[,baselineiiq_colnames[baselineiiq_colnames %in% colnames(out)]]
dailydiary<-out[,dailydiary_colnames[dailydiary_colnames %in% colnames(out)]]
iiq<-out[,iiq_colnames[iiq_colnames %in% colnames(out)]]
lab_combined<-out[,lab_combined_data_colnames[lab_combined_data_colnames %in% colnames(out)]]
subject_data<-out[,subject_data_colnames[subject_data_colnames %in% colnames(out)]]
diary<-dailydiary
subject<-subject_data

## match subjects from diary and subject data
# keep only the subjects that were selected (and filled in diaries)
subject<-subject[subject$memid != ".",]
subject<-subject[subject$memid %in% diary$memid,]



####################
lab_combined<-lab_combined[lab_combined$memid %in% diary$memid,]
cases<-lab_combined %>% group_by(memid) %>% dplyr::summarize(daysPos=sum(sarscov2_pcr_result==1),n=sum(sarscov2_pcr_result<=1),case=max(sarscov2_pcr_result==1))
## Add baseline information
cases<-merge(cases,subject,by="memid")
cases$index<-substr(cases$memid, 7, 7)=="1"

## SUBSETS OF DATA
# all households
cases_total<-cases
tmp<-data.frame(table(cases$hhid,cases$case))
hhid_2cases<-tmp$Var1[tmp$Var2==1 & tmp$Freq>1]
df_infected_hh<-data.frame("hh_id"=cases$hhid)
df_infected_hh$secondary_case<-0
df_infected_hh$secondary_case[df_infected_hh$hh_id %in%hhid_2cases]<-1
write.csv(df_infected_hh,"infected_households_full_data.csv")
write.csv(subject,"subject_data_for_analysis_of_infection.csv") ## Done to save some info for other analysis

k=2 ## Definition for debugging
subject$noNA_hh_size<-NA
output<-NULL
for (k in 1:length(unique(subject$hhid))){ 
  print(k)
  hhid.sel=(subject$hhid==unique(subject$hhid)[k])
  hhsize=cases_total$hhmem_num[cases_total$hhid==unique(subject$hhid)[k]]
  if (length(hhsize)==0){hhsize=NA}
  subject$noNA_hh_size[hhid.sel]<-hhsize
  l<-1 ## Debug
  for (l in 1:length(subject$memid[hhid.sel])){
    memid.sel=subject$memid==subject$memid[hhid.sel][l]
    ## potential contacts
    memid.potential.contacts=c(subject$cnt1_id[memid.sel],subject$cnt2_id[memid.sel],subject$cnt3_id[memid.sel],subject$cnt4_id[memid.sel],subject$cnt5_id[memid.sel],subject$cnt6_id[memid.sel],subject$cnt7_id[memid.sel],subject$cnt8_id[memid.sel],subject$cnt9_id[memid.sel],subject$cnt10_id[memid.sel])
    if(all(is.na(memid.potential.contacts))==FALSE){
      ## add intimacy levels here in a vector by contact and cbind in the next command
      memid.contacts.yest=c(subject$cnt1_yest[memid.sel],subject$cnt2_yest[memid.sel],subject$cnt3_yest[memid.sel],subject$cnt4_yest[memid.sel],subject$cnt5_yest[memid.sel],subject$cnt6_yest[memid.sel],subject$cnt7_yest[memid.sel],subject$cnt8_yest[memid.sel],subject$cnt9_yest[memid.sel],subject$cnt10_yest[memid.sel])
      memid.contacts.onset=c(subject$cnt1_onset[memid.sel],subject$cnt2_onset[memid.sel],subject$cnt3_onset[memid.sel],subject$cnt4_onset[memid.sel],subject$cnt5_onset[memid.sel],subject$cnt6_onset[memid.sel],subject$cnt7_onset[memid.sel],subject$cnt8_onset[memid.sel],subject$cnt9_onset[memid.sel],subject$cnt10_onset[memid.sel])
      memid.contacts.final=c(subject$cnt1_final[memid.sel],subject$cnt2_final[memid.sel],subject$cnt3_final[memid.sel],subject$cnt4_final[memid.sel],subject$cnt5_final[memid.sel],subject$cnt6_final[memid.sel],subject$cnt7_final[memid.sel],subject$cnt8_final[memid.sel],subject$cnt9_final[memid.sel],subject$cnt10_final[memid.sel])
      memid.contacts.yestphys=c(subject$cnt1_yestphys[memid.sel],subject$cnt2_yestphys[memid.sel],subject$cnt3_yestphys[memid.sel],subject$cnt4_yestphys[memid.sel],subject$cnt5_yestphys[memid.sel],subject$cnt6_yestphys[memid.sel],subject$cnt7_yestphys[memid.sel],subject$cnt8_yestphys[memid.sel],subject$cnt9_yestphys[memid.sel],subject$cnt10_yestphys[memid.sel])
      memid.contacts.onsetphys=c(subject$cnt1_onsetphys[memid.sel],subject$cnt2_onsetphys[memid.sel],subject$cnt3_onsetphys[memid.sel],subject$cnt4_onsetphys[memid.sel],subject$cnt5_onsetphys[memid.sel],subject$cnt6_onsetphys[memid.sel],subject$cnt7_onsetphys[memid.sel],subject$cnt8_onsetphys[memid.sel],subject$cnt9_onsetphys[memid.sel],subject$cnt10_onsetphys[memid.sel])
      memid.contacts.finalphys=c(subject$cnt1_finalphys[memid.sel],subject$cnt2_finalphys[memid.sel],subject$cnt3_finalphys[memid.sel],subject$cnt4_finalphys[memid.sel],subject$cnt5_finalphys[memid.sel],subject$cnt6_finalphys[memid.sel],subject$cnt7_finalphys[memid.sel],subject$cnt8_finalphys[memid.sel],subject$cnt9_finalphys[memid.sel],subject$cnt10_finalphys[memid.sel])
      memid.contacts.onsetbed=c(subject$cnt1_onsetbed[memid.sel],subject$cnt2_onsetbed[memid.sel],subject$cnt3_onsetbed[memid.sel],subject$cnt4_onsetbed[memid.sel],subject$cnt5_onsetbed[memid.sel],subject$cnt6_onsetbed[memid.sel],subject$cnt7_onsetbed[memid.sel],subject$cnt8_onsetbed[memid.sel],subject$cnt9_onsetbed[memid.sel],subject$cnt10_onsetbed[memid.sel])
      memid.contacts.yestbed=c(subject$cnt1_yestbed[memid.sel],subject$cnt2_yestbed[memid.sel],subject$cnt3_yestbed[memid.sel],subject$cnt4_yestbed[memid.sel],subject$cnt5_yestbed[memid.sel],subject$cnt6_yestbed[memid.sel],subject$cnt7_yestbed[memid.sel],subject$cnt8_yestbed[memid.sel],subject$cnt9_yestbed[memid.sel],subject$cnt10_yestbed[memid.sel])
      memid.contacts.finalbed=c(subject$cnt1_finalbed[memid.sel],subject$cnt2_finalbed[memid.sel],subject$cnt3_finalbed[memid.sel],subject$cnt4_finalbed[memid.sel],subject$cnt5_finalbed[memid.sel],subject$cnt6_finalbed[memid.sel],subject$cnt7_finalbed[memid.sel],subject$cnt8_finalbed[memid.sel],subject$cnt9_finalbed[memid.sel],subject$cnt10_finalbed[memid.sel])
      if(all(sapply(c(subject$hhid[memid.sel],hhsize,subject$memid[hhid.sel][l],
                  memid.potential.contacts[!is.na(memid.potential.contacts)],
                  memid.contacts.onset[!is.na(memid.potential.contacts)],
                  memid.contacts.yest[!is.na(memid.potential.contacts)],
                  memid.contacts.final[!is.na(memid.potential.contacts)],
                  memid.contacts.onsetphys[!is.na(memid.potential.contacts)],
                  memid.contacts.yestphys[!is.na(memid.potential.contacts)],
                  memid.contacts.finalphys[!is.na(memid.potential.contacts)],
                  memid.contacts.onsetbed[!is.na(memid.potential.contacts)],
                  memid.contacts.yestbed[!is.na(memid.potential.contacts)],
                  memid.contacts.finalbed[!is.na(memid.potential.contacts)]
      ),length)==length(c(subject$hhid[memid.sel],hhsize,subject$memid[hhid.sel][l],
                          memid.potential.contacts[!is.na(memid.potential.contacts)],
                          memid.contacts.onset[!is.na(memid.potential.contacts)],
                          memid.contacts.yest[!is.na(memid.potential.contacts)],
                          memid.contacts.final[!is.na(memid.potential.contacts)],
                          memid.contacts.onsetphys[!is.na(memid.potential.contacts)],
                          memid.contacts.yestphys[!is.na(memid.potential.contacts)],
                          memid.contacts.finalphys[!is.na(memid.potential.contacts)],
                          memid.contacts.onsetbed[!is.na(memid.potential.contacts)],
                          memid.contacts.yestbed[!is.na(memid.potential.contacts)],
                          memid.contacts.finalbed[!is.na(memid.potential.contacts)]
      )[[1]]))==FALSE){print(paste0("Problem"))}
      
      output=rbind(output,cbind(subject$hhid[memid.sel],unique(hhsize),subject$memid[hhid.sel][l],
                          memid.potential.contacts[!is.na(memid.potential.contacts)],
                          memid.contacts.onset[!is.na(memid.potential.contacts)],
                          memid.contacts.yest[!is.na(memid.potential.contacts)],
                          memid.contacts.final[!is.na(memid.potential.contacts)],
                          memid.contacts.onsetphys[!is.na(memid.potential.contacts)],
                          memid.contacts.yestphys[!is.na(memid.potential.contacts)],
                          memid.contacts.finalphys[!is.na(memid.potential.contacts)],
                          memid.contacts.onsetbed[!is.na(memid.potential.contacts)],
                          memid.contacts.yestbed[!is.na(memid.potential.contacts)],
                          memid.contacts.finalbed[!is.na(memid.potential.contacts)]
                        )
      )
    }
  } #for (l in 1:length(subject$memid[hhid.sel])){
}
edges=data.frame(output)
names(edges)=c("hhid","hhsize","part","cont","onset","yest","final","onsetphys","yestphys","finalphys","onsetbed","yestbed","finalbed")#,"phys")

### Removing contacts of participants with themselves
edges<-edges[edges$part!=edges$cont,]



### Info on non enrolled (not relevant for analysis in the paper)
non_enrolled_data<-read.csv("../Data_Belgium_202201/info_non_enrolled.csv")
colnames(non_enrolled_data)<-c("memid","hce_age","hce_sex")
### Removing non-ernolled participants for which we don't know the age
non_enrolled<-subset(non_enrolled_data,!is.na(non_enrolled_data$hce_age))
## Setting <1 values to 0
non_enrolled$hce_age[non_enrolled$hce_age<1]<-0




## Storing edges object
## Sharing a bedroom ("XXXbed") is considered as a physical contact.
## Networks without bed sharing are tagged as "no_bed" (so, they include only physical contacts).
## The networks of only bed sharing are tagged as "XXXbed".
 

edges.full=edges
edges.finalphys=subset(edges,(edges$final!="4"&edges$finalphys==1)  | edges$finalbed==1 )  
edges.yestphys=subset(edges,(edges$yest!="4"&edges$yestphys==1)  | edges$yestbed==1 )  
edges.onsetphys=subset(edges,(edges$onset!="4"& edges$onsetphys==1) | edges$onsetbed==1 )  



write.csv(edges.full, "edges_FULL.csv", row.names = FALSE)
write.csv(edges.finalphys, "edges_final_phys.csv", row.names = FALSE)
write.csv(edges.yestphys, "edges_yest_phys.csv", row.names = FALSE)
write.csv(edges.onsetphys, "edges_onset_phys.csv", row.names = FALSE)


edges.onsetbed=subset(edges,edges$onsetbed==1)
edges.yestbed=subset(edges,edges$yestbed==1)
edges.finalbed=subset(edges,edges$finalbed==1)

edges.onsetphys_nobed =subset(edges,edges$onset!="4"&edges$onsetphys==1)  
edges.yestphys_nobed  =subset(edges,edges$onset!="4"&edges$yestphys==1)  
edges.finalphys_nobed =subset(edges,edges$onset!="4"&edges$finalphys==1)  





source("func_generate_net_object.R")

################################################################################
#############    ONSET NETWORK    ##########
################################################################################
generate_net_obj_full_data(edges.onsetphys,subject,diary,
                           "full_dataset_onset_phys",tag_plot=FALSE)

################################################################################
#############    YEST NETWORK    ##########
################################################################################
generate_net_obj_full_data(edges.yestphys,subject,diary,
                           "full_dataset_yest_phys",tag_plot=FALSE)

################################################################################
#############    Final NETWORK    ##########
################################################################################
generate_net_obj_full_data(edges.finalphys,subject,diary,
                           "full_dataset_final_phys",tag_plot=FALSE)

################################################################################
#############    ONSET/YEST NETWORKs in common    ##########
### Note that from here onward there is the need of the "epidemic_info" files
### generated by the R file "epi_fitting/generate_epidemic_info.R"
################################################################################
dum<-read.csv(file='../epi_fitting/epidemic_info/epidemic_info_net_onset_yest.csv')  
memids_in_common_yest_onset<-dum$memid
net_onset<-readRDS("./net_objects/graph_full_dataset_onset_phys_onlypart_symmetric.RDS")
vids_to_keep<-which(get.vertex.attribute(net_onset,"vertex.names") %in%  memids_in_common_yest_onset)
net_onset_common_yest<-get.inducedSubgraph(net_onset, vids_to_keep)
saveRDS(object = net_onset_common_yest,"net_objects/graph_full_dataset_onset_phys_in_common_yest_onlypart_symmetric.RDS")
net_onset_common_yest_by_hh<-ergm.multi::uncombine_network(net_onset_common_yest,split.vattr = "hhid",names.vattr="hhid")
saveRDS(object = net_onset_common_yest_by_hh,"net_objects/graph_full_dataset_onset_phys_in_common_yest_onlypart_symmetric_splitted_by_HH.RDS")
## New way of doing this, to avoid error in fitting ("network_modelling_full_data_recursive.R")
dum<-read.csv(file='../epi_fitting/epidemic_info/epidemic_info_net_onset_yest.csv')  
hhids_in_common_yest_onset<-dum$hhid
edges.common_yest_onset=subset(edges.full,edges.full$hhid %in% hhids_in_common_yest_onset)
generate_net_obj_full_data(edges.common_yest_onset,subject %>% filter(hhid %in% hhids_in_common_yest_onset),
                           diary %>% filter(hhid %in% hhids_in_common_yest_onset),
                           "full_dataset_onset_phys_in_common_yest_bis",tag_plot=FALSE)



net_yest<-readRDS("./net_objects/graph_full_dataset_yest_phys_onlypart_symmetric.RDS")
vids_to_keep<-which(get.vertex.attribute(net_yest,"vertex.names") %in%  memids_in_common_yest_onset)
net_yest_common_onset<-get.inducedSubgraph(net_yest, vids_to_keep)
saveRDS(object = net_yest_common_onset,"net_objects/graph_full_dataset_yest_phys_in_common_onset_onlypart_symmetric.RDS")

net_yest_common_onset_by_hh<-ergm.multi::uncombine_network(net_yest_common_onset,split.vattr = "hhid",names.vattr="hhid")
saveRDS(object = net_yest_common_onset_by_hh,"net_objects/graph_full_dataset_yest_phys_in_common_onset_onlypart_symmetric_splitted_by_HH.RDS")


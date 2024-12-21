# Rcode to generate several epidemic information (so e.g. infection times, etc)
# Last updat by Pietro on 02/06/2023
# library(igraph)
library(ergm)
library(ergm.multi)
library(dplyr)
library(lubridate)
library(latex2exp)
library(ggplot2)
library(gridExtra)
rm(list = ls())
########################################################
`%notin%` <- function(a,b) ! a %in% b  
### Automatically set working directory
if(require(rstudioapi) && isAvailable()){
  current_path <- getActiveDocumentContext()$path 
  setwd(dirname(current_path ))
}

load("../full_data_set/FluTESC_shortversion_26_01_2023.RData")

#### Generate a list of memid from the networks (and the corresponding list of tags)
net_onset<-readRDS("../network_obj_generation/net_objects/graph_full_dataset_onset_phys_onlypart_symmetric.RDS")
net_yest<-readRDS("../network_obj_generation/net_objects/graph_full_dataset_yest_phys_onlypart_symmetric.RDS")
net_final<-readRDS("../network_obj_generation/net_objects/graph_full_dataset_final_phys_onlypart_symmetric.RDS")
memid_onset<-get.vertex.attribute(net_onset,"vertex.names")
memid_yest<-get.vertex.attribute(net_yest,"vertex.names")
memid_final<-get.vertex.attribute(net_final,"vertex.names")
memid_onset_yest<-intersect(memid_onset,memid_yest)
memid_onset_yest_final<-intersect(memid_onset_yest,memid_final)

list_memids<-list(memid_onset,memid_yest,memid_final,memid_onset_yest,memid_onset_yest_final)
list_tags<-c("net_onset","net_yest","net_final","net_onset_yest","net_onset_yest_final")




# Define time of infection of secondary cases
out$dt_infection<-difftime(out$firststudypos_dt,out$cdc_primary_onsetdt,units = "days") #Email suggestion Alexandra (mail1.png)
# Define duration of positive test
out$duration_pos<-difftime(out$lastpos_dt,out$firststudypos_dt,units = "days")
out$duration_pos[out$duration_pos>14]<-0


# Define time from symtpoms onset of index case to enrollment
out$delay_isolation<-difftime(out$doe_dt,out$cdc_primary_onsetdt,units = "days") #Yuwei suggestion from call 20/01/2023
out$delay_isolation[out$delay_isolation<0]<-0
# Define time from symptoms onset to last day of test (max dt of the simulation)
out$max_t_simulation<-difftime(out$lastswab_dt,out$cdc_primary_onsetdt,units = "days")



# defining incidence vector
non_index_out<-out %>%filter(index==0)
incidence<-table(non_index_out$dt_infection)
incidence<-as.integer(incidence)
write.csv(x=incidence,file='incidence_full.csv',row.names = FALSE)


epi_info<-out %>% group_by(memid)%>% 
  dplyr::summarize(case=max(sarscov2_pcr_result==1),
                   index=index,
                   hhid=hhid,
                   cdc_primary,
                   delay_isolation=delay_isolation,
                   duration_pos=duration_pos,
                   dt_infection=dt_infection,
                   max_t_simulation=max_t_simulation,
                   site=site,
                   cdc_primary_onsetdt=cdc_primary_onsetdt
                   )
write.csv(x=epi_info,file='./epidemic_info/epidemic_info_full.csv',row.names = FALSE)


for(counter in 1:length(list_memids)){
  memids_to_save<-unlist(list_memids[counter])
  epi_info_net_specific<-epi_info %>% filter(memid%in% memids_to_save)
  epi_info_net_specific<-epi_info_net_specific[order(match(epi_info_net_specific$memid,memids_to_save)),]
  fname<-paste0("./epidemic_info/epidemic_info_",list_tags[counter], ".csv")
  write.csv(x=epi_info_net_specific,file=fname,row.names = FALSE)
  ### Restricting to non-index cases
  non_index_out_net_specific<-epi_info_net_specific %>%filter(index==0)
  incidence_net_specific<-table(non_index_out_net_specific$dt_infection)
  incidence_net_specific<-as.integer(incidence_net_specific)
  fname<-paste0("incidence_",list_tags[counter], ".csv")
  write.csv(x=incidence,file=fname,row.names = FALSE)
}

### Generate the matrix of infection
# (this is a matrix with rows=participant, columns=time  and value is infected/not infected)
Ncols=max(as.integer(out$duration_pos)+as.integer(out$dt_infection),na.rm = TRUE)
time_inf_matrix<-matrix(0,nrow=length(out$memid),ncol = Ncols)
rownames(time_inf_matrix)<-out$memid
for(counter in (1:length(out$memid))){
  this_memid<-out$memid[counter]
  this_inf_dt<-as.integer(out$dt_infection[counter])
  this_inf_duration<-as.integer(out$duration_pos[counter])
  if(!is.na(this_inf_dt)){
    time_inf_matrix[counter,this_inf_dt:(this_inf_dt+this_inf_duration)]<-1
    }
}
write.csv(time_inf_matrix,"time_infection_matrix_full.csv")
### Write network specific matrix of infection
for(counter_case in 1:length(list_memids)){
  print(counter_case)
  memids_to_save<-unlist(list_memids[counter_case])
  net_specific_out<-out %>% filter(memid%in% memids_to_save)
  net_specific_out<-net_specific_out[order(match(net_specific_out$memid,memids_to_save)),]
  Ncols=max(as.integer(net_specific_out$duration_pos)+as.integer(net_specific_out$dt_infection),na.rm = TRUE)
  time_inf_matrix_net_specific<-matrix(0,nrow=length(net_specific_out$memid),ncol = Ncols)
  rownames(time_inf_matrix_net_specific)<-net_specific_out$memid
  for(counter in (1:length(net_specific_out$memid))){
    this_memid<-net_specific_out$memid[counter]
    this_inf_dt<-as.integer(net_specific_out$dt_infection[counter])
    this_inf_duration<-as.integer(net_specific_out$duration_pos[counter])
    if(!is.na(this_inf_dt)){
      time_inf_matrix_net_specific[counter,this_inf_dt:(this_inf_dt+this_inf_duration)]<-1
    }
  }
  fname<-paste0("time_infection_matrix_",list_tags[counter_case], ".csv")
  write.csv(time_inf_matrix_net_specific,fname)
}





################################################################################
## Generate community incidence  #
################################################################################
## Define minimum and maximum date to get from the commmunity incidence
min_date<-min(min(out$doe_dt),min(out$cdc_primary_onsetdt),min(out$firstpos_dt,na.rm = TRUE))-1
max_date<-max(max(out$lastswab_dt,na.rm = TRUE))+50
list_of_dates<-seq(as.Date(min_date), as.Date(max_date), by = "day")

clean_weekly_cases<-function(df_data,int_list_of_dates){
  df_data$Date<-as.Date(df_data$Date)
  df_data<- df_data %>% filter(Date %in% int_list_of_dates)
  df_data<-df_data %>% filter(Weekly.Cases>0)
  return(df_data)
}


generate_comm_transmission<-function(df_data){ ## Old: for county data
  ## Selecting interesting columns
  df_commm_transm<-data.frame(date=as.Date(df_data$Date),perc_pos_test=as.numeric(df_data$X7.day.Percent.Positivity),weekly_cases=df_data$Weekly.Cases)
  ## Sorting by date
  df_commm_transm[order(df_commm_transm$date),]
  ## Including "ratios"
  df_commm_transm$ratio_weekly_cases<-df_commm_transm$weekly_cases/df_commm_transm$weekly_cases[1]
  df_commm_transm$ratio_perc_pos_test<-df_commm_transm$perc_pos_test/df_commm_transm$perc_pos_test[1]
  return(df_commm_transm)

}



generate_comm_transmission_state<-function(df_data){
  ## Selecting interesting columns
  df_commm_transm<-data.frame(date=as.Date(x = df_data$Date,format="%b %d %Y"),weekly_cases=df_data$Weekly.Cases)
  ## Sorting by date
  df_commm_transm[order(df_commm_transm$date),]
  ## Including "ratios"
  df_commm_transm$ratio_weekly_cases<-df_commm_transm$weekly_cases/df_commm_transm$weekly_cases[1]
  return(df_commm_transm)
  
}



dum_state_tennessee<-read.csv("data_table_for_weekly_case_trends_tennessee.csv",skip = 2)  
dum_state_wisconsin<-read.csv("data_table_for_weekly_case_trends_wisconsin.csv",skip = 2)  
df_transm_tennessee_state<-generate_comm_transmission_state(dum_state_tennessee)
df_transm_wisconsin_state<-generate_comm_transmission_state(dum_state_wisconsin)


df_transm_tennessee_state$location<-"VUMC"
df_transm_wisconsin_state$location<-"MCRI"
df_transm_state<-rbind(df_transm_tennessee_state,df_transm_wisconsin_state)
write.csv(df_transm_state,"community_transmission_state.csv",row.names = FALSE)

### add "a)","b)" etc
plot_a<-ggplot(df_transm_tennessee_state)+geom_point(mapping=aes(x=date,y=ratio_weekly_cases))+ylim(0,25)+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+scale_x_date(limits =c(as.Date("2020-04-01"),as.Date("2021-05-01")) ,date_breaks = "months" , date_labels = "%b-%y")+
  xlab("")+ylab((TeX("$\\frac{{Incidence[date]}}{Incidence[April-2020]}$")))#+ggtitle("Tennessee") #  i.e. (Vanderbilt University Medical Center)

plot_b<-ggplot(out %>% filter(site=="VUMC"))+geom_histogram(mapping=aes(x=(ceiling_date(doe_dt, "week")+3)))+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+scale_x_date(limits =c(as.Date("2020-04-01"),as.Date("2021-05-01")) ,date_breaks = "months" , date_labels = "%b-%y")+
xlab("Date")+ylab("Number \n of participants ")

#plot_nashville<-grid.arrange(plot_a, plot_b, nrow = 2,ncol=1)



plot_c<-ggplot(df_transm_wisconsin_state)+geom_point(mapping=aes(x=date,y=ratio_weekly_cases))+ylim(0,25)+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+scale_x_date(limits =c(as.Date("2020-04-01"),as.Date("2021-05-01")) ,date_breaks = "months" , date_labels = "%b-%y")+
  xlab("")+ylab((TeX("$\\frac{{Incidence[date]}}{Incidence[April-2020]}$")))#+ggtitle("Wisconsin") # i.e. Marshfield Clinic Research Institute

plot_d<-ggplot(out %>% filter(site=="MCRI"))+geom_histogram(mapping=aes(x=(ceiling_date(doe_dt, "week")+3)))+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+scale_x_date(limits =c(as.Date("2020-04-01"),as.Date("2021-05-01")) ,date_breaks = "months" , date_labels = "%b-%y")+
  xlab("Date")+ylab("Number \n of participants ")

#plot_marshville<-grid.arrange(plot_c, plot_d, nrow = 2,ncol=1)


png("participation_vs_state_incidence.png", width = 1200, height = 800) # Open a new pdf file
grid.arrange(plot_a,plot_c,plot_b, plot_d, nrow = 2,ncol=2)
dev.off() # Close the file

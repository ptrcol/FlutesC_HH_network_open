# Rcode household contact networks
# Last update: v20211123
# Last update by Pietro on 23/03/2022
rm(list = ls())
########################################################
`%notin%` <- function(a,b) ! a %in% b  
### Automatically set working directory
if(require(rstudioapi) && isAvailable()){
  current_path <- getActiveDocumentContext()$path 
  setwd(dirname(current_path ))
}
library(ggplot2)

fit_tag<-"main_analysis"  
N_infected_by_hh_cat<-read.csv(paste0("./data/perc_infections_by_hh_cat_",fit_tag,".csv"))
N_infected_by_hh_cat$net<-factor(N_infected_by_hh_cat$net,levels = c("Isolation after symptoms","Isolation after test","No isolation"),
                                 labels=c("Physical distancing \n after symptoms","Physical distancing \n after test","No physical distancing"))
N_infected_by_hh_cat$hhmem_num_cat<-factor(N_infected_by_hh_cat$hhmem_num_cat, levels=c("HH size 2","HH size 3","HH size 4","HH size 5+"),
                                           labels=c("Household size 2","Household size 3","Household size 4","Household size 5 and more"))



ggplot(data=N_infected_by_hh_cat)+geom_boxplot(outlier.size = -1,mapping = aes(y=Perc_infected,x=net))+ylab("Probability of \n secondary infection")+xlab("Scenario")+
  facet_wrap(~hhmem_num_cat)+theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
ggsave(paste0("Figure_sec_cases_by_hh.png"),width = 10,height=6)

ggplot(data=N_infected_by_hh_cat)+geom_violin(outlier.size = -1,mapping = aes(y=Perc_infected,x=net))+ylab("Probability of \n secondary infection")+xlab("Scenario")+
  facet_wrap(~hhmem_num_cat)+theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))

library(dplyr)
table<- N_infected_by_hh_cat %>% group_by(hhmem_num_cat,net) %>% summarise(
  mean=mean(Perc_infected),
  perc_2.5 = quantile(Perc_infected, probs = 0.025),
  perc_25 = quantile(Perc_infected, probs = 0.25),
  perc_50 = quantile(Perc_infected, probs = 0.50),
  perc_75 = quantile(Perc_infected, probs = 0.75),
  perc_97.5 = quantile(Perc_infected, probs = 0.975)
)

table_b<-table %>% mutate(
  ,ratio_mean = mean / mean[net ==  "No physical distancing"]
  ,ratio_p2_5 = perc_2.5 / mean[net ==  "No physical distancing"]
  ,ratio_p97_5 = perc_97.5 / mean[net ==  "No physical distancing"]
)
table_b<-table_b[,c(1,2,grep("ratio",colnames(table_b)))]


N_infected_by_agecat<-read.csv(paste0("./data/perc_infections_by_age_",fit_tag,".csv"))
N_infected_by_agecat$net<-factor(N_infected_by_agecat$net,levels = c("Isolation after symptoms","Isolation after test","No isolation"),
                                 labels=c("Physical distancing \n after symptoms","Physical distancing \n after test","No physical distancing"))


ggplot(data=N_infected_by_agecat)+geom_boxplot(outlier.size = -1,mapping = aes(y=net,x=Perc_infected))+xlab("Probability of \n secondary infection")+ylab("Scenario")+
  facet_wrap(~age_cat,scales = "free_y")+theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))#+  geom_point(data=data_N_infected_by_agecat,mapping=aes(x="Isolation after test",y=Perc_infected),color="red")
ggsave(paste0("Figure_sec_cases_by_age.png"),width = 10,height=6)






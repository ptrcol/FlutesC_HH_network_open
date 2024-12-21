# Rcode household contact networks
# Last update by Pietro on 26/11/2024
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
library(latex2exp)
library(modelsummary)
library(statnet)
library(ParBayesianOptimization)
library(doParallel)
library(tidyverse)
source("../sim_libs/sim_models.R")

# if it does not exist, create subfolder "BayOpt_onset" in current working directory
if(!dir.exists("BayOpt_onset")){dir.create("BayOpt_onset")}

HH_network_no_isolation<-readRDS("../network_obj_generation/net_objects/graph_full_dataset_onset_phys_in_common_yest_onlypart_symmetric.RDS")
HH_network_w_isolation<-readRDS("../network_obj_generation/net_objects/graph_full_dataset_yest_phys_in_common_onset_onlypart_symmetric.RDS")
score_tag<-"main_analysis"

DO_PRE_GRID<-FALSE
DO_100<-FALSE
DO_MORE<-FALSE
### Type parallel can be either "FORK" or "PSOCK". Keep in mind that "FORK" is not available on Windows but is faster.
TYPE_PARALLEL<-"FORK"
n.cores <- parallel::detectCores() - 2
n.cores <- 6  ## To set manually to reduce the parallelization

set.seed(24112022)
beta1_min<-5e-2
beta1_max<-4e-1
beta2_min<-1e-4
beta2_max<-1e-1

if(DO_PRE_GRID){
  Nsteps_grid<-30

  starting_grid<-expand.grid(beta1=seq(from=beta1_min,to=beta1_max,by=(beta1_max-beta1_min)/Nsteps_grid),
                             beta2=seq(from=beta2_min,to=beta2_max,by=(beta2_max-beta2_min)/Nsteps_grid))
  
  list_scores<-list()
  library(foreach)
  n.cores <- parallel::detectCores() - 2
  my.cluster <- parallel::makeCluster(n.cores,type = "PSOCK")
  doParallel::registerDoParallel(cl = my.cluster)
  foreach::getDoParRegistered()
  list_scores <- foreach(i = 1:dim(starting_grid)[1],
                         .packages = "ergm.multi"
                         ) %dopar% {
                           source("../sim_libs/sim_models.R",local=TRUE)
    SEIR_simulate_two_networks_forBOptim_Binloglik_Ncases_by_HH_averageN50(starting_grid$beta1[i],
                                                                   starting_grid$beta2[i])
  }
  parallel::stopCluster(cl = my.cluster)
  
  starting_grid$scores<-unlist(list_scores)
  ggplot(starting_grid)+geom_raster(aes(x=beta1,y=beta2,fill=scores))
  saveRDS(list_scores,file=paste0("Scores_grid_",score_tag, ".RDS"))
  saveRDS(starting_grid,file=paste0("grid_",score_tag,".RDS"))
  # 
  list_scores<-readRDS(file=paste0("Scores_grid_",score_tag, ".RDS"))
  starting_grid<-readRDS(file=paste0("grid_",score_tag,".RDS"))
  starting_grid$scores<-unlist(list_scores)
  starting_grid<-starting_grid[order(starting_grid$scores,decreasing = TRUE),]
  
  #remove from starting_grid all lines with beta1=beta1_min and beta2=beta2_min
  starting_grid<-starting_grid[!(starting_grid$beta1==beta1_min | starting_grid$beta2==beta2_min),]
  #remove from starting_grid all lines with beta1=beta1_max and beta2=beta2_max
  starting_grid<-starting_grid[!(starting_grid$beta1==beta1_max | starting_grid$beta2==beta2_max),]
  
  # Select top 10
  top_starting_grid<-starting_grid[1:10,]
  ggplot(top_starting_grid)+geom_tile(aes(x=beta1,y=beta2,fill=scores))
  saveRDS(top_starting_grid,file=paste0("top_grid_",score_tag,".RDS"))
  top_starting_grid<-top_starting_grid[,1:2]
  }else
    {
starting_grid<-readRDS(file=paste0("../epi_fitting/top_grid_",score_tag,".RDS"))
# Discarding score column
top_starting_grid<-starting_grid[,1:2]
# Select top n.cores-1 (so that including the best fit from previous analysis there are n.cores in total)
top_starting_grid<-top_starting_grid[1:(n.cores-1),]
}
top_starting_grid<-rbind(top_starting_grid,data.frame(beta1=0.163,beta2=0.00118))
## Run Bayesian optimization
if(TYPE_PARALLEL=="PSOCK"){
cl <- makeCluster(n.cores,type="PSOCK")
parallel::clusterEvalQ(cl, library("ergm.multi"))
parallel::clusterEvalQ(cl, source("../sim_libs/sim_models.R"))
parallel::clusterExport(cl, list('HH_network_no_isolation','HH_network_w_isolation')
                        , env = environment()) 
registerDoParallel(cl)
}
if(TYPE_PARALLEL=="FORK"){
  cl <- makeCluster(n.cores,type="FORK")
  registerDoParallel(cl)
  
}
set.seed(24112022)
## Random boundaries based on top 10 best values
bounds=list(
  beta1=c(min(top_starting_grid$beta1)*(1-sample(x=1:30,size=1)/100),
          max(top_starting_grid$beta1)*(1+sample(x=1:30,size=1)/100)),    
  beta2=c(min(top_starting_grid$beta2)*(1-sample(x=1:30,size=1)/100),
          max(top_starting_grid$beta2)*(1+sample(x=1:30,size=1)/100))
)

if(DO_100){
optObj <- bayesOpt(
      FUN=SEIR_simulate_two_networks_forBOptim_Binloglik_Ncases_by_HH_averageN50
      , bounds = bounds
      ,initGrid=top_starting_grid
      , iters.n = n.cores
      , iters.k = n.cores
      ,parallel = TRUE
    )

ggplot(optObj$scoreSummary)+geom_point(aes(x=beta1,y=beta2,color=Score))
best_params<-getBestPars(optObj)
ggplot(optObj$scoreSummary)+geom_point(aes(x=beta1,y=beta2,color=1))
ggplot(optObj$scoreSummary)+geom_point(aes(x=beta1,y=beta2,color=Score))
fname<-paste0("best_params_1_runs_switch_isolation_",score_tag,".RDS")
## Save results
saveRDS(optObj,fname)
optObj<-readRDS(fname)
for(k in 1:100){
  optObj <- addIterations(optObj,iters.n=1,verbose=1,parallel=TRUE)
  ## Save output with step counter into name
  if((k%%10)==0){
    print(getBestPars(optObj))
    saveRDS(optObj,file=paste0("./BayOpt_onset/switch_isolation_",score_tag,"_k_",k,".rds"))
    }
}
ggplot(optObj$scoreSummary)+geom_point(aes(x=beta1,y=beta2,color=Score))
fname<-paste0("best_params_100_runs_switch_isolation_",score_tag,".RDS")
ggplot(optObj$scoreSummary)+geom_point(aes(x=beta1,y=beta2,color=Score))
saveRDS(optObj,fname)
}
stopCluster(cl)

set.seed(24112022)
if(TYPE_PARALLEL=="PSOCK"){
  cl <- makeCluster(n.cores,type="PSOCK")
  parallel::clusterEvalQ(cl, library("ergm.multi"))
  parallel::clusterEvalQ(cl, source("../sim_libs/sim_models.R"))
  parallel::clusterExport(cl, list('HH_network_no_isolation','HH_network_w_isolation')
                          , env = environment()) 
  registerDoParallel(cl)
}
if(TYPE_PARALLEL=="FORK"){
  cl <- makeCluster(n.cores,type="FORK")
  registerDoParallel(cl)
  
}

N_more<-200 # Number of additional steps
N_more<- N_more + 100 # For how is defined the 100 steps already done need to be taken into account
if(DO_MORE){
  fname<-paste0("best_params_100_runs_switch_isolation_",score_tag,".RDS")
  optObj<-readRDS(fname)
  for(k in 101:N_more){
    optObj <- addIterations(optObj,iters.n=1,verbose=1,parallel=TRUE)
    ## Save output with step counter into name
    if((k%%10)==0){
      print(getBestPars(optObj))
      saveRDS(optObj,file=paste0("./BayOpt_onset/switch_isolation_",score_tag,"_k_",k,".rds"))
    }
  }
  fname<-paste0("best_params_final_switch_isolation_",score_tag,".RDS")
  saveRDS(optObj,fname)
  optObj<-readRDS(fname)
  
}
stopCluster(cl)


################################################################################
#### Profile likelihood
################################################################################
fname<-paste0("best_params_final_switch_isolation_",score_tag,".RDS")
optObj<-readRDS(fname)
best_params<-getBestPars(optObj)

### Profile for b1
set.seed(24112022)
SIM_PROFILE_BETA1<-FALSE
Nvalues_profile<-10
if(!Nvalues_profile%%2){Nvalues_profile<-Nvalues_profile+1}
beta1_values<-seq(from=best_params$beta1-0.05,to=best_params$beta1+0.05,length.out=Nvalues_profile)
if(SIM_PROFILE_BETA1){
if(exists("prof_b1_df")){rm(prof_b1_df)}
for(i in (1:Nvalues_profile)){
  print(i)
  dum_list<-SEIR_simulate_two_networks_forBOptim_Binloglik_Ncases_by_HH_averageN200(
                                       beta1=beta1_values[i],
                                       beta2=best_params$beta2
  )
  
  dum_list$beta1<-beta1_values[i]
  dum_list<-as.data.frame(dum_list)
  rownames(dum_list)<-NULL
  if(!exists("prof_b1_df")){prof_b1_df<-dum_list}else{prof_b1_df<-rbind(prof_b1_df,dum_list)}
}
ggplot(prof_b1_df)+geom_point(mapping=aes(x=beta1,y=Score))
write.csv(prof_b1_df,paste0("Profile_likelihood_beta1_200sims_",score_tag,".csv"),row.names = FALSE)
}else{
  prof_b1_df<-read.csv(paste0("Profile_likelihood_beta1_200sims_",score_tag,".csv"))
}
prof_b1_df$Score_rel<-prof_b1_df$Score-max(prof_b1_df$Score)
# Create example dataset
x <- prof_b1_df$beta1
y <- prof_b1_df$Score_rel
df_prof_beta1 <- data.frame(x, y)

# Fit 4-degree polynomial
fit <- lm(y ~ poly(x, 2), data = df_prof_beta1)

# Create sequence of x values for plotting
x_seq <- seq(min(df_prof_beta1$x), max(df_prof_beta1$x), length.out = 1000)
# Compute predicted y values using the fitted polynomial
y_pred <- predict(fit, newdata = data.frame(x = x_seq))
df_pred<-data.frame(x=x_seq,y=y_pred)
write.csv(df_pred,paste0("Profile_likelihood_beta1_200sims_",score_tag,"_predicted.csv"),row.names = FALSE)
# Plot data points and polynomial fit
ggplot()+geom_line(data=df_pred,mapping=aes(x=x,y=y))+geom_point(data=df_prof_beta1,mapping=aes(x=x,y=y),color="red")+
  labs(
    x = TeX(r"( $\beta_{hh}$ )"),
    y = TeX(r"( LogLikelihood($\beta_{hh}$)-LogLikelihood($\hat{\beta_{hh}}$) )")
    ) 
ggsave(paste0("Profile_likelihood_beta_1_",score_tag,".png"),width = 10,height=6)  


new_y<-2*(max(y_pred)-y_pred)
df_error_b1<-data.frame(x=x_seq,Wilks_stat=new_y,p0=new_y<qchisq(p = 0.95, df = 1))
# Find the x of last index where condition is TRUE
first_true_index <- which(df_error_b1$p0)[1]
lower <- df_error_b1$x[first_true_index]
# Find the x of last index where condition is TRUE
last_true_index <- tail(which(df_error_b1$p0), 1)
upper <- df_error_b1$x[last_true_index]
# Find the x of the max value
max_index<-which(df_error_b1$Wilks_stat==0)
value <- df_error_b1$x[max_index]
# print values
lower
upper 
value
# lower=0.1316324
# value=0.1650658
# upper=0.1985994
best_value_b1<-value
lower_value_b1<-lower
upper_value_b1<-upper

### Profile for b2##############################################################
################################################################################
fname<-paste0("best_params_final_switch_isolation_",score_tag,".RDS")
optObj<-readRDS(fname)
best_params<-getBestPars(optObj)
SIM_PROFILE_BETA2<-FALSE
set.seed(24112022)

Nvalues_profile<-10
beta2_values<-seq(from=0.0008,to=0.0012,length.out=Nvalues_profile)
if(!Nvalues_profile%%2){Nvalues_profile<-Nvalues_profile+1}

if(best_params$beta2 %notin% beta2_values){beta2_values<-c(beta2_values,best_params$beta2)}
if(SIM_PROFILE_BETA2){
if(exists("prof_b2_df")){rm(prof_b2_df)}
for(i in (1:Nvalues_profile)){
  print(i)
  dum_list<-SEIR_simulate_two_networks_forBOptim_Binloglik_Ncases_by_HH_averageN200(
    beta1=best_params$beta1,
    beta2=beta2_values[i]
    )
  dum_list$beta2<-beta2_values[i]
  dum_list<-as.data.frame(dum_list)
  rownames(dum_list)<-NULL
  if(!exists("prof_b2_df")){prof_b2_df<-dum_list}else{prof_b2_df<-rbind(prof_b2_df,dum_list)}
}

ggplot(prof_b2_df)+geom_point(mapping=aes(x=beta2,y=Score))
write.csv(prof_b2_df,paste0("Profile_likelihood_beta2_200sims_",score_tag,".csv"),row.names = FALSE)

}else{
  prof_b2_df<-read.csv(paste0("Profile_likelihood_beta2_200sims_",score_tag,".csv"))
}



prof_b2_df$Score_rel<-prof_b2_df$Score-max(prof_b2_df$Score)


x <- prof_b2_df$beta2
y <- prof_b2_df$Score_rel
df_prof_beta2 <- data.frame(x, y)

# Fit 2-degree polynomial
fit <- lm(y ~ poly(x, 2), data = df_prof_beta2)
# Create sequence of x values for plotting
x_seq <- seq(best_params$beta2*0.75, best_params$beta2*1.25, length.out = 1000)
x_seq <- seq(0.0006, 0.0014, length.out = 1000)
# Compute predicted y values using the fitted polynomial
y_pred <- predict(fit, newdata = data.frame(x = x_seq))
df_pred<-data.frame(x=x_seq,y=y_pred)
write.csv(df_pred,paste0("Profile_likelihood_beta2_200sims_",score_tag,"_predicted.csv"),row.names = FALSE)
# Plot data points and polynomial fit
ggplot()+geom_line(data=df_pred,mapping=aes(x=x,y=y))+geom_point(data=df_prof_beta2,mapping=aes(x=x,y=y),color="red")+
  labs(
    x = TeX(r"( $\beta_{hh}$ )"),
    y = TeX(r"( LogLikelihood($\beta_{c}$)-LogLikelihood($\hat{\beta_{c}}$) )")
  ) +xlim(min(df_prof_beta2$x),max(df_prof_beta2$x))+ylim(-1,0.2)
ggsave(paste0("Profile_likelihood_beta_2_",score_tag,".png"),width = 10,height=6)  
## Compute statistics
new_y<-2*(max(y_pred)-y_pred)
df_error_b2<-data.frame(x=x_seq,Wilks_stat=new_y,p0=new_y<qchisq(p = 0.95, df = 1))

# Find the x of last index where condition is TRUE
first_true_index <- which(df_error_b2$p0)[1]
lower <- df_error_b2$x[first_true_index]
# Find the x of last index where condition is TRUE
last_true_index <- tail(which(df_error_b2$p0), 1)
upper <- df_error_b2$x[last_true_index]
# Find the x of the max value
max_index<-which(df_error_b2$Wilks_stat==0)
value <- df_error_b2$x[max_index]
lower
upper 
value
best_value_b2<-value
lower_value_b2<-lower
upper_value_b2<-upper
# lower=0.0009896173  
# upper=0.001000577  
# value= 0.0009950972
print(paste0("Beta 1:  ",best_value_b1," ( 95% CI[",lower_value_b1,"-",upper_value_b1,"])"))
print(paste0("Beta 2:  ",best_value_b2," ( 95% CI[",lower_value_b2,"-",upper_value_b2,"])"))



################################################################################
#### "Goodness of fit"
################################################################################

set.seed(24112022) # re-set seed
### Infections by HH_size INCLUDING parameters variability 
{ 
  subject_info<-read.csv("../epi_simulation//subject_data_for_analysis_of_infection.csv")
  data_to_fit<-read.csv("../epi_fitting/time_infection_matrix_net_onset_yest.csv",header=TRUE)
  memids<-data_to_fit$X
  data_to_fit<-as.matrix(data_to_fit[,2:dim(data_to_fit)[2]]) ## Dropping the Ids name column
  data_to_fit<-base::rowSums(data_to_fit,na.rm = TRUE) ## Changing to vector of ever infected
  data_to_fit[data_to_fit>0]<-1
  df_ever_infected<-data.frame(memid=memids,ever_infected=data_to_fit)
  df_ever_infected<-merge(subject_info[,c("memid","age","hhmem_num","index")],df_ever_infected,by="memid")
  df_ever_infected<-df_ever_infected %>% filter(index==0)
  
  HH_network_no_isolation<-readRDS("../network_obj_generation/net_objects/graph_full_dataset_onset_phys_in_common_yest_onlypart_symmetric.RDS")
  HH_network_w_isolation<-readRDS("../network_obj_generation/net_objects/graph_full_dataset_yest_phys_in_common_onset_onlypart_symmetric.RDS")
  if(exists("gof_inf_df")){rm(gof_inf_df)}
  N_sims<-500
  # generate N_sims values from normal distribution with average = best_value_b2/1 and sigma= (upper_value_b2/1-lower_value_b2/1)/3
  b1_values<-rnorm(N_sims,mean=best_value_b1,sd=(upper_value_b1-lower_value_b1)/3)
  b2_values<-rnorm(N_sims,mean=best_value_b2,sd=(upper_value_b2-lower_value_b2)/3)
  
  for(i in 1:N_sims){
    print(i)
    dum_list<-SEIR_simulate_two_networks(HH_network_no_isolation,
                                         HH_network_w_isolation,
                                         beta1=b1_values[i],#beta1=best_value_b1,
                                         beta2=b2_values[i]#beta2=best_value_b2
    )
    
    dum_list$n<-i
    if(!exists("gof_inf_df")){gof_inf_df<-dum_list}else{gof_inf_df<-rbind(gof_inf_df,dum_list)}
  }
  
  fname<-paste0("GOF_including_param_uncertanty_",score_tag,".csv")
  write.csv(gof_inf_df,fname)
  gof_inf_df<-merge(subject_info[,c("memid","age","hhmem_num")],gof_inf_df,by="memid")
  gof_inf_df<-gof_inf_df %>% filter(index==0)
  N_infected_gof_hhN<-gof_inf_df %>% group_by(n,hhmem_num) %>% 
    summarise(Sum_infected = sum(ever_infected),Perc_infected=sum(ever_infected)/length(ever_infected),Potential_infected=length(ever_infected))
  
  
  
  data_N_infected_by_hhN<-df_ever_infected %>% group_by(hhmem_num) %>% 
    summarise(Sum_infected = sum(ever_infected),Perc_infected=sum(ever_infected)/length(ever_infected),Potential_infected=length(ever_infected))
  lab_hh<-paste0(seq(2,7,by=1)," (",data_N_infected_by_hhN$Potential_infected/seq(1,6,by=1),")")
  lab_Npart<-paste0(seq(2,7,by=1)," (",data_N_infected_by_hhN$Potential_infected,")")
  ggplot(data=N_infected_gof_hhN)+geom_boxplot(,outlier.size = -1,mapping = aes(y=Perc_infected,group=hhmem_num,x=hhmem_num))+
    ylab("Probability of being \n a secondary case")+xlab("Household size")+
    geom_point(data=data_N_infected_by_hhN,mapping=aes(x=hhmem_num,y=Perc_infected,color="red"))+
    scale_x_continuous("Household size (# potential infected)",labels=lab_Npart,breaks=seq(2,7,by=1))+
    scale_color_manual(values = c(red = "red"), label = c(red = "Infection Data"))
  ggsave(paste0("gof_sec_by_hh_size_",score_tag,"_including_params_var.png"),width = 10,height=6)
  
}

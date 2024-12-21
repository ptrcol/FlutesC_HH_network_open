# Rcode household contact networks
# Last update by Pietro on 11/01/2024
# Last run by Pietro on 25/11/2024
rm(list = ls())
########################################################
`%notin%` <- function(a,b) ! a %in% b  
### Automatically set working directory
if(require(rstudioapi) && isAvailable()){
current_path <- getActiveDocumentContext()$path 
setwd(dirname(current_path ))
}

library(ergm)
library(ergm.multi)
library(dplyr)
library(modelsummary)
library(kableExtra)
library(ggplot2)  ## Needed for GoF plots
source("functions_model_analysis.R")
source("functions_recursive.R")

N_parallel=8

#### Load network object ####
net_onset<-readRDS("../network_obj_generation/net_objects/graph_full_dataset_onset_phys_onlypart_symmetric.RDS")
net_yest<-readRDS("../network_obj_generation/net_objects/graph_full_dataset_yest_phys_onlypart_symmetric.RDS")
net_final<-readRDS("../network_obj_generation/net_objects/graph_full_dataset_final_phys_onlypart_symmetric.RDS")
##########################################################
##########################################################
set.seed(24112022) ## Re-set seed
DO_FIT<-FALSE
OVERWRITE<-FALSE
DO_GOF<-FALSE
tag_form_name<-"full"
names_pows<-c("none","factor","linear","quadratic")
names_pows_edge<-c("none","factor","quadratic")


#### Fitting and Gof of Onset ####
tag_dataset<-"net_onset"
net_object<-net_onset
if(DO_FIT){
for(i_edge in 3:4){
term_edges<-names_pows[i_edge]
for(i_2stars in 2:i_edge ){
  term_2stars<-names_pows[i_2stars]
  for(i_tri in 2:i_2stars){
    term_triangles<-names_pows[i_tri]
    for(interaction in c(FALSE,TRUE)){
    #for(interaction in c(FALSE)){
      ## Debug
    # term_edges="quadratic"
    # term_2stars="quadratic"
    # term_triangles="linear"
    term_interaction<-NULL
    tag_term_interaction=""
    if(interaction){
      tag_term_interaction<-"with_interaction"
      term_interaction<-c("index","age_cat")
      }
    tag_form_name<-paste0("edges_",term_edges,"_2stars_",term_2stars,"_triangles_",term_triangles,tag_term_interaction)
    name_output_RDS<-paste0("model_results/model_full_data_",tag_dataset,"_",tag_form_name,".RDS")
    formula_generated<-generate_formula(tag_dataset,
                                        int_pow_degree=term_edges,
                                        int_terms_degree=NULL,#c("as.integer(has_kid==TRUE)"), #,"as.integer(n>2)"
                                        int_pow_2_stars=term_2stars,
                                        int_terms_2_stars=NULL,
                                        int_pow_triangles=term_triangles,
                                        int_terms_triangles=NULL,
                                        int_terms_nodefactors=c("index"),
                                        int_terms_nodefactors_interaction=term_interaction
    )
    if((file.exists(name_output_RDS) & (!OVERWRITE) )){warning(paste0("File ",name_output_RDS," already exist and overwrite option is ",OVERWRITE, "\n >>> AVOIDING RUNNING THE MODE <<<"))}else{
      model<- ergm(formula_generated,control = control.ergm(MCMLE.maxit = 2,parallel=N_parallel,MCMC.interval=5000,MCMC.burnin=5000,parallel.type="PSOCK"))
      a<-summary(model)
      if(!any(is.na(a$coefficients))){
        print(">>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>")
        print(name_output_RDS)
        print(">>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>")
        set.seed(24112322)
        model<- ergm(formula_generated,control = control.ergm(MCMLE.maxit = 200,MCMC.interval=50000,MCMC.burnin=10000,MCMLE.method ="Nelder-Mead")) #
        saveRDS(model,name_output_RDS)
      }
    }
    }
  }
}
}
}

if(exists("df_aic_onset")){rm(df_aic_onset)}
for(i_edge in 2:4){
  term_edges<-names_pows[i_edge]
  names_pows
  for(i_2stars in 2:i_edge ){
    term_2stars<-names_pows[i_2stars]
    for(i_tri in 2:i_2stars ){
      term_triangles<-names_pows[i_tri]
      for(interaction in c(FALSE,TRUE)){
        tag_term_interaction=""
        if(interaction){
          tag_term_interaction<-"with_interaction"
        }
        tag_form_name<-paste0("edges_",term_edges,"_2stars_",term_2stars,"_triangles_",term_triangles,tag_term_interaction)
      if(file.exists(paste0("model_results/model_full_data_",tag_dataset,"_",tag_form_name,".RDS"))){
        model<-readRDS(paste0("model_results/model_full_data_",tag_dataset,"_",tag_form_name,".RDS"))
        a<-summary(model)
        if(!any(is.na(a$coefficients))){
        df_dum<-data.frame( "Edges"=term_edges,"2-stars"=term_2stars,"Triangles"=term_triangles,"Interaction"=interaction ,"AIC"=a$aic)
        if(!exists("df_aic_onset")){df_aic_onset<-df_dum}else{df_aic_onset<-rbind(df_aic_onset,df_dum)}
        }
      }else{print(paste0(" >>> model_results/model_full_data_",tag_dataset,"_",tag_form_name,".RDS <<< does not exist") )}
      }      
    }
  }
}


best_model_onset<-readRDS(load_best_model_name(df_aic_onset,tag_dataset))
if(DO_GOF){generate_gof(net_object,best_model_onset,1000,tag_dataset,tag_form_name)}


#### Fitting and Gof of Yesterday ####
set.seed(24112022) ## Re-set seed
tag_dataset<-"net_yest"
net_object<-net_yest
tag_form_name<-"full"
names_pows<-c("none","factor","linear","quadratic")
names_pows_edge<-c("none","factor","quadratic")
if(DO_FIT){for(i_edge in 3:4){
    term_edges<-names_pows[i_edge]
    for(i_2stars in 2:i_edge ){
      term_2stars<-names_pows[i_2stars]
      for(i_tri in 2:i_2stars ){ 
        term_triangles<-names_pows[i_tri]
        for(interaction in c(FALSE,TRUE)){
        # term_edges="quadratic"  # DEBUG
        # term_2stars="quadratic" # DEBUG
        # term_triangles="linear" # DEBUG
        # interaction=TRUE       # DEBUG
        #   
        term_interaction<-NULL
        tag_term_interaction=""
        if(interaction){
          tag_term_interaction<-"with_interaction"
          term_interaction<-c("index","age_cat")
        }
        tag_form_name<-paste0("edges_",term_edges,"_2stars_",term_2stars,"_triangles_",term_triangles,tag_term_interaction)
        name_output_RDS<-paste0("model_results/model_full_data_",tag_dataset,"_",tag_form_name,".RDS")
        formula_generated<-generate_formula(tag_dataset,
                                            int_pow_degree=term_edges,
                                            int_terms_degree=NULL,#c("as.integer(has_kid==TRUE)"), #,"as.integer(n>2)"
                                            int_pow_2_stars=term_2stars,
                                            int_terms_2_stars=NULL,
                                            int_pow_triangles=term_triangles,
                                            int_terms_triangles=NULL,
                                            int_terms_nodefactors=c("index"),
                                            int_terms_nodefactors_interaction=term_interaction
        )
        if((file.exists(name_output_RDS) & (!OVERWRITE) )){warning(paste0("File ",name_output_RDS," already exist and overwrite option is ",OVERWRITE, "\n >>> AVOIDING RUNNING THE MODE <<<"))}else{
          model<- ergm(formula_generated,control = control.ergm(MCMLE.maxit = 2,parallel=N_parallel,MCMC.interval=5000,MCMC.burnin=5000,parallel.type="PSOCK")) #
          a<-summary(model)
          if(!any(is.na(a$coefficients))){
            model<- ergm(formula_generated,control = control.ergm(MCMLE.maxit = 200,parallel=N_parallel,MCMC.interval=10000,MCMC.burnin=10000,parallel.type="PSOCK")) #
            saveRDS(model,name_output_RDS)
            warning(paste0("Saved model output on File ",name_output_RDS))
          }
        }
      } # loop on "interaction"
        
      }
    }
  }
}

## Generate table with AIC
if(exists("df_aic_yest")){rm(df_aic_yest)}
for(i_edge in 2:4){
  term_edges<-names_pows[i_edge]
  for(i_2stars in 2:i_edge ){
    term_2stars<-names_pows[i_2stars]
    for(i_tri in 2:i_2stars ){
      term_triangles<-names_pows[i_tri]
      for(interaction in c(FALSE,TRUE)){
        tag_term_interaction=""
        if(interaction){
          tag_term_interaction<-"with_interaction"
        }
        tag_form_name<-paste0("edges_",term_edges,"_2stars_",term_2stars,"_triangles_",term_triangles,tag_term_interaction)
      if(file.exists(paste0("model_results/model_full_data_",tag_dataset,"_",tag_form_name,".RDS"))){
        model<-readRDS(paste0("model_results/model_full_data_",tag_dataset,"_",tag_form_name,".RDS"))
        a<-summary(model)
        if(!any(is.na(a$coefficients))){
          df_dum<-data.frame( "Edges"=term_edges,"2-stars"=term_2stars,"Triangles"=term_triangles,"Interaction"=interaction,"AIC"=a$aic)
          if(!exists("df_aic_yest")){df_aic_yest<-df_dum}else{df_aic_yest<-rbind(df_aic_yest,df_dum)}
        }
      }else{print(paste0(" >>> model_results/model_full_data_",tag_dataset,"_",tag_form_name,".RDS <<< does not exist") )}
    }
    }
  }
}

best_model_yest<-readRDS(load_best_model_name(df_aic_yest,tag_dataset))
if(DO_GOF){generate_gof(net_object,best_model_yest,1000,tag_dataset,tag_form_name)}

#### Fitting and Gof of follow-up ####
set.seed(24112022) ## Re-set seed
tag_dataset<-"net_final"
net_object<-net_final
tag_form_name<-"full"
names_pows<-c("none","factor","linear","quadratic")
names_pows_edge<-c("none","factor","quadratic")
if(DO_FIT){
  for(i_edge in 3:4){
    term_edges<-names_pows[i_edge]
    names_pows
    for(i_2stars in 2:i_edge ){
      term_2stars<-names_pows[i_2stars]
      for(i_tri in 2:i_2stars ){ 
        term_triangles<-names_pows[i_tri]
        for(interaction in c(FALSE,TRUE)){
        term_interaction<-NULL
        tag_term_interaction=""
        
        # term_edges="quadratic"  # DEBUG
        # term_2stars="quadratic" # DEBUG
        # term_triangles="linear" # DEBUG
        # interaction=FALSE       # DEBUG
        # 
        
        if(interaction){
          tag_term_interaction<-"with_interaction"
          term_interaction<-c("index","age_cat")
        }
        tag_form_name<-paste0("edges_",term_edges,"_2stars_",term_2stars,"_triangles_",term_triangles,tag_term_interaction)
        name_output_RDS<-paste0("model_results/model_full_data_",tag_dataset,"_",tag_form_name,".RDS")
        formula_generated<-generate_formula(tag_dataset,
                                            int_pow_degree=term_edges,
                                            int_terms_degree=NULL,#c("as.integer(has_kid==TRUE)"), #,"as.integer(n>2)"
                                            int_pow_2_stars=term_2stars,
                                            int_terms_2_stars=NULL,
                                            int_pow_triangles=term_triangles,
                                            int_terms_triangles=NULL,
                                            int_terms_nodefactors=c("index"),
                                            int_terms_nodefactors_interaction=term_interaction
        )
        if((file.exists(name_output_RDS) & (!OVERWRITE) )){warning(paste0("File ",name_output_RDS," already exist and overwrite option is ",OVERWRITE, "\n >>> AVOIDING RUNNING THE MODE <<<"))}else{
          model<- ergm(formula_generated,control = control.ergm(MCMLE.maxit = 2,parallel=N_parallel,MCMC.interval=5000,MCMC.burnin=5000,parallel.type="PSOCK")) #
          a<-summary(model)
          if(!any(is.na(a$coefficients))){
            model<- ergm(formula_generated,control = control.ergm(MCMLE.maxit = 200,parallel=N_parallel,MCMC.interval=1000,MCMC.burnin=10000,parallel.type="PSOCK")) #
            saveRDS(model,name_output_RDS)
          }
        }
      }
      }
    }
  }
}
## Generate table with AIC
if(exists("df_aic_final")){rm(df_aic_final)}
for(i_edge in 2:4){
  term_edges<-names_pows[i_edge]
  for(i_2stars in 2:i_edge ){
    term_2stars<-names_pows[i_2stars]
    for(i_tri in 2:i_2stars ){
      term_triangles<-names_pows[i_tri]
      for(interaction in c(FALSE,TRUE)){
        tag_term_interaction=""
        if(interaction){tag_term_interaction<-"with_interaction"}
        tag_form_name<-paste0("edges_",term_edges,"_2stars_",term_2stars,"_triangles_",term_triangles,tag_term_interaction)
      if(file.exists(paste0("model_results/model_full_data_",tag_dataset,"_",tag_form_name,".RDS"))){
        model<-readRDS(paste0("model_results/model_full_data_",tag_dataset,"_",tag_form_name,".RDS"))
        a<-summary(model)
        if(!any(is.na(a$coefficients))){
          df_dum<-data.frame( "Edges"=term_edges,"2-stars"=term_2stars,"Triangles"=term_triangles,"Interaction"=interaction,"AIC"=a$aic)
          if(!exists("df_aic_final")){df_aic_final<-df_dum}else{df_aic_final<-rbind(df_aic_final,df_dum)}
        }
        }else{print(paste0(" >>> model_results/model_full_data_",tag_dataset,"_",tag_form_name,".RDS <<< does not exist") )}
        
    }
  }
}
}

best_model_final<-readRDS(load_best_model_name(df_aic_final,tag_dataset))
if(DO_GOF){generate_gof(net_object,best_model_final,1000,tag_dataset,tag_form_name)}




modelsummary(list(best_model_onset,best_model_yest,best_model_final),shape = term + statistic ~ model,coef_omit = "Intercept",estimate="{estimate}{stars}",
             output="kableExtra")

load_best_model_name(df_aic_onset,"net_onset")
load_best_model_name(df_aic_yest,"net_yest")
load_best_model_name(df_aic_final,"net_final")
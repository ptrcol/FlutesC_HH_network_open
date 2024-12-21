

table_co_occurrence_finer<-function(int_network){
  a<-ergm.multi::uncombine_network(int_network)
  ### Compute co-occurrency according to finer age cat
  matr_co_occurrency<-matrix(0,nrow=5,ncol=5)
  for(i_net in 1:length(a)){
    finer_age_cat<-network::get.vertex.attribute(a[[i_net]],'age_cat_finer')
    #print(paste0("finer_age_cat= ",finer_age_cat))
    ## Counting diagonals
    for(age_cat_considered in 1:5){
      ### Diagonal term
      addendum<-0
      dum<-length(which(finer_age_cat==age_cat_considered))
      if(dum>1){addendum<-dum-1
      #print(paste0("age_cat_considered= ",age_cat_considered," addendum= ",addendum))
      }
      matr_co_occurrency[age_cat_considered,age_cat_considered]<-matr_co_occurrency[age_cat_considered,age_cat_considered]+addendum
      #################
      ### off-diag term
      for(age_cat_considered_bis in (age_cat_considered+1):5){
        if(age_cat_considered_bis<6){
          if((age_cat_considered%in% finer_age_cat) &(age_cat_considered_bis%in% finer_age_cat)){
            matr_co_occurrency[age_cat_considered,age_cat_considered_bis]<-      matr_co_occurrency[age_cat_considered,age_cat_considered_bis]+1
          }
        }
      }
      ################
    }
  }
  return(matr_co_occurrency)
}

table_co_occurrence<-function(int_network){
  a<-ergm.multi::uncombine_network(int_network)
  ### Compute co-occurrency according to  age cat
  matr_co_occurrency<-matrix(0,nrow=4,ncol=4)
  for(i_net in 1:length(a)){
    age_cat<-network::get.vertex.attribute(a[[i_net]],'age_cat')
    #print(paste0("age_cat= ",age_cat))
    ## Counting diagonals
    for(age_cat_considered in 1:4){
      ### Diagonal term
      addendum<-0
      dum<-length(which(age_cat==age_cat_considered))
      if(dum>1){addendum<-dum-1
      #print(paste0("age_cat_considered= ",age_cat_considered," addendum= ",addendum))
      }
      matr_co_occurrency[age_cat_considered,age_cat_considered]<-matr_co_occurrency[age_cat_considered,age_cat_considered]+addendum
      #################
      ### off-diag term
      for(age_cat_considered_bis in (age_cat_considered+1):5){
        if(age_cat_considered_bis<6){
          if((age_cat_considered%in%age_cat) &(age_cat_considered_bis%in%age_cat)){
            matr_co_occurrency[age_cat_considered,age_cat_considered_bis]<-      matr_co_occurrency[age_cat_considered,age_cat_considered_bis]+1
          }
        }
      }
      ################
    }
  }
  return(matr_co_occurrency)
}


plot_gof_models<-function(int_model_to_plot,int_model_name,int_label_net,int_label_mod){
  if(length(grep("edges",int_model_name))){
    gof_to_plot<-gofN(GOF=~edges, object=int_model_to_plot)
    png(paste0("./gof_plots/",int_label_net,"_",int_label_mod,"_edges_residuals_vs_fitted.png"))
    plot(gof_to_plot,which=1)
    dev.off()
    png(paste0("./gof_plots/",int_label_net,"_",int_label_mod,"_edges_scale_location.png"))
    plot(gof_to_plot,which=2)
    dev.off()
    png(paste0("./gof_plots/",int_label_net,"_",int_label_mod,"_edges_QQ.png"))
    plot(gof_to_plot,which=3)
    dev.off()
  }
  if(length(grep("2star",int_model_name))){
    gof_to_plot<-gofN(GOF=~kstar(2), object=int_model_to_plot)
    png(paste0("./gof_plots/",int_label_net,"_",int_label_mod,"_2star_residuals_vs_fitted.png"))
    plot(gof_to_plot,which=1)
    dev.off()
    png(paste0("./gof_plots/",int_label_net,"_",int_label_mod,"_2star_scale_location.png"))
    plot(gof_to_plot,which=2)
    dev.off()
    png(paste0("./gof_plots/",int_label_net,"_",int_label_mod,"_2star_QQ.png"))
    plot(gof_to_plot,which=3)
    dev.off()
  }
  if(length(grep("triangles",int_model_name))){
    gof_to_plot<-gofN(GOF=~triangles, object=int_model_to_plot)
    png(paste0("./gof_plots/",int_label_net,"_",int_label_mod,"_triangles_residuals_vs_fitted.png"))
    plot(gof_to_plot,which=1)
    dev.off()
    png(paste0("./gof_plots/",int_label_net,"_",int_label_mod,"_triangles_scale_location.png"))
    plot(gof_to_plot,which=2)
    dev.off()
    png(paste0("./gof_plots/",int_label_net,"_",int_label_mod,"_triangles_QQ.png"))
    plot(gof_to_plot,which=3)
    dev.off()
  }
  
  gof_degree_model_onset<-gof(int_model_to_plot~degree)
  png(paste0("./gof_plots/",int_label_net,"_",int_label_mod,"_degree.png"))
  plot(gof_degree_model_onset)
  dev.off()
  
}


plot_ergm_model<-function(ergm_model,var_names){
  summary_ergm<-summary(ergm_model)
  coeff_df<-data.frame(summary_ergm$coefficients)
  colnames(coeff_df)<-c("Estimate","Std.Error","MCMCperc","z.value","Pval")
  coeff_df$name<-rownames(coeff_df)
  ggplotObj<-ggplot(data=coeff_df)+geom_errorbar(aes(xmin=Estimate-1*Std.Error,xmax=Estimate+1*Std.Error,y=1:length(name)))+
    geom_point(aes(x=Estimate,y=1:length(name)))+scale_y_continuous(breaks=1:length(coeff_df$name),labels=var_names)+ylab("")+
    xlab("Log odds")+geom_vline(xintercept = 1,color="black")
  return(ggplotObj)
}


# plot_ergm_models<-function(ergm_model_list,int_mod_names,var_names){
#   summary_ergm<-summary(ergm_model)
#   coeff_df<-data.frame(summary_ergm$coefficients)
#   
#   
#   colnames(coeff_df)<-c("Estimate","Std.Error","MCMCperc","z.value","Pval")
#   coeff_df$name<-rownames(coeff_df)
#   
#   ggplotObj<-ggplot(data=coeff_df)+geom_errorbar(aes(xmin=Estimate-1*Std.Error,xmax=Estimate+1*Std.Error,y=1:length(name)))+
#     geom_point(aes(x=Estimate,y=1:length(name)))+scale_y_continuous(breaks=1:length(coeff_df$name),labels=var_names)+ylab("")+
#     xlab("Log odds")+geom_vline(xintercept = 1,color="black")
#   
#   return(ggplotObj)
# }



library(latex2exp)
plot_ergm_models<-function(ergm_model_list,int_mod_names,int_var_names){
  # ergm_model_list<-list(model_onset_full,model_yest_full,model_final_full)
  # int_mod_names<-c("Onset","Enrollment","Follow-up")
  # int_var_names<-var_names_plot
  for(i_model in 1:length(ergm_model_list)){
    summary_ergm<-summary(ergm_model_list[[i_model]])
    coeff_df<-data.frame(summary_ergm$coefficients)
    colnames(coeff_df)<-c("Estimate","Std.Error","MCMCperc","z.value","Pval")
    coeff_df$model<-int_mod_names[[i_model]]
    coeff_df$name<-rownames(coeff_df)
    coeff_df$y_for_plot<-1:length(coeff_df$name)
    if(i_model==1){full_coeff_df<-coeff_df}else{full_coeff_df<-rbind(full_coeff_df,coeff_df)}
  }
  full_coeff_df$color<-"black"
  full_coeff_df$color[full_coeff_df$Pval<0.05]<-"red"
  
  full_coeff_df$significance<-FALSE
  full_coeff_df$significance[full_coeff_df$Pval<0.05]<-TRUE
  
  full_coeff_df$model<-factor( full_coeff_df$model,levels=int_mod_names)
  ggplotObj<-ggplot(data=full_coeff_df)+geom_errorbar(aes(xmin=Estimate-1*Std.Error,xmax=Estimate+1*Std.Error,y=y_for_plot,color=significance))+
    geom_point(aes(x=Estimate,y=y_for_plot,color=significance))+scale_y_continuous(breaks=1:length(coeff_df$name),labels=TeX(int_var_names))+ylab("")+
    xlab("Log odds")+geom_vline(xintercept = 1,color="black")+facet_grid(.~model)+scale_color_manual(values = c("FALSE" = "black", "TRUE" = "red"))+
    theme(legend.position = "none")
  return(ggplotObj)
}
plot_ergm_models_log<-function(ergm_model_list,int_mod_names,int_var_names){
  # ergm_model_list<-list(model_onset_full,model_yest_full,model_final_full)
  # int_mod_names<-c("Onset","Enrollment","Follow-up")
  # int_var_names<-var_names_plot
  
  
  for(i_model in 1:length(ergm_model_list)){
    summary_ergm<-summary(ergm_model_list[[i_model]])
    coeff_df<-data.frame(summary_ergm$coefficients)
    colnames(coeff_df)<-c("Estimate","Std.Error","MCMCperc","z.value","Pval")
    coeff_df$model<-int_mod_names[[i_model]]
    coeff_df$name<-rownames(coeff_df)
    coeff_df$y_for_plot<-1:length(coeff_df$name)
    if(i_model==1){full_coeff_df<-coeff_df}else{full_coeff_df<-rbind(full_coeff_df,coeff_df)}
  }
  full_coeff_df$color<-"black"
  full_coeff_df$color[full_coeff_df$Pval<0.05]<-"red"
  
  full_coeff_df$significance<-FALSE
  full_coeff_df$significance[full_coeff_df$Pval<0.05]<-TRUE
  
  full_coeff_df$model<-factor( full_coeff_df$model,levels=int_mod_names)
  min_exp<-floor(min(full_coeff_df$Estimate))-1
  min_exp<-floor(min(full_coeff_df$Estimate-1.96*full_coeff_df$Std.Error))
  max_exp<-ceiling(max(full_coeff_df$Estimate))+1
  min_exp<-ceiling(min(full_coeff_df$Estimate+1.96*full_coeff_df$Std.Error))
  list_values<-`^`(10, (min_exp:max_exp))
  #list_values<-c(0.00001,0.0001,0.001,0.01,0.1,1,10,100,1000,10000)
  #list_values<-exp(list_values)
  ggplotObj<-ggplot(data=full_coeff_df)+geom_errorbar(aes(xmin=Estimate-1.96*Std.Error,xmax=Estimate+1.96*Std.Error,y=y_for_plot,color=significance))+
    geom_point(aes(x=Estimate,y=y_for_plot,color=significance))+scale_y_continuous(breaks=1:length(coeff_df$name),labels=TeX(int_var_names))+ylab("")+
    xlab("Log odds")+geom_vline(xintercept = 0,color="black")+facet_grid(.~model,scales="free_x")+scale_color_manual(values = c("FALSE" = "black", "TRUE" = "red"))+
    theme(legend.position = "none")+
    scale_x_continuous(name="Odds of a contact", breaks=log(list_values), labels=list_values)+theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
  return(ggplotObj)
}


AC_levels <- c("1","2","3","4")

wchld <- matrix(c("", "", "","S:!S",
                  NA, "", "","S:!S",
                  NA, NA, "","S:!S",
                  NA, NA, NA,"S:!S"
), 4, 4, byrow=TRUE)

selwchld <- function(x) !is.na(wchld) & wchld==x

full_select <- matrix(c("", "", "","",
                  NA, "", "","",
                  NA, NA, "","",
                  NA, NA, NA,NA
), 4, 4, byrow=TRUE)

sel_full <- function(x) !is.na(full_select) & full_select==x


generate_formula<-function(int_dataset_name,int_pow_degree,int_terms_degree,int_pow_2_stars,int_terms_2_stars,int_pow_triangles,int_terms_triangles,int_terms_nodefactors,int_terms_nodefactors_interaction){
  form_to_return<-""
  form_to_return<-paste0(form_to_return,int_dataset_name," ~ ")
  ### Edges terms
  if(int_pow_degree=="linear"){form_to_return<-paste0(form_to_return, "N(~edges,~0+log(n)")}
  if(int_pow_degree=="quadratic"){form_to_return<-paste0(form_to_return, "N(~edges,~0+ log(n)+I(log(n)^2)")}
  #if(int_pow_degree=="factor"){form_to_return<-paste0(form_to_return, "N(~edges,~ as.factor(n)")}
  if(int_pow_degree=="factor"){form_to_return<-paste0(form_to_return, "N(~edges,~ cut(n, c(0, 2, 3, 4, 7),labels=c('2','3','4','5+'))")}
  
  if(length(int_terms_degree)>0){
    for(i in 1:length(int_terms_degree)){form_to_return<-paste0(form_to_return,"+",int_terms_degree[i])}
  }
  ### Closing brackets
  form_to_return<-paste0(form_to_return, ")")
  ### 2 stars terms
  if(int_pow_2_stars!="none"){
    if(int_pow_2_stars=="linear"){form_to_return<-paste0(form_to_return, "+N(~kstar(2),~ log(n)")}
    if(int_pow_2_stars=="quadratic"){form_to_return<-paste0(form_to_return, "+N(~kstar(2),~ log(n)+I(log(n)^2)")}
    #if(int_pow_2_stars=="factor"){form_to_return<-paste0(form_to_return, "+N(~kstar(2),~ as.factor(n)")}
    if(int_pow_2_stars=="factor"){form_to_return<-paste0(form_to_return, "+N(~kstar(2),~ cut(n, c(0, 2, 3, 4, 7),labels=c('2','3','4','5+'))")}
    if(length(int_terms_2_stars)>0){
      for(i in 1:length(int_terms_2_stars)){form_to_return<-paste0(form_to_return,"+",int_terms_2_stars[i])}
    }
    ### Closing brackets
    form_to_return<-paste0(form_to_return, ")")
  }
  ### Triangles terms
  if(int_pow_triangles!="none"){
    if(int_pow_triangles=="linear"){form_to_return<-paste0(form_to_return, "+N(~triangles,~ log(n)")}
    if(int_pow_triangles=="quadratic"){form_to_return<-paste0(form_to_return, "+N(~triangles,~ log(n)+I(log(n)^2)")}
    #if(int_pow_triangles=="factor"){form_to_return<-paste0(form_to_return, "+N(~triangles,~ as.factor(n)")}
    if(int_pow_triangles=="factor"){form_to_return<-paste0(form_to_return, "+N(~triangles,~ cut(n, c(0, 2, 3, 4, 7),labels=c('2','3','4','5+'))")}
    if(length(int_terms_triangles)>0){
      for(i in 1:length(int_terms_triangles)){form_to_return<-paste0(form_to_return,"+",int_terms_triangles[i])}
    }
    ### Closing brackets
    form_to_return<-paste0(form_to_return, ")")
  }
  
  ### Mixing terms
  form_to_return<-paste0(form_to_return, "+N(~mm('age_cat',levels=AC_levels,levels2=selwchld('')))")
  
  
  ### Node factor
  if(length(int_terms_nodefactors)>0){
    for(i in 1:length(int_terms_nodefactors)){form_to_return<-paste0(form_to_return,"+nodefactor('",int_terms_nodefactors[i],"')")}
  }
  
  ### INTERACTION among Node factor 2 nodefactor
  if(length(int_terms_nodefactors_interaction)>0){
    if(length(int_terms_nodefactors_interaction)!=2){warning("Interaction term is not of length 2")}
#    for(i in 1:length(int_terms_nodefactors_interaction)){
      form_to_return<-paste0(form_to_return,"+nodefactor('",int_terms_nodefactors_interaction[1],"'):nodefactor('",int_terms_nodefactors_interaction[2],"')")
#     }
  }
  
  return(as.formula(form_to_return))
}

generate_formula_no_log<-function(int_dataset_name,int_pow_degree,int_terms_degree,int_pow_2_stars,int_terms_2_stars,int_pow_triangles,int_terms_triangles,int_terms_nodefactors){
  form_to_return<-""
  form_to_return<-paste0(form_to_return,int_dataset_name," ~ ")
  ### Edges terms
  if(int_pow_degree=="linear"){form_to_return<-paste0(form_to_return, "N(~edges,~0+n")}
  if(int_pow_degree=="quadratic"){form_to_return<-paste0(form_to_return, "N(~edges,~0+ n+I(n^2)")}
  #if(int_pow_degree=="factor"){form_to_return<-paste0(form_to_return, "N(~edges,~ as.factor(n)")}
  if(int_pow_degree=="factor"){form_to_return<-paste0(form_to_return, "N(~edges,~ cut(n, c(0, 2, 3, 4, 7),labels=c('2','3','4','5+'))")}
  
  if(length(int_terms_degree)>0){
    for(i in 1:length(int_terms_degree)){form_to_return<-paste0(form_to_return,"+",int_terms_degree[i])}
  }
  ### Closing brackets
  form_to_return<-paste0(form_to_return, ")")
  ### 2 stars terms
  if(int_pow_2_stars!="none"){
    if(int_pow_2_stars=="linear"){form_to_return<-paste0(form_to_return, "+N(~kstar(2),~ n")}
    if(int_pow_2_stars=="quadratic"){form_to_return<-paste0(form_to_return, "+N(~kstar(2),~ n+I(n^2)")}
    #if(int_pow_2_stars=="factor"){form_to_return<-paste0(form_to_return, "+N(~kstar(2),~ as.factor(n)")}
    if(int_pow_2_stars=="factor"){form_to_return<-paste0(form_to_return, "+N(~kstar(2),~ cut(n, c(0, 2, 3, 4, 7),labels=c('2','3','4','5+'))")}
    if(length(int_terms_2_stars)>0){
      for(i in 1:length(int_terms_2_stars)){form_to_return<-paste0(form_to_return,"+",int_terms_2_stars[i])}
    }
    ### Closing brackets
    form_to_return<-paste0(form_to_return, ")")
  }
  ### Triangles terms
  if(int_pow_triangles!="none"){
    if(int_pow_triangles=="linear"){form_to_return<-paste0(form_to_return, "+N(~triangles,~ n")}
    if(int_pow_triangles=="quadratic"){form_to_return<-paste0(form_to_return, "+N(~triangles,~ n+I(n^2)")}
    #if(int_pow_triangles=="factor"){form_to_return<-paste0(form_to_return, "+N(~triangles,~ as.factor(n)")}
    if(int_pow_triangles=="factor"){form_to_return<-paste0(form_to_return, "+N(~triangles,~ cut(n, c(0, 2, 3, 4, 7),labels=c('2','3','4','5+'))")}
    if(length(int_terms_triangles)>0){
      for(i in 1:length(int_terms_triangles)){form_to_return<-paste0(form_to_return,"+",int_terms_triangles[i])}
    }
    ### Closing brackets
    form_to_return<-paste0(form_to_return, ")")
  }
  
  ### Mixing terms
  form_to_return<-paste0(form_to_return, "+N(~mm('age_cat',levels=AC_levels,levels2=selwchld('')))")
  
  ### Node factor
  if(length(int_terms_nodefactors)>0){
    for(i in 1:length(int_terms_nodefactors)){form_to_return<-paste0(form_to_return,"+nodefactor('",int_terms_nodefactors[i],"')")}
  }
  return(as.formula(form_to_return))
}

generate_formula_wsenior<-function(int_dataset_name,int_pow_degree,int_terms_degree,int_pow_2_stars,int_terms_2_stars,int_pow_triangles,int_terms_triangles,int_terms_nodefactors){
  form_to_return<-""
  form_to_return<-paste0(form_to_return,int_dataset_name," ~ ")
  ### Edges terms
  if(int_pow_degree=="linear"){form_to_return<-paste0(form_to_return, "N(~edges,~0+log(n)")}
  if(int_pow_degree=="quadratic"){form_to_return<-paste0(form_to_return, "N(~edges,~0+ log(n)+I(log(n)^2)")}
  #if(int_pow_degree=="factor"){form_to_return<-paste0(form_to_return, "N(~edges,~ as.factor(n)")}
  if(int_pow_degree=="factor"){form_to_return<-paste0(form_to_return, "N(~edges,~ cut(n, c(0, 2, 3, 4, 7),labels=c('2','3','4','5+'))")}
  
  if(length(int_terms_degree)>0){
    for(i in 1:length(int_terms_degree)){form_to_return<-paste0(form_to_return,"+",int_terms_degree[i])}
  }
  ### Closing brackets
  form_to_return<-paste0(form_to_return, ")")
  ### 2 stars terms
  if(int_pow_2_stars!="none"){
    if(int_pow_2_stars=="linear"){form_to_return<-paste0(form_to_return, "+N(~kstar(2),~ log(n)")}
    if(int_pow_2_stars=="quadratic"){form_to_return<-paste0(form_to_return, "+N(~kstar(2),~ log(n)+I(log(n)^2)")}
    #if(int_pow_2_stars=="factor"){form_to_return<-paste0(form_to_return, "+N(~kstar(2),~ as.factor(n)")}
    if(int_pow_2_stars=="factor"){form_to_return<-paste0(form_to_return, "+N(~kstar(2),~ cut(n, c(0, 2, 3, 4, 7),labels=c('2','3','4','5+'))")}
    if(length(int_terms_2_stars)>0){
      for(i in 1:length(int_terms_2_stars)){form_to_return<-paste0(form_to_return,"+",int_terms_2_stars[i])}
    }
    ### Closing brackets
    form_to_return<-paste0(form_to_return, ")")
  }
  ### Triangles terms
  if(int_pow_triangles!="none"){
    if(int_pow_triangles=="linear"){form_to_return<-paste0(form_to_return, "+N(~triangles,~ log(n)")}
    if(int_pow_triangles=="quadratic"){form_to_return<-paste0(form_to_return, "+N(~triangles,~ log(n)+I(log(n)^2)")}
    #if(int_pow_triangles=="factor"){form_to_return<-paste0(form_to_return, "+N(~triangles,~ as.factor(n)")}
    if(int_pow_triangles=="factor"){form_to_return<-paste0(form_to_return, "+N(~triangles,~ cut(n, c(0, 2, 3, 4, 7),labels=c('2','3','4','5+'))")}
    if(length(int_terms_triangles)>0){
      for(i in 1:length(int_terms_triangles)){form_to_return<-paste0(form_to_return,"+",int_terms_triangles[i])}
    }
    ### Closing brackets
    form_to_return<-paste0(form_to_return, ")")
  }
  full_select
  ### Mixing terms
  form_to_return<-paste0(form_to_return, "+N(~mm('age_cat',levels=AC_levels,levels2=sel_full('')))")
  
  
  ### Node factor
  if(length(int_terms_nodefactors)>0){
    for(i in 1:length(int_terms_nodefactors)){form_to_return<-paste0(form_to_return,"+nodefactor('",int_terms_nodefactors[i],"')")}
  }
  return(as.formula(form_to_return))
}


###############################################################################
#### Plotting degree by hh size ####
###############################################################################
split_by_hh_sizes<-function(int_uncombined_net){
  int_uncombined_net<-uncombine_network(int_uncombined_net)
  int_hh_sizes<-rep(NA,length(int_uncombined_net))
  for(i in 1:length(int_uncombined_net)){
    int_hh_sizes[i]<-length(network::get.vertex.attribute(int_uncombined_net[[i]],"age_cat"))
  }
  N_unique_sizes<-length(unique(int_hh_sizes))
  max_unique_sizes<-max(unique(int_hh_sizes))
  net_list_by_size<-list(c())
  for(i_size in 2:max_unique_sizes){
    net_list_by_size<-c(net_list_by_size,list(c()))
  }
  
  int_uncombined_net_to_deplete<-int_uncombined_net
  for(i_size in 1:max_unique_sizes){
    count<-1
    i_net<-1
    while(i_net<length(int_uncombined_net_to_deplete)){
      if(length(network::get.vertex.attribute(int_uncombined_net_to_deplete[[i_net]],"age_cat"))==i_size){
        net_list_by_size[[i_size]][[count]]<-int_uncombined_net_to_deplete[[i_net]]
        int_uncombined_net_to_deplete[[i_net]]<-NULL
        count<-count+1
      }
      i_net<-i_net+1
    }
  }
  return(net_list_by_size)
}
comp_degree_by_hh_size<-function(int_net_by_hh){
  for(i_size in 1:length(int_net_by_hh)){
    if(length(int_net_by_hh[[i_size]])>0){
      a<-combine_networks(int_net_by_hh[[i_size]])
      #print(length(a))
      c<-as.data.frame(degreedist(a))
      colnames(c)<-"degree"
      b<-data.frame("perc"= c$degree,"degree"=substring(rownames(c), first = 7, last = 7),"hh_size"=i_size)
      b$perc<-b$perc/sum(b$perc)
      if(!exists("df_deg")) {df_deg<-b}else{df_deg<-rbind(df_deg,b)}
    }
  }
  return(df_deg)
}


comp_degree<-function(int_net){
  c<-as.data.frame(degreedist(int_net))
  colnames(c)<-"degree"
  dum_df<-data.frame("perc"= c$degree,"degree"=substring(rownames(c), first = 7, last = 7))
  dum_df$perc<-dum_df$perc/sum(dum_df$perc)
  return(dum_df)
}



simulate_and_degree_dists<-function(int_model,int_N_sims){
  a<-simulate(int_model,int_N_sims)

  for(i in 1:int_N_sims){
    net<-split_by_hh_sizes(a[[i]])
    dum_df<-comp_degree_by_hh_size(net)
    dum_df$sim_N<-i
    if(i==1){df_deg_sim<-dum_df}else{df_deg_sim<-rbind(df_deg_sim,dum_df)}
  }
  for(i in 1:int_N_sims){
    print(i)
    net<-a[[i]]
    c<-as.data.frame(degreedist(net))
    colnames(c)<-"degree"
    print(rownames(c))
    dum_df<-data.frame("perc"= c$degree,"degree"=substring(rownames(c), first = 7, last = 7))
    dum_df$perc<-dum_df$perc/sum(dum_df$perc)
    dum_df$sim_N<-i
    if(i==1){df_deg_sim_aggr<-dum_df}else{df_deg_sim_aggr<-rbind(df_deg_sim_aggr,dum_df)}
  }
  
  return(list(a,as.data.frame(df_deg_sim),as.data.frame(df_deg_sim_aggr)))
}


load_summary_best_model<-function(int_df_aic,int_tag_dataset){
  int_df_aic<-df_aic_onset
  index<-which(int_df_aic$AIC==min(int_df_aic$AIC))
  int_term_edges<-int_df_aic$Edges[index]
  int_term_2stars<-int_df_aic$X2.stars[index]
  int_term_triangles<-int_df_aic$Triangles[index]
  int_tag_form_name<-paste0("edges_",int_term_edges,"_2stars_",int_term_2stars,"_triangles_",int_term_triangles)
  int_model<-readRDS(paste0("model_results/model_full_data_",int_tag_dataset,"_",int_tag_form_name,".RDS"))
  return(summary(int_model))
}

load_best_model_name<-function(int_df_aic,int_tag_dataset){
  index<-which(int_df_aic$AIC==min(int_df_aic$AIC))
  int_term_edges<-int_df_aic$Edges[index]
  int_term_2stars<-int_df_aic$X2.stars[index]
  int_term_triangles<-int_df_aic$Triangles[index]
  int_interaction<-int_df_aic$Interaction[index]
  tag_term_interaction=""
  if(int_interaction){
    tag_term_interaction<-"with_interaction"
  }
  int_tag_form_name<-paste0("edges_",int_term_edges,"_2stars_",int_term_2stars,"_triangles_",int_term_triangles,tag_term_interaction)
  return(paste0("model_results/model_full_data_",int_tag_dataset,"_",int_tag_form_name,".RDS"))
}

load_best_model_name_no_log<-function(int_df_aic,int_tag_dataset){
  index<-which(int_df_aic$AIC==min(int_df_aic$AIC))
  int_term_edges<-int_df_aic$Edges[index]
  int_term_2stars<-int_df_aic$X2.stars[index]
  int_term_triangles<-int_df_aic$Triangles[index]
  int_tag_form_name<-paste0("edges_",int_term_edges,"_2stars_",int_term_2stars,"_triangles_",int_term_triangles)
  return(paste0("model_results/model_full_data2_",int_tag_dataset,"_",int_tag_form_name,".RDS"))
}

load_best_model_name_wsenior<-function(int_df_aic,int_tag_dataset){
  index<-which(int_df_aic$AIC==min(int_df_aic$AIC))
  int_term_edges<-int_df_aic$Edges[index]
  int_term_2stars<-int_df_aic$X2.stars[index]
  int_term_triangles<-int_df_aic$Triangles[index]
  int_tag_form_name<-paste0("edges_",int_term_edges,"_2stars_",int_term_2stars,"_triangles_",int_term_triangles)
  return(paste0("model_results/model_full_data3_",int_tag_dataset,"_",int_tag_form_name,".RDS"))
}

generate_gof<-function(int_net,int_model,int_N_sims,int_tag_dataset,int_tag_form_name){
  ## Generate degree distribution (by hh and aggregated) for data
  net_byhh<-split_by_hh_sizes(int_net)
  df_deg<-comp_degree_by_hh_size(net_byhh)
  df_deg_aggr<-comp_degree(int_net)
  ## Generate degree distribution (by hh and aggregated) for model
  dum_intermediate<-simulate_and_degree_dists(int_model,int_N_sims)
  df_deg_sim<-data.frame(dum_intermediate[2])
  df_deg_sim_aggr<-data.frame(dum_intermediate[3])
  
  #print(paste0("./gof_plots/",int_tag_dataset,"_",int_tag_form_name,"_degree_by_hh_size.png"))
  
  p1<-ggplot()+geom_boxplot(data=df_deg_sim,aes(x=degree,y=perc))+
    facet_grid(rows="hh_size")+geom_point(data=df_deg,aes(x=degree,y=perc),color="red")+xlab("Degree")+ylab("Proportion")
  ggsave(paste0("./gof_plots/",int_tag_dataset,"_",int_tag_form_name,"_degree_by_hh_size.png"),plot=p1,width=8,height=12)
  
  p2<-ggplot()+geom_boxplot(data=df_deg_sim_aggr,aes(x=degree,y=perc))+
    geom_point(data=df_deg_aggr,aes(x=degree,y=perc),color="red")+xlab("Degree")+ylab("Proportion")
  ggsave(filename = paste0("./gof_plots/",int_tag_dataset,"_",int_tag_form_name,"_degree.png"),plot=p2,width=8,height=4)
  
}

pivot_quantity<-function(int_binded_list,int_var_name,int_new_var_name){
  int_summary<-int_binded_list%>% 
    pivot_wider(id_cols = all_of(int_var_name), names_from = "time_point", values_from = "display")%>% 
    dplyr::rename(Group = all_of(int_var_name)) %>%
    dplyr::mutate(Category = {{int_new_var_name}})%>%
    arrange(Group)
  
}
prop_sample_net_prop<-function(int_df_attr_netr,int_time_point,int_var_name){
  int_summary <- int_df_attr_netr%>%
    dplyr::count(.data[[int_var_name]],name ="count")%>%
    dplyr::mutate(prop = prop.table(.data[["count"]]), 3) %>%
    dplyr::mutate(display = paste(count, " (", format_n(prop * 100), "%)", sep = ""))
  int_summary$time_point<-rep(int_time_point,dim(int_summary)[1]) 
  return(int_summary)
}

format_n <- function(n) {
  formatC(n, 2, format = "f")
}

inverse_logit<-function(int_val){
  return(exp(int_val)/(1+exp(int_val)))
}

combine_models<-function(int_df_list,int_mod_names,N_models,int_var_names){
  int_full_table<-data.frame(var_name=rownames(int_df_list[[1]]))
  for(mod_name in int_mod_names){
    int_full_table[[mod_name]]<-NA
  }
  for(mod_name in int_mod_names){
    int_df_coeffs<-int_df_list[[which(int_mod_names==mod_name)]]
    for(i_line in (1:dim(int_df_coeffs)[1])){
      
      string_significance<-"   "
      if(int_df_coeffs$Pr...z..[i_line]<0.05){string_significance<-"*  "}
      if(int_df_coeffs$Pr...z..[i_line]<0.01){string_significance<-"** "}
      if(int_df_coeffs$Pr...z..[i_line]<0.001){string_significance<-"***"}
      string_value<-format(round(int_df_coeffs$Estimate[i_line],2), nsmall = 2)
      if(int_df_coeffs$Estimate[i_line]>0){string_value<-paste0(" ",string_value)}
      table_entry<-paste0(string_value," (" ,format(round(int_df_coeffs$Std..Error[i_line],2), nsmall = 2) ,")",string_significance)  
      int_full_table[[mod_name]][i_line]<-table_entry
      #print(table_entry)
    }
  }
  int_full_table$var_name<-int_var_names
  return(int_full_table)
}



combine_models_prob<-function(int_df_list,int_mod_names,N_models,int_var_names){
  int_full_table<-data.frame(var_name=rownames(int_df_list[[1]]))
  for(mod_name in int_mod_names){
    int_full_table[[mod_name]]<-NA
  }
  for(mod_name in int_mod_names){
    int_df_coeffs<-int_df_list[[which(int_mod_names==mod_name)]]
    for(i_line in (1:dim(int_df_coeffs)[1])){
      dum_prob<-inverse_logit(int_df_coeffs$Estimate[i_line])
      string_significance<-"   "
      if(int_df_coeffs$Pr...z..[i_line]<0.05){string_significance<-"*  "}
      if(int_df_coeffs$Pr...z..[i_line]<0.01){string_significance<-"** "}
      if(int_df_coeffs$Pr...z..[i_line]<0.001){string_significance<-"***"}
      string_value<-format(round(dum_prob,2), nsmall = 2)
      if(int_df_coeffs$Estimate[i_line]>0){string_value<-paste0(" ",string_value)}
      table_entry<-paste0(string_value," (" ,format(round(int_df_coeffs$Std..Error[i_line]*dum_prob*(1-dum_prob),2), nsmall = 2) ,")",string_significance)  
      int_full_table[[mod_name]][i_line]<-table_entry
      #print(table_entry)
    }
  }
  int_full_table$var_name<-int_var_names
  return(int_full_table)
}


net_attr_as_df<-function(int_list_net,int_time_point){
  names_net_attributes<-list.network.attributes(int_list_net[[1]])
  names_net_attributes<-names_net_attributes[which(substr(names_net_attributes,1,1)!=".")]  ### Remove internal objects
  list_attr_not_to_be_considered<-c("ergm")   ### Remove ergm attribute
  names_net_attributes<-names_net_attributes[which(names_net_attributes %notin% list_attr_not_to_be_considered) ]  ### Remove int objects
  int_N_attributes<-length(names_net_attributes)
  int_N_networks<-length(int_list_net)
  df_net_attr<-data.frame(matrix(ncol=int_N_attributes,nrow=int_N_networks, dimnames=list(c(), names_net_attributes)))
  df_net_attr$n_edges<-rep(NA,int_N_networks)
  i_network<-1
  i_attribute<-1
  for(i_network in (1:int_N_networks)){
    for(i_attribute in (1:int_N_attributes)){
      a<-get.network.attribute(int_list_net[[i_network]],names_net_attributes[i_attribute])
      df_net_attr[i_network,i_attribute]<-a
    }
    df_net_attr$n_edges[i_network]<-network.edgecount(int_list_net[[i_network]])
    
    
  }
  
  df_net_attr$time_point<-int_time_point
  return(df_net_attr)
}




vertex_attr_as_df<-function(int_list_net,int_time_point){
  names_vertex_attributes<-network::list.vertex.attributes(int_list_net[[1]])
  int_N_attributes<-length(names_vertex_attributes)
  int_N_networks<-length(int_list_net)
  
  ## First network: creating the data frame for the vertex attributes
  int_N_vertex_this<-network::network.size(int_list_net[[1]])
  df_net_attr<-data.frame(matrix(ncol=int_N_attributes,nrow=int_N_vertex_this, dimnames=list(c(), names_vertex_attributes)))
  for(i_attribute in (1:int_N_attributes)){
    a<-network::get.vertex.attribute(int_list_net[[1]],names_vertex_attributes[i_attribute])
    df_net_attr[[names_vertex_attributes[i_attribute]]]<-a
  }
  
  i_network<-1
  i_attribute<-1
  for(i_network in (2:int_N_networks)){
    int_N_vertex_this<-network::network.size(int_list_net[[i_network]])
    df_net_attr_this<-data.frame(matrix(ncol=int_N_attributes,nrow=int_N_vertex_this, dimnames=list(c(), names_vertex_attributes)))
    for(i_attribute in (1:int_N_attributes)){
      a<-network::get.vertex.attribute(int_list_net[[i_network]],names_vertex_attributes[i_attribute])
      df_net_attr_this[[names_vertex_attributes[i_attribute]]]<-a
    }
    df_net_attr<-rbind(df_net_attr,df_net_attr_this)
    
  }
  df_net_attr$time_point<-int_time_point
  return(df_net_attr)
}


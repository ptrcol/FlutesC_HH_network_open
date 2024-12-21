source("functions_net_features.R")
generate_net_obj<-function(int_edges,int_subject,int_non_enrolled_data,int_diary,int_label,tag_plot=TRUE){
  ### Variable definition for debugging
  # tag_plot=TRUE
  # int_label<-"blabla"
  # int_edges<-edges.yestphys
  # int_subject<-subject
  # int_non_enrolled_data<-non_enrolled_data
  # int_diary<-diary

  print(paste0("Tag plot = ",tag_plot  ))
  plot_folder<-paste0("./network_plots/by_hhid/",int_label,"/")
  print(plot_folder)
  if(!dir.exists(plot_folder)){
    dir.create(plot_folder) ## If the folder does not exist: create
  }
  
  if(!dir.exists( paste0(plot_folder,"including_non_enrolled/"))){
    dir.create( paste0(plot_folder,"including_non_enrolled/")
    ) ## If the folder does not exist: create
  }
   
  int_edges<-subset(int_edges,(int_edges$cont %in% int_subject$memid)|(int_edges$cont %in% int_non_enrolled_data$memid))
  df_densities<-data.frame("hh_size"=NA,density=NA,"hh_id"=NA)
  list_unique_hhids<-unique(int_edges$hhid)
  graph_phys<-list()
  graph_phys_symmetric<-list()
  i_hh_id<-1
  i_count_hh_ne<-1
  for(i_hh_id in 1:length(list_unique_hhids)){
    hh_id<-list_unique_hhids[i_hh_id]
    hhid.sel=(int_subject$hhid==hh_id)
    list_hh_members<-int_subject$memid[hhid.sel]
    edges_this_household<-subset(int_edges,int_edges$hhid==hh_id)
    edges_this_household<-subset(edges_this_household,(edges_this_household$cont %in% int_subject$memid)|(edges_this_household$cont %in% int_non_enrolled_data$memid) )
    list_hh_members<-c(list_hh_members,edges_this_household$part,edges_this_household$cont)
    list_hh_members<-unique(list_hh_members)
    list_hh_members<-list_hh_members[list_hh_members %in% int_subject$memid| list_hh_members %in% int_non_enrolled_data$memid]
    NHHmembers<-length(list_hh_members)
    if(dim(edges_this_household)[1]>0){
      adj_matrix<-matrix(rep(0,NHHmembers*NHHmembers),nrow = NHHmembers,ncol = NHHmembers)
      colnames(adj_matrix)<-list_hh_members
      rownames(adj_matrix)<-list_hh_members
      for(i_row in 1:dim(edges_this_household)[1]){
        rownam<-as.character(edges_this_household$cont[i_row])
        colnam<-as.character(edges_this_household$part[i_row])
        adj_matrix[rownam,colnam]<-1
      }
      adj_matrix.symmetric<-pmax(adj_matrix, t(adj_matrix))
      graph_phys[[i_count_hh_ne]]<-as.network(x = adj_matrix, # the network object
                                              directed = TRUE, # specify whether the network is directed
                                              loops = FALSE, # do we allow self ties (should not allow them)
                                              matrix.type = "adjacency" # the type of input
      )
      graph_phys[[i_count_hh_ne]]<-set_age_sex_part_and_nonenrolled_graph(graph_phys[[i_count_hh_ne]],int_subject,int_non_enrolled_data)
      graph_phys[[i_count_hh_ne]] %v% "hhsize" <-  rep(edges_this_household$hhsize[1],NHHmembers)
      
      graph_phys_symmetric[[i_count_hh_ne]]<-as.network(x = adj_matrix.symmetric, # the network object
                                                        directed = FALSE, # specify whether the network is directed
                                                        loops = FALSE, # do we allow self ties (should not allow them)
                                                        matrix.type = "adjacency" # the type of input
      )
      graph_phys_symmetric[[i_count_hh_ne]]<-set_age_sex_part_and_nonenrolled_graph(graph_phys_symmetric[[i_count_hh_ne]],int_subject,int_non_enrolled_data)
      graph_phys_symmetric[[i_count_hh_ne]] %v% "hhsize" <-  rep(edges_this_household$hhsize[1],NHHmembers)
      if(tag_plot==TRUE){
        plot_ohe_HH(graph_phys_symmetric[[i_count_hh_ne]],paste0(plot_folder,"including_non_enrolled/"),hh_id,NHHmembers)
      }
      
      
      df_densities_line<-data.frame("hh_size"=edges_this_household$hhsize[1],"hh_id"=hh_id,density=network.density(graph_phys_symmetric[[i_count_hh_ne]]))
      df_densities<-rbind(df_densities,df_densities_line)
      
      i_count_hh_ne<-i_count_hh_ne+1
    }
  }
  df_densities<-na.omit(df_densities)
  
  
  HHNet <- Networks(graph_phys,directed=FALSE) # Construct a multinetwork LHS. Note that it should be treated as read-only from this point on.
  set.network.attribute(HHNet,"directed",FALSE)
  adj_matrix<-as.matrix.network.adjacency(HHNet)
  fname<-paste0("./net_objects/graph_",int_label,".RDS")
  saveRDS(HHNet,fname) ; fname<-NULL
  
  HHNet_symmetric <- Networks(graph_phys_symmetric,directed=FALSE) # Construct a multinetwork LHS. Note that it should be treated as read-only from this point on.
  set.network.attribute(HHNet_symmetric,"directed",FALSE)
  fname<-paste0("./net_objects/graph_",int_label,"_symmetric.RDS")
  saveRDS(HHNet_symmetric,fname); fname<-NULL
  
  fname<-paste0("./densities/densities_",int_label,"_symmetric.RDS")
  saveRDS(df_densities,fname) ; fname<-NULL
  
  
  ######################################################################## 
  ####### Network only with enrolled participants 
  ######################################################################## 
  df_densities<-data.frame("hh_size"=NA,density=NA,"hh_id"=NA,"has_kid"=NA)
  ### Removing contacts with non-enrolled individuals
  int_edges<-subset(int_edges,int_edges$cont %in% int_diary$memid)
  int_edges<-subset(int_edges,int_edges$part %in% int_diary$memid)
  list_unique_hhids<-unique(int_edges$hhid)
  graph_phys_onlypart<-list()
  graph_phys_onlypart_symmetric<-list()
  i_hh_id<-1
  i_count_hh_ne<-1
  for(i_hh_id in 1:length(list_unique_hhids)){
    hh_id<-list_unique_hhids[i_hh_id]
    hhid.sel=(int_subject$hhid==hh_id)
    list_hh_members<-int_subject$memid[hhid.sel]
    edges_this_household<-subset(int_edges,int_edges$hhid==hh_id)
    NHHmembers<-length(list_hh_members)
    if(dim(edges_this_household)[1]>0){
      NHHmembers<-length(list_hh_members)
      adj_matrix<-matrix(rep(0,NHHmembers*NHHmembers),nrow = NHHmembers,ncol = NHHmembers)
      colnames(adj_matrix)<-list_hh_members
      rownames(adj_matrix)<-list_hh_members
      for(i_row in 1:dim(edges_this_household)[1]){
        rownam<-as.character(edges_this_household$cont[i_row])
        colnam<-as.character(edges_this_household$part[i_row])
        adj_matrix[rownam,colnam]<-1
      }
      adj_matrix.symmetric<-pmax(adj_matrix, t(adj_matrix))
      
      graph_phys_onlypart[[i_count_hh_ne]]<-as.network(x = adj_matrix, # the network object
                                                       directed = TRUE, # specify whether the network is directed
                                                       loops = FALSE, # do we allow self ties (should not allow them)
                                                       matrix.type = "adjacency" # the type of input
      )
      source("functions_net_features.R")
      graph_phys_onlypart[[i_count_hh_ne]]<-set_features_network(graph_phys_onlypart[[i_count_hh_ne]],int_subject)
      graph_phys_onlypart[[i_count_hh_ne]] %v% "hhsize" <-  rep(edges_this_household$hhsize[1],NHHmembers)
      graph_phys_onlypart_symmetric[[i_count_hh_ne]]<-as.network(x = adj_matrix.symmetric, # the network object
                                                                 directed = FALSE, # specify whether the network is directed
                                                                 loops = FALSE, # do we allow self ties (should not allow them)
                                                                 matrix.type = "adjacency" # the type of input
      )
      graph_phys_onlypart_symmetric[[i_count_hh_ne]]<-set_features_network(graph_phys_onlypart_symmetric[[i_count_hh_ne]],int_subject)
      graph_phys_onlypart_symmetric[[i_count_hh_ne]] %v% "hhsize" <-  rep(edges_this_household$hhsize[1],NHHmembers)
      if(tag_plot==TRUE){
        plot_ohe_HH(graph_phys_onlypart_symmetric[[i_count_hh_ne]],plot_folder,hh_id,NHHmembers)
      }
      df_densities_line<-data.frame("hh_size"=edges_this_household$hhsize[1],"hh_id"=hh_id,density=network.density(graph_phys_onlypart_symmetric[[i_count_hh_ne]]),"has_kid"=get.network.attribute(graph_phys_onlypart_symmetric[[i_count_hh_ne]],"has_kid") )
      df_densities<-rbind(df_densities,df_densities_line)
      i_count_hh_ne<-i_count_hh_ne+1
    }
  }
  df_densities<-na.omit(df_densities)
  
  HHNet_onlypart<-NULL
  HHNet_onlypart <- Networks(graph_phys_onlypart,directed=FALSE) # Construct a multinetwork LHS. Note that it should be treated as read-only from this point on.
  set.network.attribute(HHNet_onlypart,"directed",FALSE)
  fname<-paste0("./net_objects/graph_",int_label,"_onlypart.RDS")
  saveRDS(HHNet_onlypart,fname) ; fname<-NULL
  
  HHNet_onlypart_symmetric<-NULL
  HHNet_onlypart_symmetric <- Networks(graph_phys_onlypart_symmetric,directed=FALSE) # Construct a multinetwork LHS. Note that it should be treated as read-only from this point on.
  set.network.attribute(HHNet_onlypart_symmetric,"directed",FALSE)
  fname<-paste0("./net_objects/graph_",int_label,"_onlypart_symmetric.RDS")
  saveRDS(HHNet_onlypart_symmetric,fname) ; fname<-NULL
  
  fname<-paste0("./densities/densities_",int_label,"_onlypart_symmetric.RDS")
  saveRDS(df_densities,fname);fname<-NULL
  
  ######################################################################## 
  ####### Network only with enrolled participants 
  ####### and restricting to HH size of three or more
  ######################################################################## 
  df_densities<-data.frame("hh_size"=NA,density=NA,"hh_id"=NA,"has_kid"=NA)
  ## Removing household of size 2
  int_edges_original<-int_edges ## Storing this to match the id of the size-larger-than-3 household with the original dataset
  #int_edges<-subset(int_edges,int_edges$hhsize>2)
  list_unique_hhids<-unique(int_edges$hhid)
  graph_phys_onlypart<-list()
  graph_phys_onlypart_symmetric<-list()
  i_hh_id<-12
  i_count_hh_ne<-1
  for(i_hh_id in 1:length(list_unique_hhids)){
    hh_id<-list_unique_hhids[i_hh_id]
    hh_id_original<-int_edges_original$hhid[which(int_edges_original$hhid==hh_id)[1]]
    hhid.sel=(int_subject$hhid==hh_id)
    list_hh_members<-int_subject$memid[hhid.sel]
    edges_this_household<-subset(int_edges,int_edges$hhid==hh_id)
    #edges_this_household<-subset(edges_this_household,(edges_this_household$cont %in% int_subject$memid)|(edges_this_household$cont %in% int_non_enrolled_data$memid) )
    #list_hh_members<-c(list_hh_members,edges_this_household$part,edges_this_household$cont)
    #list_hh_members<-unique(list_hh_members)
    #list_hh_members<-list_hh_members[list_hh_members %in% int_subject$memid| list_hh_members %in% int_non_enrolled_data$memid]
    NHHmembers<-length(list_hh_members)
    if(NHHmembers>2){
      if(dim(edges_this_household)[1]>0){
          adj_matrix<-matrix(rep(0,NHHmembers*NHHmembers),nrow = NHHmembers,ncol = NHHmembers)
          colnames(adj_matrix)<-list_hh_members
          rownames(adj_matrix)<-list_hh_members
          for(i_row in 1:dim(edges_this_household)[1]){
            rownam<-as.character(edges_this_household$cont[i_row])
            colnam<-as.character(edges_this_household$part[i_row])
            adj_matrix[rownam,colnam]<-1
          }
          adj_matrix.symmetric<-pmax(adj_matrix, t(adj_matrix))
          
          graph_phys_onlypart[[i_count_hh_ne]]<-as.network(x = adj_matrix, # the network object
                                                           directed = TRUE, # specify whether the network is directed
                                                           loops = FALSE, # do we allow self ties (should not allow them)
                                                           matrix.type = "adjacency" # the type of input
          )
          source("functions_net_features.R")
          graph_phys_onlypart[[i_count_hh_ne]]<-set_features_network(graph_phys_onlypart[[i_count_hh_ne]],int_subject)
          graph_phys_onlypart[[i_count_hh_ne]] %v% "hhsize" <-  rep(edges_this_household$hhsize[1],NHHmembers)
          graph_phys_onlypart_symmetric[[i_count_hh_ne]]<-as.network(x = adj_matrix.symmetric, # the network object
                                                                     directed = FALSE, # specify whether the network is directed
                                                                     loops = FALSE, # do we allow self ties (should not allow them)
                                                                     matrix.type = "adjacency" # the type of input
          )
          graph_phys_onlypart_symmetric[[i_count_hh_ne]]<-set_features_network(graph_phys_onlypart_symmetric[[i_count_hh_ne]],int_subject)
          graph_phys_onlypart_symmetric[[i_count_hh_ne]] %v% "hhsize" <-  rep(edges_this_household$hhsize[1],NHHmembers)
          
          if(tag_plot==TRUE){
            plot_ohe_HH(graph_phys_onlypart_symmetric[[i_count_hh_ne]],plot_folder,hh_id_original,NHHmembers,int_additional_tag="hhsize3plus")
          }
          
          df_densities_line<-data.frame("hh_size"=edges_this_household$hhsize[1],"hh_id"=hh_id,density=network.density(graph_phys_onlypart_symmetric[[i_count_hh_ne]]),"has_kid"=get.network.attribute(graph_phys_onlypart_symmetric[[i_count_hh_ne]],"has_kid") )
          df_densities<-rbind(df_densities,df_densities_line)
          i_count_hh_ne<-i_count_hh_ne+1
      }
    }
  }
  df_densities<-na.omit(df_densities)
  
  HHNet_onlypart<-NULL
  HHNet_onlypart <- Networks(graph_phys_onlypart,directed=FALSE) # Construct a multinetwork LHS. Note that it should be treated as read-only from this point on.
  set.network.attribute(HHNet_onlypart,"directed",FALSE)
  fname<-paste0("./net_objects/graph_",int_label,"_onlypart_hhsize3plus.RDS")
  saveRDS(HHNet_onlypart,fname) ; fname<-NULL
  
  HHNet_onlypart_symmetric<-NULL
  HHNet_onlypart_symmetric <- Networks(graph_phys_onlypart_symmetric,directed=FALSE) # Construct a multinetwork LHS. Note that it should be treated as read-only from this point on.
  set.network.attribute(HHNet_onlypart_symmetric,"directed",FALSE)
  fname<-paste0("./net_objects/graph_",int_label,"_onlypart_symmetric_hhsize3plus.RDS")
  saveRDS(HHNet_onlypart_symmetric,fname) ; fname<-NULL
  
  fname<-paste0("./densities/densities_",int_label,"_onlypart_symmetric_hhsize3plus.RDS")
  saveRDS(df_densities,fname);fname<-NULL
}

generate_net_obj_full_data<-function(int_edges,int_subject,int_diary,int_label,tag_plot=TRUE){
  ### Variable definition for debugging
  # tag_plot=TRUE
  # int_label<-"blabla"
  # int_edges<-edges.yestphys
  # int_subject<-subject
  # int_non_enrolled_data<-non_enrolled_data
  # int_diary<-diary

  print(paste0("Tag plot = ",tag_plot  ))
  plot_folder<-paste0("./network_plots/by_hhid/",int_label,"/")
  print(plot_folder)
  if(!dir.exists(plot_folder)){
    dir.create(plot_folder) ## If the folder does not exist: create
  }
  
  if(!dir.exists( paste0(plot_folder,"including_non_enrolled/"))){
    dir.create( paste0(plot_folder,"including_non_enrolled/")
    ) ## If the folder does not exist: create
  }
  
  
  
  ######################################################################## 
  ####### Network only with enrolled participants 
  ######################################################################## 
  df_densities<-data.frame("hh_size"=NA,density=NA,"hh_id"=NA,"has_kid"=NA)
  ### Removing contacts with non-enrolled individuals
  int_edges<-subset(int_edges,int_edges$cont %in% int_diary$memid)
  int_edges<-subset(int_edges,int_edges$part %in% int_diary$memid)
  list_unique_hhids<-unique(int_edges$hhid)
  graph_phys_onlypart<-list()
  graph_phys_onlypart_symmetric<-list()
  i_hh_id<-1
  i_count_hh_ne<-1
  for(i_hh_id in 1:length(list_unique_hhids)){
    hh_id<-list_unique_hhids[i_hh_id]
    hhid.sel=(int_subject$hhid==hh_id)
    list_hh_members<-int_subject$memid[hhid.sel]
    edges_this_household<-subset(int_edges,int_edges$hhid==hh_id)
    if(dim(edges_this_household)[1]>0){
      NHHmembers<-length(list_hh_members)
      adj_matrix<-matrix(rep(0,NHHmembers*NHHmembers),nrow = NHHmembers,ncol = NHHmembers)
      colnames(adj_matrix)<-list_hh_members
      rownames(adj_matrix)<-list_hh_members
      for(i_row in 1:dim(edges_this_household)[1]){
        rownam<-as.character(edges_this_household$cont[i_row])
        colnam<-as.character(edges_this_household$part[i_row])
        adj_matrix[rownam,colnam]<-1
      }
      adj_matrix.symmetric<-pmax(adj_matrix, t(adj_matrix))
      
      graph_phys_onlypart[[i_count_hh_ne]]<-as.network(x = adj_matrix, # the network object
                                                       directed = TRUE, # specify whether the network is directed
                                                       loops = FALSE, # do we allow self ties (should not allow them)
                                                       matrix.type = "adjacency" # the type of input
      )
      source("functions_net_features.R")
      graph_phys_onlypart[[i_count_hh_ne]]<-set_features_network(graph_phys_onlypart[[i_count_hh_ne]],int_subject)
      graph_phys_onlypart[[i_count_hh_ne]] %v% "hhsize" <-  rep(edges_this_household$hhsize[1],NHHmembers)
      graph_phys_onlypart_symmetric[[i_count_hh_ne]]<-as.network(x = adj_matrix.symmetric, # the network object
                                                                 directed = FALSE, # specify whether the network is directed
                                                                 loops = FALSE, # do we allow self ties (should not allow them)
                                                                 matrix.type = "adjacency" # the type of input
      )
      graph_phys_onlypart_symmetric[[i_count_hh_ne]]<-set_features_network(graph_phys_onlypart_symmetric[[i_count_hh_ne]],int_subject)
      graph_phys_onlypart_symmetric[[i_count_hh_ne]] %v% "hhsize" <-  rep(edges_this_household$hhsize[1],NHHmembers)
      if(tag_plot==TRUE){
        plot_ohe_HH(graph_phys_onlypart_symmetric[[i_count_hh_ne]],plot_folder,hh_id,NHHmembers)
      }
      df_densities_line<-data.frame("hh_size"=edges_this_household$hhsize[1],"hh_id"=hh_id,density=network.density(graph_phys_onlypart_symmetric[[i_count_hh_ne]]),"has_kid"=get.network.attribute(graph_phys_onlypart_symmetric[[i_count_hh_ne]],"has_kid") )
      df_densities<-rbind(df_densities,df_densities_line)
      i_count_hh_ne<-i_count_hh_ne+1
    }
  }
  df_densities<-na.omit(df_densities)
  
  HHNet_onlypart<-NULL
  HHNet_onlypart <- Networks(graph_phys_onlypart,directed=FALSE) # Construct a multinetwork LHS. Note that it should be treated as read-only from this point on.
  set.network.attribute(HHNet_onlypart,"directed",FALSE)
  fname<-paste0("./net_objects/graph_",int_label,"_onlypart.RDS")
  saveRDS(HHNet_onlypart,fname) ; fname<-NULL
  
  HHNet_onlypart_symmetric<-NULL
  HHNet_onlypart_symmetric <- Networks(graph_phys_onlypart_symmetric,directed=FALSE) # Construct a multinetwork LHS. Note that it should be treated as read-only from this point on.
  set.network.attribute(HHNet_onlypart_symmetric,"directed",FALSE)
  fname<-paste0("./net_objects/graph_",int_label,"_onlypart_symmetric.RDS")
  saveRDS(HHNet_onlypart_symmetric,fname) ; fname<-NULL
  
  fname<-paste0("./densities/densities_",int_label,"_onlypart_symmetric.RDS")
  saveRDS(df_densities,fname);fname<-NULL
  
  ######################################################################## 
  ####### Network only with enrolled participants 
  ####### and restricting to HH size of three or more
  ######################################################################## 
  df_densities<-data.frame("hh_size"=NA,density=NA,"hh_id"=NA,"has_kid"=NA)
  ## Removing household of size 2
  int_edges_original<-int_edges ## Storing this to match the id of the size-larger-than-3 household with the original dataset
  #int_edges<-subset(int_edges,int_edges$hhsize>2)
  list_unique_hhids<-unique(int_edges$hhid)
  graph_phys_onlypart<-list()
  graph_phys_onlypart_symmetric<-list()
  i_hh_id<-12
  i_count_hh_ne<-1
  for(i_hh_id in 1:length(list_unique_hhids)){
    hh_id<-list_unique_hhids[i_hh_id]
    hh_id_original<-int_edges_original$hhid[which(int_edges_original$hhid==hh_id)[1]]
    hhid.sel=(int_subject$hhid==hh_id)
    list_hh_members<-int_subject$memid[hhid.sel]
    edges_this_household<-subset(int_edges,int_edges$hhid==hh_id)
    #edges_this_household<-subset(edges_this_household,(edges_this_household$cont %in% int_subject$memid)|(edges_this_household$cont %in% int_non_enrolled_data$memid) )
    #list_hh_members<-c(list_hh_members,edges_this_household$part,edges_this_household$cont)
    #list_hh_members<-unique(list_hh_members)
    #list_hh_members<-list_hh_members[list_hh_members %in% int_subject$memid| list_hh_members %in% int_non_enrolled_data$memid]
    NHHmembers<-length(list_hh_members)
    if(NHHmembers>2){
      if(dim(edges_this_household)[1]>0){
        adj_matrix<-matrix(rep(0,NHHmembers*NHHmembers),nrow = NHHmembers,ncol = NHHmembers)
        colnames(adj_matrix)<-list_hh_members
        rownames(adj_matrix)<-list_hh_members
        for(i_row in 1:dim(edges_this_household)[1]){
          rownam<-as.character(edges_this_household$cont[i_row])
          colnam<-as.character(edges_this_household$part[i_row])
          adj_matrix[rownam,colnam]<-1
        }
        adj_matrix.symmetric<-pmax(adj_matrix, t(adj_matrix))
        
        graph_phys_onlypart[[i_count_hh_ne]]<-as.network(x = adj_matrix, # the network object
                                                         directed = TRUE, # specify whether the network is directed
                                                         loops = FALSE, # do we allow self ties (should not allow them)
                                                         matrix.type = "adjacency" # the type of input
        )
        source("functions_net_features.R")
        graph_phys_onlypart[[i_count_hh_ne]]<-set_features_network(graph_phys_onlypart[[i_count_hh_ne]],int_subject)
        graph_phys_onlypart[[i_count_hh_ne]] %v% "hhsize" <-  rep(edges_this_household$hhsize[1],NHHmembers)
        graph_phys_onlypart_symmetric[[i_count_hh_ne]]<-as.network(x = adj_matrix.symmetric, # the network object
                                                                   directed = FALSE, # specify whether the network is directed
                                                                   loops = FALSE, # do we allow self ties (should not allow them)
                                                                   matrix.type = "adjacency" # the type of input
        )
        graph_phys_onlypart_symmetric[[i_count_hh_ne]]<-set_features_network(graph_phys_onlypart_symmetric[[i_count_hh_ne]],int_subject)
        graph_phys_onlypart_symmetric[[i_count_hh_ne]] %v% "hhsize" <-  rep(edges_this_household$hhsize[1],NHHmembers)
        
        if(tag_plot==TRUE){
          plot_ohe_HH(graph_phys_onlypart_symmetric[[i_count_hh_ne]],plot_folder,hh_id_original,NHHmembers,int_additional_tag="hhsize3plus")
        }
        
        df_densities_line<-data.frame("hh_size"=edges_this_household$hhsize[1],"hh_id"=hh_id,density=network.density(graph_phys_onlypart_symmetric[[i_count_hh_ne]]),"has_kid"=get.network.attribute(graph_phys_onlypart_symmetric[[i_count_hh_ne]],"has_kid") )
        df_densities<-rbind(df_densities,df_densities_line)
        i_count_hh_ne<-i_count_hh_ne+1
      }
    }
  }
  df_densities<-na.omit(df_densities)
  
  HHNet_onlypart<-NULL
  HHNet_onlypart <- Networks(graph_phys_onlypart,directed=FALSE) # Construct a multinetwork LHS. Note that it should be treated as read-only from this point on.
  set.network.attribute(HHNet_onlypart,"directed",FALSE)
  fname<-paste0("./net_objects/graph_",int_label,"_onlypart_hhsize3plus.RDS")
  saveRDS(HHNet_onlypart,fname) ; fname<-NULL
  
  HHNet_onlypart_symmetric<-NULL
  HHNet_onlypart_symmetric <- Networks(graph_phys_onlypart_symmetric,directed=FALSE) # Construct a multinetwork LHS. Note that it should be treated as read-only from this point on.
  set.network.attribute(HHNet_onlypart_symmetric,"directed",FALSE)
  fname<-paste0("./net_objects/graph_",int_label,"_onlypart_symmetric_hhsize3plus.RDS")
  saveRDS(HHNet_onlypart_symmetric,fname) ; fname<-NULL
  
  fname<-paste0("./densities/densities_",int_label,"_onlypart_symmetric_hhsize3plus.RDS")
  saveRDS(df_densities,fname);fname<-NULL
}




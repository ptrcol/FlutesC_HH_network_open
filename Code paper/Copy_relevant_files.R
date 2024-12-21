  # Rcode to copy relevant files from other folders into the subfolder "data"
  # Last update by Pietro on 26/11/2024
  rm(list = ls())
  ########################################################
  `%notin%` <- function(a,b) ! a %in% b  
  ### Automatically set working directory
  if(require(rstudioapi) && isAvailable()){
    current_path <- getActiveDocumentContext()$path 
    setwd(dirname(current_path ))
  }
  DO_OVERWRITE<-TRUE # overwrite files if they already exist
  
  
  # Define a list of files and folders to copy from
  list_origin_file_name<-list()
  list_origin_folder<-list()
  
  # Populate the lists
  
  ## Data for figure 2
  list_origin_file_name[1]<-"densities_full_dataset_onset_phys_onlypart_symmetric.RDS"
  list_origin_folder[1]<-"../network_obj_generation/densities/"
  list_origin_file_name[2]<-"densities_full_dataset_yest_phys_onlypart_symmetric.RDS"
  list_origin_folder[2]<-list_origin_folder[1]
  list_origin_file_name[3]<-"densities_full_dataset_final_phys_onlypart_symmetric.RDS"
  list_origin_folder[3]<-list_origin_folder[1]
  
  ## Data for figure 3 and 4
  list_origin_file_name[4]<-"model_full_data_net_onset_edges_quadratic_2stars_quadratic_triangles_linear.RDS"
  list_origin_folder[4]<-"../network_fitting/model_results/"
  list_origin_file_name[5]<-"model_full_data_net_yest_edges_quadratic_2stars_quadratic_triangles_linearwith_interaction.RDS"
  list_origin_folder[5]<-"../network_fitting/model_results/"
  list_origin_file_name[6]<-"model_full_data_net_final_edges_quadratic_2stars_quadratic_triangles_quadratic.RDS"
  list_origin_folder[6]<-"../network_fitting/model_results/"
  
  ## Data for figure 5 and 6
  fit_tag<-"main_analysis"  
  list_origin_file_name[7]<-paste0("perc_infections_by_hh_cat_",fit_tag,".csv")
  list_origin_folder[7]<-"../epi_simulation/"
  list_origin_file_name[8]<-paste0("perc_infections_by_age_",fit_tag,".csv")
  list_origin_folder[8]<-"../epi_simulation/"
  
  
  for(i in 1:length(list_origin_file_name)){
    origin_file_name<-list_origin_file_name[i]
    origin_folder<-list_origin_folder[i]
    # Name variables that are fixed  
    destination_file_name<-origin_file_name
    destination_folder<-"./data/"
    origin_name<-paste0(origin_folder,origin_file_name)
    destination_name<-paste0(destination_folder,destination_file_name)
    # copy file from other folders. Name of the file "x" and name of the folder "y"
    file.copy(from = origin_name, to = destination_name, overwrite = DO_OVERWRITE)
  }
  
  
  
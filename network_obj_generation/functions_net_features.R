


set_features_network<- function( int_net_object, dataframe_from){
  
  # int_net_object<-graph_phys_onlypart[[i_count_hh_ne]]
  # dataframe_from<-int_subject
  # 
  name_vertices<-network.vertex.names(int_net_object)
  part_characteristics<-dataframe_from
  part_characteristics <- part_characteristics[part_characteristics$memid %in% name_vertices,]
  part_characteristics_ordered <- part_characteristics[match(name_vertices, part_characteristics$memid), ]
  int_net_object %v% "hhid"<-part_characteristics_ordered$hhid
  int_net_object %v% "gender"<-part_characteristics_ordered$sex
  int_net_object %v% "age"<-part_characteristics_ordered$age
  int_net_object %v% "site"<-part_characteristics_ordered$site
  # out<-paste0("Ordered age",part_characteristics_ordered$age)
  # print(out)
  out<-paste0("Ordered cut",cut(part_characteristics_ordered$age,breaks=c(-1,17,34,64,100),labels=c(1,2,3,4))) #,labels=c("(0,18)","[18,35)","[35,65)","[65,100]")
  #print(out)
  
  int_net_object %v% "age_cat"    <- as.integer(cut(part_characteristics_ordered$age,breaks=c(-1,17,34,64,100),labels=c(1,2,3,4)))
  int_net_object %v% "age_cat_finer"    <- as.integer(cut(part_characteristics_ordered$age,breaks=c(-1,11,17,34,64,100),labels=c(1,2,3,4,5)))
  #print(cut(part_characteristics_ordered$age,breaks=c(0,17,34,64,100),labels=c("child","young adult","older adult","senior")))
   # int_net_object %v% "age_cat"    <- cut(part_characteristics_ordered$age,breaks=c(0,17,34,64,100),labels=c("child","young adult","older adult","senior"))
  hh_size_aggregated<-part_characteristics_ordered$noNA_hh_size
  int_net_object %v% "hhsize_original"<-hh_size_aggregated
  int_net_object %n% "hhsize_original"<-hh_size_aggregated[1]
  hh_size_aggregated[hh_size_aggregated>=5]<-"5+"
  int_net_object %v% "hhsize"<-hh_size_aggregated
  int_net_object %n% "hhsize"<-hh_size_aggregated[1]
  int_net_object %n% "hh_size_not2"<-hh_size_aggregated[1]!=2
  int_net_object %n% "hh_size_1or2"<-hh_size_aggregated[1]<=2
  int_net_object %n% "hh_size_3or4"<-hh_size_aggregated[1]==3|hh_size_aggregated[1]==4
  int_net_object %n% "hh_size_5plus"<-hh_size_aggregated[1]>=5
  
  int_net_object %v% "bedroom"    <-part_characteristics_ordered$bedrooms[1]
  int_net_object %n% "bedroom_density"    <-part_characteristics_ordered$bedrooms[1]/part_characteristics_ordered$noNA_hh_size[1]
  int_net_object %v% "hhtype"    <-part_characteristics_ordered$hhtype[1]
  
  
  
  int_net_object %v% "index"<-part_characteristics_ordered$index
  int_net_object %v% "school"<-part_characteristics_ordered$school
  
  #####  Vaccination:
  temp<-part_characteristics_ordered$sarsvx
  temp[is.na(temp)]<-0
  int_net_object %v% "vaccinated"<-temp
  int_net_object %v% "flu_vaccinated"<-part_characteristics_ordered$fluvx
  int_net_object %v% "flu_vaccinated_lastyr"<-part_characteristics_ordered$fluvx_lastyr

  #####  Symptoms : cough
  int_net_object %v% "cough_prior"<-part_characteristics_ordered$cough_prior
  int_net_object %v% "cough_diary_onset"<-part_characteristics_ordered$cough_diary_onset
  int_net_object %v% "cough_diary_1aft"<-part_characteristics_ordered$cough_diary_1aft
  int_net_object %v% "cough_diary_2aft"<-part_characteristics_ordered$cough_diary_2aft
  int_net_object %v% "cough_diary_3aft"<-part_characteristics_ordered$cough_diary_3aft
  int_net_object %v% "cough_diary_4aft"<-part_characteristics_ordered$cough_diary_4aft
  int_net_object %v% "cough_diary_5aft"<-part_characteristics_ordered$cough_diary_5aft
  #####  Symptoms : fever
  int_net_object %v% "fever_prior"<-part_characteristics_ordered$fever_prior
  int_net_object %v% "fever_diary_onset"<-part_characteristics_ordered$fever_diary_onset
  int_net_object %v% "fever_diary_1aft"<-part_characteristics_ordered$fever_diary_1aft
  int_net_object %v% "fever_diary_2aft"<-part_characteristics_ordered$fever_diary_2aft
  int_net_object %v% "fever_diary_3aft"<-part_characteristics_ordered$fever_diary_3aft
  int_net_object %v% "fever_diary_4aft"<-part_characteristics_ordered$fever_diary_4aft
  int_net_object %v% "fever_diary_5aft"<-part_characteristics_ordered$fever_diary_5aft
  #####  Symptoms : test_smell loss
  int_net_object %v% "taste_smell_loss_prior"<-part_characteristics_ordered$taste_smell_loss_prior
  int_net_object %v% "taste_smell_loss_onset"<-part_characteristics_ordered$taste_smell_loss_onset
  int_net_object %v% "taste_smell_loss_diary_1aft"<-part_characteristics_ordered$taste_smell_loss_diary_1aft
  int_net_object %v% "taste_smell_loss_diary_2aft"<-part_characteristics_ordered$taste_smell_loss_diary_2aft
  int_net_object %v% "taste_smell_loss_diary_3aft"<-part_characteristics_ordered$taste_smell_loss_diary_3aft
  int_net_object %v% "taste_smell_loss_diary_4aft"<-part_characteristics_ordered$taste_smell_loss_diary_4aft
  int_net_object %v% "taste_smell_loss_diary_5aft"<-part_characteristics_ordered$taste_smell_loss_diary_5aft
  
  a<-part_characteristics_ordered$fever_diary_onset
  a[is.na(a)]<-0
  b<-part_characteristics_ordered$cough_diary_onset
  b[is.na(b)]<-0
  c<-part_characteristics_ordered$taste_smell_loss_onset
  c[is.na(c)]<-0
  int_net_object %v% "any_sympt_onset"<-a|b|c
  a<-part_characteristics_ordered$fever_prior
  a[is.na(a)]<-0
  b<-part_characteristics_ordered$cough_prior
  b[is.na(b)]<-0
  c<-part_characteristics_ordered$taste_smell_loss_prior
  c[is.na(c)]<-0
  dum_any_symp_prior<-a|b|c
  dum_any_symp_prior[is.na(dum_any_symp_prior)]<-0
  int_net_object %v% "any_sympt_prior"<-dum_any_symp_prior
  
  if(any(part_characteristics_ordered$age<18) ){
    int_net_object %n% "has_kid"    <-TRUE
  }else{int_net_object %n% "has_kid"    <-FALSE}
  #### Risk status
  dum<-part_characteristics_ordered$immuno_cond
  dum[is.na(dum)]<-0
  dum[dum==9]<-0
  int_net_object %v% "immuno_compromised"<-dum
  
  dum<-part_characteristics_ordered$extreme_obese
  dum[is.na(dum)]<-0
  int_net_object %v% "extreme_obese"<-dum
  
  dum<-part_characteristics_ordered$pregnant
  dum[is.na(dum)]<-0
  int_net_object %v% "pregnant"<-dum
  
  dum<-part_characteristics_ordered$diab
  dum[is.na(dum)]<-0
  int_net_object %v% "diabetes"<-dum
  
  dum<-part_characteristics_ordered$kidney
  dum[is.na(dum)]<-0
  int_net_object %v% "kidney_condition"<-dum
  
  dum<-part_characteristics_ordered$cancer
  dum[is.na(dum)]<-0
  int_net_object %v% "cancer"<-dum
  
  return(int_net_object)
}

set_age_sex_part_and_nonenrolled_graph<-function(int_net_object,data_frame_enrolled_participants,data_frame_non_enrolled){
  name_vertices<-network.vertex.names(int_net_object)
  characteristics<-data.frame(memid=name_vertices )
  part_characteristics<-data_frame_enrolled_participants
  part_characteristics <- part_characteristics[part_characteristics$memid %in% name_vertices,]
  part_characteristics_ordered <- part_characteristics[match(name_vertices, part_characteristics$memid), ]
  
  non_enrolled_characteristics<-data_frame_non_enrolled[data_frame_non_enrolled$memid %in% name_vertices,]
  non_enrolled_characteristics_ordered <- non_enrolled_characteristics[match(name_vertices, non_enrolled_characteristics$memid), ]
  
  part_characteristics_ordered$sex[is.na(part_characteristics_ordered$sex)] = non_enrolled_characteristics_ordered$hce_sex[is.na(part_characteristics_ordered$sex)]
  part_characteristics_ordered$age[is.na(part_characteristics_ordered$age)] = non_enrolled_characteristics_ordered$hce_age[is.na(part_characteristics_ordered$age)]
  part_characteristics_ordered$index[is.na(part_characteristics_ordered$index)] = 0
  
  part_characteristics_ordered$age<-as.integer(part_characteristics_ordered$age)
  
  int_net_object %v% "hhid"   <- part_characteristics_ordered$hhid
  
  hh_size_aggregated<-part_characteristics_ordered$noNA_hh_size
  int_net_object %n% "hhsize_original"<-hh_size_aggregated[1]
  int_net_object %v% "hhsize_original"<-hh_size_aggregated
  hh_size_aggregated[hh_size_aggregated>=5]<-"5+"
  int_net_object %v% "hhsize"<-hh_size_aggregated
  int_net_object %n% "hhsize"<-hh_size_aggregated[1]
  
  int_net_object %n% "hh_size_not2"<-hh_size_aggregated[1]!=2
  int_net_object %n% "hh_size_1or2"<-hh_size_aggregated[1]<=2
  int_net_object %n% "hh_size_3or4"<-hh_size_aggregated[1]==3|hh_size_aggregated[1]==4
  int_net_object %n% "hh_size_5plus"<-hh_size_aggregated[1]>=5
  
  hh_size_cat<-1
  if(hh_size_aggregated[1]==3|hh_size_aggregated[1]==4){hh_size_cat<-2}
  if(hh_size_aggregated[1]>=5){hh_size_cat<-3}
  int_net_object %n% "hh_size_cat"<-hh_size_cat
  int_net_object %v% "bedroom"    <-part_characteristics_ordered$bedroom[1]
  int_net_object %n% "bedroom_density"    <-part_characteristics_ordered$bedroom[1]/part_characteristics_ordered$noNA_hh_size[1]
  
  int_net_object %v% "gender" <-part_characteristics_ordered$sex
  int_net_object %v% "age"    <-part_characteristics_ordered$age
    #int_net_object %v% "age_cat"    <- cut(part_characteristics_ordered$age,breaks=c(0,17,34,64,100),labels=c("(0,18)","[18,35)","[35,65)","[65,100]"))
  int_net_object %v% "age_cat"    <- as.integer(cut(part_characteristics_ordered$age,breaks=c(-1,17,34,64,100),labels=c(1,2,3,4)))
  #print( as.integer(cut(part_characteristics_ordered$age,breaks=c(-1,17,34,64,100),labels=c(1,2,3,4))))
  int_net_object %v% "age_cat_finer"    <- as.integer(cut(part_characteristics_ordered$age,breaks=c(0,11,17,34,64,100),labels=c(1,2,3,4,5)))
  #warning(as.vector(part_characteristics_ordered$sex))
  #warning(cut(part_characteristics_ordered$age,breaks=c(0,17,34,64,100),labels=c("(0,18)","[18,35)","[35,65)","[65,100]")))
  int_net_object %v% "bedroom"    <-part_characteristics_ordered$bedroom[1]
  int_net_object %v% "hhtype"    <-part_characteristics_ordered$hhtype[1]
  int_net_object %v% "index"    <-part_characteristics_ordered$index
  
  
  if(any(part_characteristics_ordered$age<18 )){
    int_net_object %n% "has_kid"    <-TRUE
    #print(part_characteristics_ordered$age)
    #print("Household has kid")
  }else{int_net_object %n% "has_kid"    <-FALSE
  #print("Household DOES NOT have kid")
  }
  
  
  return(int_net_object)
}


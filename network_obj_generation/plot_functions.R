plot_full_network<-function(filename,net_object){
  png(filename,width= 3500,height = 3500)
  plot(net_object, 
       edge.arrow.size=.01,
       vertex.size=.8
       ,layout=layout.fruchterman.reingold
  )
  dev.off()
  
}


plot_single_households<-function(graph_object,folder,nametag){
  sub_graphs  <- decompose.graph(graph_object)
  list_densities<-list()
  list_hh_size<-list()
  for(i in 1:length(sub_graphs)){
    png(paste0(folder,"hhid_",nametag,"_",i,".png"))
    a<-V(sub_graphs[[i]])$hhid[1]
    b<-edge_density(sub_graphs[[i]], loops = FALSE)
    c<-gorder(sub_graphs[[i]])
    list_densities[i]<-b
    list_hh_size[i]<-c
    plot(sub_graphs[[i]]
         ,main=paste0("HH_id= ",a," density= ",b)
    )
    dev.off()
  }
  return(df_dens<-data.frame("hh_size"=unlist(list_hh_size),"density"=unlist(list_densities)))
}






plot_ohe_HH<-function(graph_object_one_hh,int_folder,int_hhid,int_hh_size,int_additional_tag=FALSE){
  
  fname<-paste0(int_folder,"hhid_",int_hhid,".png")
  if(int_additional_tag!=FALSE){fname<-paste0(int_folder,"hhid_",int_hhid,"_",int_additional_tag,".png")}
  png(fname)
  #print(fname)
  
  
  plot.network(graph_object_one_hh,label=network.vertex.names(graph_object_one_hh)) 
  title(paste("hhid",int_hhid,";",
              "hhsize",int_hh_size
              )
  )
  dev.off()
}

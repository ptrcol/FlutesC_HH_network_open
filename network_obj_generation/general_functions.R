library(EpiModel)
dyads.counting=function(tmp)
{
  np=sum(apply(tmp,1,sum)!=0)
  nh=sum(apply(tmp,1,sum)==0)
  pos.dyads=np*(np+nh)-np
  obs.dyads=sum(tmp)
  return(list(pos.dyads=pos.dyads,obs.dyads=obs.dyads,prop.dyads=obs.dyads/pos.dyads))
}

df_with_net_attributes<-function(net_object){
  names_vertex_attributes<-network::list.vertex.attributes(net_object)
  df_net<-data.frame("id"=network::get.vertex.attribute(net_object,names_vertex_attributes[1]))
  for (i in 2:length(names_vertex_attributes)){
    name_col<-names_vertex_attributes[i]
    df_net[[names_vertex_attributes[i]]]<-network::get.vertex.attribute(net_object,names_vertex_attributes[i])
  }
  df_net$degree<-get_degree(net_object)
  return(df_net)
}

# Rcode household contact networks
# Last update: v20211123
# Last update by Pietro on 26/11/2024
rm(list = ls())
########################################################
`%notin%` <- function(a,b) ! a %in% b  
### Automatically set working directory
if(require(rstudioapi) && isAvailable()){
  current_path <- getActiveDocumentContext()$path 
  setwd(dirname(current_path ))
}
library(purrr)
library(ggpubr)
library(ggplot2)
source("functions_model_analysis.R")
source("functions_IO.R")


model_onset_full<-readRDS("./data/model_full_data_net_onset_edges_quadratic_2stars_quadratic_triangles_linear.RDS")
model_yest_full<-readRDS("./data/model_full_data_net_yest_edges_quadratic_2stars_quadratic_triangles_linearwith_interaction.RDS")
model_final_full<-readRDS("./data/model_full_data_net_final_edges_quadratic_2stars_quadratic_triangles_quadratic.RDS")


##### Size effects fucntions ####
nseff <- function(coef, interval=c(2,7),intercept=TRUE,...){
  
  f <- function(x) exp(map_dbl(x, ~sum(log(.)^((1-intercept):d)*coef)))
  par(mar=c(3,4,3,0))
  plot(xr, f(xr), xlab="$n$", ylab="Conditional odds ratio",...)
  curve(f(x), from=interval[1], to=interval[2], add=TRUE)
}
make_df_size_effect <- function(coef, interval=c(2,7),intercept=TRUE,df_tag,int_term,...){
  d <- length(coef)-intercept
  xr <- interval[1]:interval[2]
  f <- function(x) exp(map_dbl(x, ~sum(log(.)^((1-intercept):d)*coef)))
  
  int_dum_df<-data.frame("n"=xr, "cond_odds"=f(xr),"label"=df_tag,"term"=int_term)
  return(int_dum_df)
}


poly_ns_eff_edges = c("N(log(n))~edges", "N(I(log(n)^2))~edges")
poly_ns_eff_2stars = c("N(1)~kstar2", "N(log(n))~kstar2", "N(I(log(n)^2))~kstar2")
poly_ns_eff_triangles_lin = c("N(1)~triangle", "N(log(n))~triangle")
poly_ns_eff_triangles_quadr = c("N(1)~triangle", "N(log(n))~triangle", "N(I(log(n)^2))~triangle")


df_size_effect<-make_df_size_effect(coef(model_onset_full)[poly_ns_eff_edges], main="edges", intercept=FALSE,df_tag = "Onset",int_term="Edges")
df_size_effect<-rbind(df_size_effect,make_df_size_effect(coef(model_yest_full)[poly_ns_eff_edges], main="edges", intercept=FALSE,df_tag = "Enrollment",int_term="Edges"))
df_size_effect<-rbind(df_size_effect,make_df_size_effect(coef(model_final_full)[poly_ns_eff_edges], main="edges", intercept=FALSE,df_tag = "Follow-up",int_term="Edges"))


df_size_effect<-rbind(df_size_effect,make_df_size_effect(coef(model_onset_full)[poly_ns_eff_2stars],  main="2-stars", interval=c(3,7),df_tag = "Onset",int_term="2-stars"))
df_size_effect<-rbind(df_size_effect,make_df_size_effect(coef(model_yest_full)[poly_ns_eff_2stars],  main="2-stars", interval=c(3,7),df_tag = "Enrollment",int_term="2-stars"))
df_size_effect<-rbind(df_size_effect,make_df_size_effect(coef(model_final_full)[poly_ns_eff_2stars],  main="2-stars", interval=c(3,7),df_tag = "Follow-up",int_term="2-stars"))

df_size_effect<-rbind(df_size_effect,make_df_size_effect(coef(model_onset_full)[poly_ns_eff_triangles_lin],  main="triangles", interval=c(3,7),df_tag = "Onset",int_term="Triangles"))
df_size_effect<-rbind(df_size_effect,make_df_size_effect(coef(model_yest_full)[poly_ns_eff_triangles_lin],  main="triangles", interval=c(3,7),df_tag = "Enrollment",int_term="Triangles"))
df_size_effect<-rbind(df_size_effect,make_df_size_effect(coef(model_final_full)[poly_ns_eff_triangles_quadr],  main="triangles", interval=c(3,7),df_tag = "Follow-up",int_term="Triangles"))

df_size_effect$label<-factor(df_size_effect$label,levels=c("Onset","Enrollment","Follow-up"))
df_size_effect$term<-factor(df_size_effect$term,levels=c("Edges","2-stars","Triangles"))
ggplot(df_size_effect)+ geom_point(aes(x=n,y=cond_odds))  +facet_grid(cols=vars(label),rows = vars(term))


### Comment when generating figure "interactions.png"
# df_size_effect$term<-factor(df_size_effect$term,levels=c("Edges","2-stars","Triangles"), labels=c(
#   "Physical contact between two people (edge)",
#   "Two people have physical contact with a third person  \n but not with each other (2-star)",
#   "Three people each have physical contact \n with each other (triangle)"
# ))


ggplot(df_size_effect)+ geom_point(aes(x=n,y=cond_odds,color=label))+ geom_line(aes(x=n,y=cond_odds,color=label))   +facet_wrap(~term,nrow=1,ncol=3,scales="free")+
  xlab("Hosehold size")+ylab("Conditional odds ratio")+ labs(color = "Time-point")
ggsave("Figure4.png",height=4,width=12)
p_size<-ggplot(df_size_effect)+ geom_point(aes(x=n,y=cond_odds,color=label))+ geom_line(aes(x=n,y=cond_odds,color=label))   +facet_grid(cols=vars(term))+
  xlab("Hosehold size")+ylab("Conditional odds ratio")+ labs(color = "Time-point")


plot_ag_interactions<-function(int_model,int_range){
  matr_cats <- c("Child \n (0-17)","Young Adult \n (18-34) ","Older Adult \n (35-65)")
  matr_cn <- matrix(c(
    "N(1)~mm[age_cat=1,age_cat=1]", "N(1)~mm[age_cat=1,age_cat=2]",       "N(1)~mm[age_cat=1,age_cat=3]",
    NA,      "N(1)~mm[age_cat=2,age_cat=2]",       "N(1)~mm[age_cat=2,age_cat=3]",
    NA,      NA,       "N(1)~mm[age_cat=3,age_cat=3]"
  ), 3, 3, byrow=TRUE) %>% t() %>% `[`(,3:1)
  matr_coef <- coef(int_model)[c(matr_cn)] %>% matrix(3, 3, dimnames=list(matr_cats,rev(matr_cats)))
  
  matr_bl <- na.omit(unique(c(matr_cn)))
  matr_bl_rows <- lapply(matr_bl, function(bl) na.omit(row(matr_cn)[c(matr_cn==bl)]))
  matr_bl_cols <- lapply(matr_bl, function(bl) na.omit(col(matr_cn)[c(matr_cn==bl)]))
  bl_df <- mapply(function(r,c){
    c(xmin=min(r)-0.5, xmax=max(r)+0.5, ymin=min(c)-0.5, ymax=max(c)+0.5)
  }, matr_bl_rows, matr_bl_cols) %>% t() %>% as.data.frame()
  
  df<-as.data.frame(as.table(matr_coef))
  colnames(df)<-c("ag1","ag2","Coefficient")
  df$text_color<-ifelse(df$Coefficient>mean(range(int_range, na.rm=TRUE)),"up" ,"low" )
  df$coef_text <- format(df$Coefficient,digits=1,nsmall=1)
  gobj<-ggplot(df,aes(x=ag1,y=ag2,fill=Coefficient)) +geom_tile()+
    geom_text(aes(label=coef_text,color=text_color),show.legend = FALSE)+scale_colour_manual(values=c("cadetblue3","black"))+
    scale_fill_gradient(limits=c(floor(int_range[1]), ceiling(int_range[2])))+
    scale_x_discrete(position = "top",guide = guide_axis(angle = 90)) +
    xlab("") + ylab("")
  return(gobj)
}

plot_ag_interactions_homogeneous_mixing<-function(int_model,int_range){
  matr_cats <- c("Child \n (0-17)","Young Adult \n (18-34) ","Older Adult \n (35-65)")
  matr_cn <- matrix(c(
    "N(1)~mm[age_cat=1,age_cat=1]", "N(1)~mm[age_cat=1,age_cat=2]",       "N(1)~mm[age_cat=1,age_cat=3]",
    NA,      "N(1)~mm[age_cat=2,age_cat=2]",       "N(1)~mm[age_cat=2,age_cat=3]",
    NA,      NA,       "N(1)~mm[age_cat=3,age_cat=3]"
  ), 3, 3, byrow=TRUE) %>% t() %>% `[`(,3:1)
  matr_coef <- coef(int_model)[c(matr_cn)] %>% matrix(3, 3, dimnames=list(matr_cats,rev(matr_cats)))
  
  matr_bl <- na.omit(unique(c(matr_cn)))
  matr_bl_rows <- lapply(matr_bl, function(bl) na.omit(row(matr_cn)[c(matr_cn==bl)]))
  matr_bl_cols <- lapply(matr_bl, function(bl) na.omit(col(matr_cn)[c(matr_cn==bl)]))
  bl_df <- mapply(function(r,c){
    c(xmin=min(r)-0.5, xmax=max(r)+0.5, ymin=min(c)-0.5, ymax=max(c)+0.5)
  }, matr_bl_rows, matr_bl_cols) %>% t() %>% as.data.frame()
  
  df<-as.data.frame(as.table(matr_coef))
  colnames(df)<-c("ag1","ag2","Coefficient")
  print(df$Coefficient)
  df$Coefficient[!is.na(df$Coefficient)]<-mean(df$Coefficient,na.rm=TRUE)
  df$text_color<-ifelse(df$Coefficient>mean(range(int_range, na.rm=TRUE)),"up" ,"low" )
  df$coef_text <- format(df$Coefficient,digits=1,nsmall=1)
  gobj<-ggplot(df,aes(x=ag1,y=ag2,fill=Coefficient)) +geom_tile()+
    geom_text(aes(label=coef_text,color=text_color),show.legend = FALSE)+scale_colour_manual(values=c("cadetblue3","black"))+
    scale_fill_gradient(limits=c(floor(int_range[1]), ceiling(int_range[2])))+
    scale_x_discrete(position = "top",guide = guide_axis(angle = 90)) +
    xlab("") + ylab("")
  return(gobj)
}

matr_cn <- matrix(c(
  "N(1)~mm[age_cat=1,age_cat=1]", "N(1)~mm[age_cat=1,age_cat=2]",       "N(1)~mm[age_cat=1,age_cat=3]",
  NA,      "N(1)~mm[age_cat=2,age_cat=2]",       "N(1)~mm[age_cat=2,age_cat=3]",
  NA,      NA,       "N(1)~mm[age_cat=3,age_cat=3]"
), 3, 3, byrow=TRUE) %>% t() %>% `[`(,3:1)
rng = range(c(coef(model_onset_full)[c(matr_cn)],coef(model_yest_full)[c(matr_cn)],coef(model_final_full)[c(matr_cn)] ),na.rm = TRUE)



p1<-plot_ag_interactions(model_onset_full,rng)
p2<-plot_ag_interactions(model_yest_full,rng)
p3<-plot_ag_interactions(model_final_full,rng)

ggarrange(p1+ guides(fill="none"),
          p2+ guides(fill="none")+theme(axis.ticks.y = element_blank(),axis.text.y = element_blank()),
          p3+ theme(axis.ticks.y = element_blank(),axis.text.y = element_blank())+ theme(legend.position="right") + labs(fill = "Model coefficient"), 
          labels = c("Onset", "Enrollment", "Follow-up"),
          ncol = 3, nrow = 1,
          widths = c(1.1,1,1.3),
          heights =c(1,1,1)  ## Sono qui, adjust plot
)
ggsave(filename="Figure3.png",width=15,height = 4)


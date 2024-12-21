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

df_densities_onset_phys_onlypart<-readRDS("./data/densities_full_dataset_onset_phys_onlypart_symmetric.RDS")
df_densities_yest_phys_onlypart<-readRDS("./data/densities_full_dataset_yest_phys_onlypart_symmetric.RDS")
df_densities_final_phys_onlypart<-readRDS("./data/densities_full_dataset_final_phys_onlypart_symmetric.RDS")

df_densities_onset_phys_onlypart$dataset<-"onset"
df_densities_yest_phys_onlypart$dataset<-"yest"
df_densities_final_phys_onlypart$dataset<-"final"
all_hh_ids<-union(df_densities_onset_phys_onlypart$hh_id,df_densities_yest_phys_onlypart$hh_id)
all_hh_ids<-union(all_hh_ids,df_densities_final_phys_onlypart$hh_id)

densities_all_wide<-data.frame('hh_id'=all_hh_ids)


df_density_onlypart_all<-rbind(df_densities_onset_phys_onlypart,df_densities_yest_phys_onlypart,df_densities_final_phys_onlypart)
df_density_onlypart_all$dataset<-factor(df_density_onlypart_all$dataset,levels=c("onset","yest","final"))
df_density_onlypart_all$hh_size[df_density_onlypart_all$hh_size %in% c(5,6,7)]<-"5+"




fac_names <- c(
  `2` = "Household Size 2",
  `3` = "Household Size 3",
  `4` = "Household Size 4",
  `5+` = "Household Size 5\n and above"
)
png("Figure2.png",width=1200,height=800)
ggplot(df_density_onlypart_all)+geom_violin(scale = "count",aes(x=dataset,y=density), width=.2, position=position_dodge(0.05))+
  geom_violin(aes(x=dataset,y=density), width=.2, position=position_dodge(0.05))+
  facet_grid(rows=~hh_size,labeller = as_labeller(fac_names))+ 
  stat_summary(aes(x=dataset,y=density),fun = median, fun.min = median, fun.max = median,
               geom = "crossbar", width = 0.3, color = "red", size = 0.5)+
  scale_shape_manual("", values=c("Median"="x"))+ylab("Probability of a physical contact \n with another household member")+xlab("Time point")+
  scale_x_discrete(labels=c("onset" = "Onset", "yest" = "Enrollment","final" = "Follow-up"))+
  theme(axis.text.x = element_text( angle=90),text = element_text(size = 20))+ ggbeeswarm::geom_quasirandom(aes(x=dataset,y=density),alpha=0.1)
dev.off()


# divide them by household size

### Household size 2  -> not significant is commented
dum<-df_density_onlypart_all %>% filter(hh_size==2)
# onset vs yest
# a<-ks.test(dum$density[dum$dataset=="onset"],dum$density[dum$dataset=="yest"],alternative="two.sided")
# print(paste("Household size 2, onset vs yest, p-value=",a$p.value) )
# # onset vs final
# a<-ks.test(dum$density[dum$dataset=="onset"],dum$density[dum$dataset=="final"],alternative="two.sided")
# print(paste("Household size 2, onset vs final, p-value=",a$p.value) )
# # yest vs final
# a<-ks.test(dum$density[dum$dataset=="yest"],dum$density[dum$dataset=="final"],alternative="two.sided")
# print(paste("Household size 2, yest vs final, p-value=",a$p.value) )


### Household size 3  -> not significant is commented
dum<-df_density_onlypart_all %>% filter(hh_size==3)
# onset vs yest
# a<-ks.test(dum$density[dum$dataset=="onset"],dum$density[dum$dataset=="yest"],alternative="two.sided")
# print(paste("Household size 3, onset vs yest, p-value=",a$p.value) )
# onset vs final
# a<-ks.test(dum$density[dum$dataset=="onset"],dum$density[dum$dataset=="final"],alternative="two.sided")
# print(paste("Household size 3, onset vs final, p-value=",a$p.value) )
# yest vs final
# a<-ks.test(dum$density[dum$dataset=="yest"],dum$density[dum$dataset=="final"],alternative="two.sided")
# print(paste("Household size 3, yest vs final, p-value=",a$p.value) )

### Household size 4   -> not significant is commented
dum<-df_density_onlypart_all %>% filter(hh_size==4)
# onset vs yest
a<-ks.test(dum$density[dum$dataset=="onset"],dum$density[dum$dataset=="yest"],alternative="two.sided")
print(paste("Household size 4, onset vs yest, p-value=",a$p.value) )
# onset vs final
# a<-ks.test(dum$density[dum$dataset=="onset"],dum$density[dum$dataset=="final"],alternative="two.sided")
# print(paste("Household size 4, onset vs final, p-value=",a$p.value) )
# yest vs final
a<-ks.test(dum$density[dum$dataset=="yest"],dum$density[dum$dataset=="final"],alternative="two.sided")
print(paste("Household size 4, yest vs final, p-value=",a$p.value) )


### Household size 5+
dum<-df_density_onlypart_all %>% filter(hh_size=="5+")
# onset vs yest
a<-ks.test(dum$density[dum$dataset=="onset"],dum$density[dum$dataset=="yest"],alternative="two.sided")
print(paste("Household size 5+, onset vs yest, p-value=",a$p.value) )
# onset vs final
# a<-ks.test(dum$density[dum$dataset=="onset"],dum$density[dum$dataset=="final"],alternative="two.sided")
# print(paste("Household size 5+, onset vs final, p-value=",a$p.value) )
# yest vs final
# a<-ks.test(dum$density[dum$dataset=="yest"],dum$density[dum$dataset=="final"],alternative="two.sided")
# print(paste("Household size 5+, yest vs final, p-value=",a$p.value) )







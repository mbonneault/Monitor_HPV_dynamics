rm(list = ls())
library(ggplot2)

#setwd("Repositeries")
JeuComplet<-read.csv2("JeuComplet_tot_Prev_NV_V_part.csv")
JeuComplet<-cbind(JeuComplet,dif_vac_unvac=JeuComplet$val_prevNV_Vac-JeuComplet$val_prevNV_unVac)
JeuComplet<-cbind(JeuComplet,dif_AvAp=JeuComplet$val_prevNVtot-JeuComplet$val_prevNV_5)

JeuComplet0_65<-JeuComplet[which(JeuComplet$year<25 & JeuComplet$partner=='all'),]

table<-data.frame()
taille<-c()
for(i in 1:length(unique(JeuComplet0_65$year))){

  taille<-c(taille,mean(JeuComplet0_65$Ntot[which(JeuComplet0_65$year==unique(JeuComplet0_65$year)[i])])/mean(JeuComplet0_65$Ntot[which(JeuComplet0_65$year==20)])*100)
}
table<-data.frame(year=as.factor(unique(JeuComplet0_65$year)),age=c("15","15-17","15-19","15-24","15-29","15-29"),Proportion=taille)

taille=14
ggplot(data=table, aes(x=year, y=Proportion)) +
  geom_bar(stat="identity", fill="steelblue")+labs(x = "Years after vaccine introduction", y = "Women concerned by the vaccination (%)")+
  geom_text(aes(label=age), vjust=-0.3, size=4)+theme_minimal()+theme(legend.title=element_text(size=taille),legend.text=element_text(size=taille-1),axis.text.y = element_text(size =taille),axis.text.x = element_text(size =taille),axis.title.x = element_text(size = taille),axis.title.y = element_text(size =taille),legend.position = "none",title=element_text(size=taille))


                                                              
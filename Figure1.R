rm(list = ls())

library(ggplot2)
library("gridExtra")
setwd("C:/Users/mel/Documents/HPV/modele/Analyse_A2/02092020/donnees")

JeuComplet<-read.csv2("JeuComplet_tot_Prev_NV_V_tot.csv")

#### Data formatting ---- All women

Jeu<-JeuComplet[which(JeuComplet$vac==0.60 & (JeuComplet$interaction ==0.5 |JeuComplet$interaction ==1 |JeuComplet$interaction ==1.5) ),]
JeuCompletbis<-data.frame()

#Tested values for interaction vaccine coverage, year after vaccine introduction and studing all population 

JeuP2=data.frame()

interac=c(0.5,1.0,1.5)
vac=c(0.6)
year=c(-10,-5,1,3,5,10,15,20,25,30,35,40)
o=1

for(n in 1:length(interac))
{
  for(m in 1:length(vac))
  {
    for(v in 1:length(year)){
      JeuP2 = rbind(JeuP2,data.frame(interac=interac[n],vac=vac[m],Year=year[v],ref=o))
    }
    o=o+1
  }
}

# calculation of median values and empirical interval
for(int in 1:length((JeuP2$ref))){
  Jeubis<-Jeu[which(Jeu$interaction==JeuP2$interac[int] & Jeu$vac==JeuP2$vac[int] &Jeu$year==JeuP2$Year[int]),]
  
  JeuCompletbis<-rbind(JeuCompletbis,data.frame(Interaction=JeuP2$interac[int],vac=JeuP2$vac[int],year=JeuP2$Year[int],prevNVtotM=median(Jeubis$val_prevNVtot)*100,prevNV_VacM=median(Jeubis$val_prevNV_Vac)*100,prevNV_unVacM=median(Jeubis$val_prevNV_unVac)*100,prevVtotM=median(Jeubis$val_prevVtot)*100,prevV_VacM=median(Jeubis$val_prevV_Vac)*100,prevV_unVacM=median(Jeubis$val_prevV_unVac)*100))
  
}


JeuCompletbis$Interaction<-as.factor(JeuCompletbis$Interaction)

taille=14
# figure 1A
all<-ggplot(JeuCompletbis) +
  geom_line( aes(year,prevNVtotM, color=Interaction,linetype="NV"),size=1.5)+theme_bw()+geom_line(aes(year,prevVtotM, color=Interaction, linetype="V"),size=1.5)+labs(title= "A) All women")+ xlab("Years before and after vaccine introduction") +   ylab(" V- and NV- genotype prevalence (%)")+ geom_vline(xintercept=c(0,15), linetype="dotted")+ scale_linetype_manual("Prevalence", values=c("NV"=1,"V"=2),labels=c("NV","V"))+scale_color_manual("Interaction",values=c('lightcoral','black','mediumseagreen'))+ylim(0,40)+theme(legend.title=element_text(size=taille),legend.text=element_text(size=taille-1),axis.text.y = element_text(size =taille),axis.text.x = element_text(size =taille),axis.title.x = element_text(size = taille),axis.title.y = element_text(size =taille),legend.position = "none",title=element_text(size=taille))+xlim(-5,40)





JeuVac<-read.csv2("C:/Users/mel/Documents/HPV/modele/Analyse_A2/03072020/donnees/JeuComplet_tot_Prev_NV_V_part.csv")


#### Data formatting ---- vaccinated and unvaccinated women
JeuVac<-JeuVac[which(JeuVac$vac==0.60 & (JeuVac$interaction ==0.5 |JeuVac$interaction ==1 |JeuVac$interaction ==1.5) ),]
JeuCompletVac<-data.frame()

JeuP=data.frame()
#Tested values for interaction vaccine coverage, year after vaccine introduction and studing all population 
interac=c(0.5,1.0,1.5)
vac=c(0.6)
year=c(1,3,5,10,15,20,25,30,35,40)
o=1

for(n in 1:length(interac))
{
  for(m in 1:length(vac))
  {
    for(v in 1:length(year)){
      JeuP = rbind(JeuP,data.frame(interac=interac[n],vac=vac[m],Year=year[v],ref=o))
    }
    o=o+1
  }
}

# calculation of median values and empirical interval
for(int in 1:length((JeuP$ref))){
  Jeubis<-JeuVac[which(JeuVac$interaction==JeuP$interac[int] & JeuVac$vac==JeuP$vac[int] &JeuVac$year==JeuP$Year[int]),]
  
  JeuCompletVac<-rbind(JeuCompletVac,data.frame(Interaction=JeuP$interac[int],vac=JeuP$vac[int],year=JeuP$Year[int],prevNVtotM=median(Jeubis$val_prevNVtot)*100,prevNV_VacM=median(Jeubis$val_prevNV_Vac)*100,prevNV_unVacM=median(Jeubis$val_prevNV_unVac)*100,prevVtotM=median(Jeubis$val_prevVtot)*100,prevV_VacM=median(Jeubis$val_prevV_Vac)*100,prevV_unVacM=median(Jeubis$val_prevV_unVac)*100))
  
}


JeuCompletVac$Interaction<-as.factor(JeuCompletVac$Interaction)

JeuCompletVacC<-JeuCompletVac[which(JeuCompletVac$year>0),]

# figure 1B
Vac<-ggplot(JeuCompletVacC) +
  geom_line( aes(year,prevNV_VacM, color=Interaction,linetype="NV"),size=1.5)+theme_bw()+geom_line(aes(year,prevV_VacM, color=Interaction, linetype="V"),size=1.5)+labs(title= "B) Vaccinated women")+ xlab("Years after vaccine introduction") +   ylab(" V- and NV- genotype prevalence (%)")+ geom_vline(xintercept=c(0,15), linetype="dotted")+ scale_linetype_manual("Prevalence", values=c("NV"=1,"V"=2),labels=c("NV","V"))+scale_color_manual("Interaction",values=c('lightcoral','black','mediumseagreen'))+ylim(0,40)+xlim(-5,40)+theme(legend.title=element_text(size=taille),legend.text=element_text(size=taille-1),axis.text.y = element_text(size =taille),axis.text.x = element_text(size =taille),axis.title.x = element_text(size = taille),axis.title.y = element_text(size =taille),legend.position = "none",title=element_text(size=taille))

# figure 1C
unVac<-ggplot(JeuCompletVacC) +
  geom_line( aes(year,prevNV_unVacM, color=Interaction,linetype="NV"),size=1.5)+theme_bw()+geom_line(aes(year,prevV_unVacM, color=Interaction, linetype="V"),size=1.5)+labs(title= "C) Unvaccinated women")+ xlab("Years after vaccine introduction") +   ylab(" V- and NV- genotype prevalence (%)")+ geom_vline(xintercept=c(0,15), linetype="dotted")+ scale_linetype_manual("Prevalence", values=c("NV"=1,"V"=2),labels=c("NV","V"))+scale_color_manual("Interaction",values=c('lightcoral','black','mediumseagreen'))+ylim(0,40)+xlim(-5,40)+theme(legend.title=element_text(size=taille),legend.text=element_text(size=taille-1),axis.text.y = element_text(size =taille),axis.text.x = element_text(size =taille),axis.title.x = element_text(size = taille),axis.title.y = element_text(size =taille),legend.position = "none",title=element_text(size=taille))


grid.arrange(all,Vac,unVac, nrow=3)

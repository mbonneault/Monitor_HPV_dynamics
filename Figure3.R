rm(list = ls())

library(ggplot2)
library(viridis)

#setwd("Repositeries")
Revue<-read.csv2("review.csv")

#### Data formatting ----
Revue<-cbind(Revue,taille=Revue$taille_ref+Revue$taille_compa,tailleCat=Revue$taille_re,nameCat=c(rep("<500",25),rep("500-999",25),rep("1000- 4 999",25),rep("4 999- 9 999",25),rep(">9 999",18)))
Revue<-cbind(Revue,Nametaille=paste(Revue$taille,Revue$Name,sep=" "))
Revue$tailleCat[which(Revue$taille<500)]<-1
Revue$tailleCat[which(Revue$taille<1000 & Revue$taille>=500)]<-2
Revue$tailleCat[which(Revue$taille<5000 & Revue$taille>=1000)]<-3
Revue$tailleCat[which(Revue$taille<10000 & Revue$taille>=5000)]<-4
Revue$tailleCat[which(Revue$taille>=10000)]<-5

Revue$nameCat[which(Revue$taille<500)]<-"<500"
Revue$nameCat[which(Revue$taille<1000 & Revue$taille>=500)]="500-999"
Revue$nameCat[which(Revue$taille<5000 & Revue$taille>=1000)]<-"1000- 4 999"
Revue$nameCat[which(Revue$taille<10000 & Revue$taille>=5000)]<-"4 999- 9 999"
Revue$nameCat[which(Revue$taille>=10000)]<-">9 999"
Revue$vac<-as.numeric(as.character(Revue$vac))
Revue$decal<-as.numeric(as.character(Revue$decal))
Revue$Difference<-as.numeric(as.character(Revue$Difference))
Revue$IC_min<-as.numeric(as.character(Revue$IC_min))
Revue$IC_max<-as.numeric(as.character(Revue$IC_max))

Revue$nameCat <- factor(Revue$nameCat , levels = c("<500","500-999","1000- 4 999","4 999- 9 999",">9 999"))

revue_prePost<-Revue[which(Revue$design=="Pre/post" & Revue$genotype=="NV" ),]
revue_prePost<-as.data.frame(revue_prePost)

#figure 3C
ggplot()+ scale_y_continuous(limits=c(-26,28),breaks=seq(-30,30,10))+geom_hline(yintercept=c(0),color='red', linetype="dashed")+geom_segment(data=revue_prePost,aes(x=anne_min,xend=anne_max,y=Difference,yend=Difference, color=vac))+scale_color_viridis(name="Vaccine \n  coverage(%)",option = "D", direction = -1)+ theme_bw()+labs(title="Post-vs-prevaccine comparison - NV genotypes", x = "Years after vaccine introduction", y = "NV prevalence difference (%)")+geom_errorbar(data=revue_prePost,aes(x=decal,ymin=IC_min,ymax=IC_max, color=vac), width=.3)+geom_point(data=revue_prePost,aes(x=decal,y=Difference, color=vac,size=nameCat))+ scale_x_continuous(limits=c(0,10),breaks=seq(0,10,2))+scale_size_manual(name="Sample size",values=c(1,2,3,4,5),labels=c("<500","500-999","1000- 4 999","4 999- 9 999",">9 999"))


revue_prePost<-Revue[which(Revue$design=="Pre/post" & Revue$genotype=="V" ),]
revue_prePost<-as.data.frame(revue_prePost)
#figure 3A
ggplot()+ scale_y_continuous(limits=c(-30,6.5),breaks=seq(-30,5,5))+geom_hline(yintercept=c(0),color='red', linetype="dashed")+geom_segment(data=revue_prePost,aes(x=anne_min,xend=anne_max,y=Difference,yend=Difference, color=vac))+scale_color_viridis(name="Vaccine \n  coverage (%)",option = "D", direction = -1)+ theme_bw()+labs(title="Post-vs-prevaccine comparison - V genotypes", x = "Years after vaccine introduction", y = "V prevalence difference (%)")+geom_errorbar(data=revue_prePost,aes(x=decal,ymin=IC_min,ymax=IC_max, color=vac), width=.3)+geom_point(data=revue_prePost,aes(x=decal,y=Difference, color=vac,size=nameCat))+ scale_x_continuous(limits=c(0,10),breaks=seq(0,10,2))+scale_size_manual(name="Sample size",values=c(1,2,3,4,5),labels=c("<500","500-999","1000- 4 999","4 999- 9 999",">9 999"))



revue_prePost<-Revue[which(Revue$design=="Vac/Unvac" & Revue$genotype=="NV" ),]
revue_prePost<-as.data.frame(revue_prePost)

#figure 3D
ggplot()+ scale_y_continuous(limits=c(-26,28),breaks=seq(-30,30,10))+geom_hline(yintercept=c(0),color='red', linetype="dashed")+geom_segment(data=revue_prePost,aes(x=anne_min,xend=anne_max,y=Difference,yend=Difference, color=vac))+scale_color_viridis(name="Vaccine \n  coverage (%)",option = "D", direction = -1)+ theme_bw()+labs(title="Vac-vs-Unvac comparison - NV genotypes", x = "Years after vaccine introduction", y = "NV prevalence difference (%)")+geom_errorbar(data=revue_prePost,aes(x=decal,ymin=IC_min,ymax=IC_max, color=vac), width=.3)+geom_point(data=revue_prePost,aes(x=decal,y=Difference, color=vac,size=nameCat))+ scale_x_continuous(limits=c(0,10),breaks=seq(0,10,2))+scale_size_manual(name="Sample size",values=c(1,2,3,4,5),labels=c("<500","500-999","1000- 4 999","4 999- 9 999",">9 999"))


revue_prePost<-Revue[which(Revue$design=="Vac/Unvac" & Revue$genotype=="V" ),]
revue_prePost<-as.data.frame(revue_prePost)

#to see all points
revue_prePost$decal[1]<-2.55 
revue_prePost$decal[6]<-2.6

#figure 3B
ggplot()+ scale_y_continuous(limits=c(-30,6.5),breaks=seq(-30,5,5))+geom_hline(yintercept=c(0),color='red', linetype="dashed")+geom_segment(data=revue_prePost,aes(x=anne_min,xend=anne_max,y=Difference,yend=Difference, color=vac))+scale_color_viridis(name="Vaccine \n  coverage (%)",option = "D", direction = -1)+ theme_bw()+labs(title="Vac-vs-Unvac comparison - V genotypes", x = "Years after vaccine introduction", y = "V prevalence difference (%)")+geom_errorbar(data=revue_prePost,aes(x=decal,ymin=IC_min,ymax=IC_max, color=vac), width=.3)+geom_point(data=revue_prePost,aes(x=decal,y=Difference, color=vac,size=nameCat))+ scale_x_continuous(limits=c(0,10),breaks=seq(0,10,2))+scale_size_manual(name="Sample size",values=c(1,2,3,4,5),labels=c("<500","500-999","1000- 4 999","4 999- 9 999",">9 999"))


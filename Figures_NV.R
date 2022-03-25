rm(list = ls())


library(ggplot2)
library(ggnewscale)


#setwd("Repositeries")
JeuComplet<-read.csv2("JeuComplet_tot_Prev_NV_V_part.csv")

#### Data formatting ----
JeuComplet<-cbind(JeuComplet,dif_vac_unvac=JeuComplet$val_prevNV_Vac-JeuComplet$val_prevNV_unVac)
JeuComplet<-cbind(JeuComplet,dif_AvAp=JeuComplet$val_prevNVtot-JeuComplet$val_prevNV_5)


#Tested values for interaction vaccine coverage, year after vaccine introduction and studing all population or population according the number of partner(s) during the year( 1 to 3 or more)


JeuP2<-data.frame()

interaction=c(0.4,0.5,0.6,0.7,0.8,0.9,1.0,1.1,1.2,1.3,1.4,1.5)
vac=as.numeric(seq(0.2,0.9,0.1))
kvac=c(4,2.3,1.5,1,0.667,0.43,0.25,0.11)
year=c(1,3,seq(5,40,5))
partner=c('all','1-3','>3')
o=1

for(n in 1:length(interaction)){
  for(m in 1:length(vac)){
    for(v in 1:length(year)){
      for(l in 1:length(partner)){
        JeuP2 = rbind(JeuP2,data.frame(interaction=interaction[n],vac=vac[m],year=year[v],partner=partner[l],ref=o))
      }
      o=o+1
    }
  }
}

#calculation of the number of sample size required to show prevalences differences in Pre-vs-post(nsujetAvAptot), Prevaccine sample size required nsujetAvAp_min and Vaccinated-vs-unvaccinated (nsujetvacUnvac)
vac<-as.numeric(as.character(vac))
JeuP2$vac<-as.numeric(as.character(JeuP2$vac))
nsujetvac=c()
nsujetunvac=c()
nsujetvactot=c()
nsujetAvAptot=c()

nsujetAvAp_min=c()
JeuComplet$vac<-as.numeric(as.character(JeuComplet$vac))
for(int in 1:length(JeuComplet$vac)){
  pos=which(vac==as.numeric(as.character(JeuComplet$vac[int])))
  if(is.na(pos)){
    print(int)
  }
  nsujetvac=c(nsujetvac,(kvac[pos]+1)/kvac[pos]*2.802^2/(2*(asin(sqrt(JeuComplet$val_prevNV_Vac[int]))-asin(sqrt(JeuComplet$val_prevNV_unVac[int]))))^2)
  nsujetunvac=c(nsujetunvac,kvac[pos]*nsujetvac[int])
  nsujetvactot=c(nsujetvactot,4*2.802^2/(2*(asin(sqrt(JeuComplet$val_prevNV_Vac[int]))-asin(sqrt(JeuComplet$val_prevNV_unVac[int]))))^2)
  nsujetAvAptot=c(nsujetAvAptot,4*2.802^2/(2*(asin(sqrt(JeuComplet$val_prevNVtot[int]))-asin(sqrt(JeuComplet$val_prevNV_5[int]))))^2)
  nsujetAvAp_min=c(nsujetAvAp_min,2.802^2/4/(asin(sqrt(JeuComplet$val_prevNVtot[int]))-asin(sqrt(JeuComplet$val_prevNV_5[int])))^2)
}

JeuComplet<-cbind(JeuComplet,nsujetvac=nsujetvac,nsujetunvac=nsujetunvac,nsujetvactot=nsujetvactot,nsujetvacUnvac=nsujetvac+nsujetunvac,nsujetAvAptot=nsujetAvAptot,nsujetAvAp_min=nsujetAvAp_min)

JeuCompletMIC<-data.frame()

# calculation of median values and empirical interval
for(int in 1:length((JeuP2$ref))){
  Jeubis<-JeuComplet[which(JeuComplet$interaction==JeuP2$interaction[int] & JeuComplet$vac==JeuP2$vac[int] & JeuComplet$year==JeuP2$year[int] & JeuComplet$partner==JeuP2$partner[int]),]
  
  JeuCompletMIC<-rbind(JeuCompletMIC,data.frame(Interaction=JeuP2$interaction[int],vac=JeuP2$vac[int],year=JeuP2$year[int],partner=as.character(JeuP2$partner[int]),dif_vac_unvacM=median(Jeubis$dif_vac_unvac),dif_vac_unvacIm=Jeubis$dif_vac_unvac[order(Jeubis$dif_vac_unvac)[5]],dif_vac_unvacIp=Jeubis$dif_vac_unvac[order(Jeubis$dif_vac_unvac)[95]],dif_AvApM=median(Jeubis$dif_AvAp),dif_AvApIm=Jeubis$dif_AvAp[order(Jeubis$dif_AvAp)[5]],dif_AvApIp=Jeubis$dif_AvAp[order(Jeubis$dif_AvAp)[95]],
dif_vac_unvacM_abs=median(abs(Jeubis$dif_vac_unvac)),dif_vac_unvacIm_abs=abs(Jeubis$dif_vac_unvac[order(abs(Jeubis$dif_vac_unvac))[5]]),dif_vac_unvacIp_abs=abs(Jeubis$dif_vac_unvac[order(abs(Jeubis$dif_vac_unvac))[95]]),dif_AvApM_abs=median(abs(Jeubis$dif_AvAp)),dif_AvApIm_abs=abs(Jeubis$dif_AvAp[order(abs(Jeubis$dif_AvAp))[5]]),dif_AvApIp_abs=abs(Jeubis$dif_AvAp[order(abs(Jeubis$dif_AvAp))[95]]), nsujetvacM=median(Jeubis$nsujetvac),nsujetvacIm=Jeubis$nsujetvac[order(Jeubis$nsujetvac)[5]],nsujetvacIp=Jeubis$nsujetvac[order(Jeubis$nsujetvac)[95]],nsujetunvacM=median(Jeubis$nsujetunvac), nsujetunvacIm=Jeubis$nsujetunvac[order(Jeubis$nsujetunvac)[5]],nsujetunvacIp=Jeubis$nsujetunvac[order(Jeubis$nsujetunvac)[95]],nsujetvactotM=median(Jeubis$nsujetvactot),nsujetvactotIm=Jeubis$nsujetvactot[order(Jeubis$nsujetvactot)[5]],nsujetvactotIp=Jeubis$nsujetvactot[order(Jeubis$nsujetvactot)[95]],nsujetvacUnvacM=median(Jeubis$nsujetvacUnvac),nsujetvacUnvacIm=Jeubis$nsujetvacUnvac[order(Jeubis$nsujetvacUnvac)[5]],nsujetvacUnvacIp=Jeubis$nsujetvacUnvac[order(Jeubis$nsujetvacUnvac)[95]],sujetAvAptotM=median(Jeubis$nsujetAvAptot),nsujetAvAptotIm=Jeubis$nsujetAvAptot[order(Jeubis$nsujetAvAptot)[5]],nsujetAvAptotIp=Jeubis$nsujetAvAptot[order(Jeubis$nsujetAvAptot)[95]],nsujetAvAp_minM=median(Jeubis$nsujetAvAp_min),nsujetAvAp_minIm=Jeubis$nsujetAvAp_min[order(Jeubis$nsujetAvAp_min)[5]],nsujetAvAp_minIp=Jeubis$nsujetAvAp_min[order(Jeubis$nsujetAvAp_min)[95]]))
}


#################################################################################
## Figures plot ---- ###############################
################################################################################
paletteFunc <- colorRampPalette(c('darkmagenta','mediumvioletred','pink'),alpha=T)
palette1 <- paletteFunc(4)
barplot(seq(1,4,1), col=palette1)

paletteFunc <- colorRampPalette(c('mediumturquoise','dodgerblue','darkblue'),alpha=T)
palette2 <- paletteFunc(4)
barplot(seq(1,4,1), col=palette2) 

#selection of interaction values and vaccine coverage
JeuCompletMIC<-JeuCompletMIC[which(JeuCompletMIC$Interaction>0.4 & (JeuCompletMIC$vac==0.2 |JeuCompletMIC$vac==0.4| JeuCompletMIC$vac==0.6|JeuCompletMIC$vac==0.8)),]


JeuCompletT10<-JeuCompletMIC[which((JeuCompletMIC$Interaction==0.5 |JeuCompletMIC$Interaction==0.9|JeuCompletMIC$Interaction==1 |JeuCompletMIC$Interaction==1.5 |JeuCompletMIC$Interaction==1.1) & (JeuCompletMIC$partner=="all")),]
JeuCompletT10$Interaction<-as.factor(JeuCompletT10$Interaction)

JeuCompletT10<-cbind(JeuCompletT10,Interac=c(rep("Strong interaction (0.5)",40),rep("Weak interaction (0.9)",40),rep("No interaction (1)",40),rep("Weak interaction (1.1)",40),rep("Strong interaction (1.5)",40)))
JeuCompletT10$Interac<-as.factor(JeuCompletT10$Interac)
JeuCompletT10$Interac<-factor(JeuCompletT10$Interac,level=c("Strong interaction (0.5)","Weak interaction (0.9)","No interaction (1)","Weak interaction (1.1)","Strong interaction (1.5)"))

JeuCompletT10$vac<-as.factor(JeuCompletT10$vac*100)
JeuCompletT10$dif_vac_unvacM_abs<-JeuCompletT10$dif_vac_unvacM_abs*100
JeuCompletT10$dif_vac_unvacIm_abs<-JeuCompletT10$dif_vac_unvacIm_abs*100
JeuCompletT10$dif_vac_unvacIp_abs<-JeuCompletT10$dif_vac_unvacIp_abs*100
JeuCompletT10$dif_AvApM_abs<-JeuCompletT10$dif_AvApM_abs*100
JeuCompletT10$dif_AvApIm_abs<-JeuCompletT10$dif_AvApIm_abs*100
JeuCompletT10$dif_AvApIp_abs<-JeuCompletT10$dif_AvApIp_abs*100

JeuCompletT10$dif_vac_unvacM<-JeuCompletT10$dif_vac_unvacM*100
JeuCompletT10$dif_vac_unvacIm<-JeuCompletT10$dif_vac_unvacIm*100
JeuCompletT10$dif_vac_unvacIp<-JeuCompletT10$dif_vac_unvacIp*100
JeuCompletT10$dif_AvApM<-JeuCompletT10$dif_AvApM*100
JeuCompletT10$dif_AvApIm<-JeuCompletT10$dif_AvApIm*100
JeuCompletT10$dif_AvApIp<-JeuCompletT10$dif_AvApIp*100

#figure 5 - Difference of prevalence 
ggplot(JeuCompletT10)+geom_hline(yintercept=c(0),color='red', linetype="dashed")+geom_vline(xintercept=15,color='grey',size = 0.2,linetype="dashed")+geom_point(aes(year,dif_vac_unvacM_abs,color=vac,shape=vac))+geom_line(aes(year,dif_vac_unvacM_abs,color=vac))+geom_hline(yintercept=seq(2,8,2),color='grey',size = 0.2)+  theme_classic() +labs(x = "Years after vaccine introduction", y = "Absolute NV-prevalence difference (%)") +geom_errorbar(aes(year,ymin = dif_vac_unvacIm_abs, ymax = dif_vac_unvacIp_abs, color=vac),width = 0.1)+ scale_color_manual(name="Vaccine  \n coverage (%) \n \n Vac-vs-Unvac",values=palette2,guide=guide_legend(override.aes=list(shape=c(1,2,0,6)),order=1))+ labs(title= "Competitive interactions      Neutral interaction                Synergistic interactions")+ scale_y_continuous(limits=c(-0.0, 9.1),breaks=seq(0,8,2))+scale_shape_manual(values=c(1,2,0,6),guide=F)+ new_scale_color() +new_scale("shape")+
  geom_point(aes(year,dif_AvApM_abs,color=vac,shape=vac,fill=vac))+geom_line(aes(year,dif_AvApM_abs,color=vac))+ geom_errorbar(aes(year,ymin = dif_AvApIm_abs, ymax = dif_AvApIp_abs, color=vac),width = 0.1)+scale_shape_manual(values=c(21,24,22,25),guide=FALSE)+ scale_color_manual(name="Post-vs-Pre",values=palette1[4:1],guide=guide_legend(order=2))+scale_fill_manual(values=palette1[4:1],guide=FALSE)+
  guides(color=guide_legend(override.aes=list(shape=c(21,24,22,25),fill=palette1[4:1])))+ facet_grid(. ~ Interac)+
  theme_bw()


#figure 5 - Sample size required
ggplot(JeuCompletT10) + 
  geom_hline(yintercept=c(10^3,10^5,10^7,10^9),color='grey',size = 0.2)+geom_vline(xintercept=15,color='grey',size = 0.2,linetype="dashed")+geom_point(aes(year,nsujetvacUnvacM,color=vac,shape=vac))+geom_line(aes(year,nsujetvacUnvacM,color=vac))+theme_classic() +labs(x = "Years after vaccine introduction", y = "Sample size required") +geom_errorbar(aes(year,ymin = nsujetvacUnvacIm, ymax = nsujetvacUnvacIp, color=vac),width = 0.02)+ scale_color_manual(name="Vaccine  \n coverage(%) \n \n Vac-vs-Unvac",values=palette2,guide=guide_legend(override.aes=list(shape=c(1,2,0,6)),order=1))+labs(title= "Competitive interactions               Neutral interaction          Synergistic interactions")+ scale_y_log10(limits=c(5*10^2,5*10^9),breaks=c(10^3,10^5,10^7,10^9))+scale_shape_manual(values=c(1,2,0,6),guide=F)+ new_scale_color() +new_scale("shape")+
  geom_point( aes(year,sujetAvAptotM,color=vac,shape=vac,fill=vac))+geom_line( aes(year,sujetAvAptotM,color=vac))+ geom_errorbar(aes(year,ymin = nsujetAvAptotIm, ymax = nsujetAvAptotIp, color=vac),width = 0.02)+scale_color_manual(name="Post-vs-Pre",values=palette1[4:1],guide=guide_legend(order=2))+scale_shape_manual(values=c(21,24,22,25),guide=FALSE)+scale_fill_manual(values=palette1[4:1],guide=FALSE)+
  guides(color=guide_legend(override.aes=list(shape=c(21,24,22,25),fill=palette1[4:1])))+ facet_grid(. ~ Interac)+
  theme_bw()


#Figure S2 - Prevaccine sample size required
ggplot(JeuCompletT10) + 
  geom_hline(yintercept=c(10^3,10^5,10^7,10^9),color='grey',size = 0.2)+geom_vline(xintercept=15,color='grey',size = 0.2,linetype="dashed")+ labs(title= "Competitive interactions               Neutral interaction          Synergistic interactions")+labs(x = "Years after vaccine introduction", y = "Prevaccine sample size required")+ scale_y_log10(limits=c(10^2,1.1*10^9),breaks=c(10^3,10^5,10^7,10^9))+  geom_point( aes(year,nsujetAvAp_minM,color=vac,shape=vac,fill=vac))+geom_line( aes(year,nsujetAvAp_minM,color=vac))+ geom_errorbar(aes(year,ymin = nsujetAvAp_minIm, ymax = nsujetAvAp_minIp, color=vac),width = 0.02)+scale_color_manual(name="Vaccine  \n coverage(%) \n \n Post-vs-Pre",values=palette1[4:1],guide=guide_legend(order=2))+scale_shape_manual(values=c(21,24,22,25),guide=FALSE)+scale_fill_manual(values=palette1[4:1],guide=FALSE)+
  guides(color=guide_legend(override.aes=list(shape=c(21,24,22,25),fill=palette1[4:1])))+ facet_grid(. ~ Interac)+
  theme_bw()



#selection of one synergistic interaction value and population according sexual activity 
JeuCompletT10<-JeuCompletMIC[which(JeuCompletMIC$Interaction==1.5 & !(JeuCompletMIC$partner=="all")),]

JeuCompletT10$Interaction<-as.factor(JeuCompletT10$Interaction)
JeuCompletT10<-cbind(JeuCompletT10,Part=rep(c("1-3 partners",">3 partners"),40))
JeuCompletT10$Part<-factor(JeuCompletT10$Part,levels=c("1-3 partners",">3 partners"))



JeuCompletT10$vac<-as.factor(JeuCompletT10$vac*100)


JeuCompletT10$dif_vac_unvacM<-JeuCompletT10$dif_vac_unvacM*100
JeuCompletT10$dif_vac_unvacIm<-JeuCompletT10$dif_vac_unvacIm*100
JeuCompletT10$dif_vac_unvacIp<-JeuCompletT10$dif_vac_unvacIp*100
JeuCompletT10$dif_AvApM<-JeuCompletT10$dif_AvApM*100
JeuCompletT10$dif_AvApIm<-JeuCompletT10$dif_AvApIm*100
JeuCompletT10$dif_AvApIp<-JeuCompletT10$dif_AvApIp*100

#figure 6 B Nv genotypes
ggplot(JeuCompletT10)+geom_point(aes(year,dif_vac_unvacM,color=vac,shape=vac))+geom_vline(xintercept=15,color='grey',size = 0.2,linetype="dashed")+geom_line(aes(year,dif_vac_unvacM,color=vac))+geom_hline(yintercept=c(seq(-100,10,10)),color='grey',size = 0.2)+    theme_bw() +labs(x = "Years after vaccine introduction", y =  "Prevalence of difference for NV (%)") +geom_errorbar(aes(year,ymin = dif_vac_unvacIm, ymax = dif_vac_unvacIp, color=vac),width = 0.1)+ scale_color_manual(name="Vaccine  \n coverage (%) \n \n Vac-Unvac",values=palette2,guide=guide_legend(override.aes=list(shape=c(1,2,0,6)),order=1))+ labs(title= "B. NV genotypes ")+ scale_y_continuous(limits=c(-100,11),breaks=seq(-100,10,20))+geom_hline(yintercept=c(0),color='red', linetype="dashed")+scale_shape_manual(values=c(1,2,0,6),guide=F)+ new_scale_color() +new_scale("shape")+
  geom_point(aes(year,dif_AvApM,color=vac,shape=vac,fill=vac))+geom_line(aes(year,dif_AvApM,color=vac))+ geom_errorbar(aes(year,ymin = dif_AvApIm, ymax = dif_AvApIp, color=vac),width = 0.1)+scale_shape_manual(values=c(21,24,22,25),guide=FALSE)+ scale_color_manual(name="Post-Pre",values=palette1[4:1],guide=guide_legend(order=2))+scale_fill_manual(values=palette1[4:1],guide=FALSE)+
  guides(color=guide_legend(override.aes=list(shape=c(21,24,22,25),fill=palette1[4:1])))+ facet_grid(. ~ Part)

#Figure S4 -5 years after vaccine introduction

JeuCompletT10<-JeuCompletMIC[which(JeuCompletMIC$year==5 & JeuCompletMIC$partner=="all"),]
JeuCompletT10<-cbind(JeuCompletT10,essai=JeuCompletT10$Interaction+JeuCompletT10$vac*0.1-0.05)
JeuCompletT10$vac<-as.factor(JeuCompletT10$vac*100)

JeuCompletT10$dif_vac_unvacM_abs<-JeuCompletT10$dif_vac_unvacM_abs*100
JeuCompletT10$dif_vac_unvacIm_abs<-JeuCompletT10$dif_vac_unvacIm_abs*100
JeuCompletT10$dif_vac_unvacIp_abs<-JeuCompletT10$dif_vac_unvacIp_abs*100
JeuCompletT10$dif_AvApM_abs<-JeuCompletT10$dif_AvApM_abs*100
JeuCompletT10$dif_AvApIm_abs<-JeuCompletT10$dif_AvApIm_abs*100
JeuCompletT10$dif_AvApIp_abs<-JeuCompletT10$dif_AvApIp_abs*100

ggplot(JeuCompletT10)+geom_point(aes(essai,dif_vac_unvacM_abs,color=vac,shape=vac))+geom_hline(yintercept=seq(2,8,2),color='grey',size = 0.2)+ geom_vline(xintercept=seq(0.45,1.45,0.10),color='grey',size = 0.2)+  theme_classic() +geom_errorbar(aes(essai,ymin = dif_vac_unvacIm_abs, ymax = dif_vac_unvacIp_abs, color=vac),width = 0.02)+scale_shape_manual(values=c(1,2,0,6),guide="none")+ scale_color_manual(name="Vaccine  \n coverage \n \n Vac-vs-Unvac",values=palette2,guide=guide_legend(order=1))+guides(color=guide_legend(override.aes=list(shape=c(1,2,0,6),color=palette2)))+ geom_hline(yintercept=c(0),color='darkgrey', linetype="dashed")+labs(title= "5 years after vaccine introduction",x = "Strength of interaction", y = "Absolute NV-prevalence difference (%)")+ scale_y_continuous(limits=c(-0.0, 9),breaks=seq(0,8,2))+ new_scale_color() +new_scale("shape")+ 
  geom_point(aes(essai,dif_AvApM_abs,color=vac,shape=vac,fill=vac))+ geom_errorbar(aes(essai,ymin = dif_AvApIm_abs, ymax = dif_AvApIp_abs, color=vac),width = 0.02)+scale_shape_manual(values=c(21,24,22,25),guide="none")+ scale_color_manual(name="Post-vs-Pre",values=palette1[4:1],guide=guide_legend(order=2))+scale_fill_manual(values=palette1[4:1],guide="none")+ guides(color=guide_legend(override.aes=list(shape=c(21,24,22,25),fill=palette1[4:1])))+ geom_vline(xintercept=seq(0.45,1.45,0.10),color='grey',size = 0.2)



ggplot(JeuCompletT10) + 
  geom_vline(xintercept=seq(0.45,1.45,0.10),color='grey',size = 0.2)+geom_hline(yintercept=c(10^3,10^5,10^7,10^9),color='grey',size = 0.2)+theme_classic() +labs(x = "Strength of interaction", y = "")+labs(title= "5 years after vaccine introduction") + scale_x_continuous(limits=c(0.45, 1.55),breaks=seq(0.5,1.5,0.1))+ geom_point( aes(essai,nsujetvacUnvacM,color=vac,shape=vac,fill=vac))+ geom_errorbar(aes(essai,ymin = nsujetvacUnvacIm, ymax =nsujetvacUnvacIp, color=vac),width = 0.02)+ scale_color_manual(name="Vaccine  \n coverage (%) \n \n Vac-vs-Unvac",values=palette2,guide=guide_legend(override.aes=list(shape=c(1,2,0,6)),order=1))+scale_shape_manual(values=c(1,2,0,6),guide="none")+scale_y_log10(limits=c(5*10^2,5*10^9),breaks=c(10^3,10^5,10^7,10^9))+
  guides(color=guide_legend(override.aes=list(shape=c(1,2,0,65),color=palette2)))+new_scale_color() +new_scale("shape")+ 
  geom_point( aes(essai,sujetAvAptotM,color=vac,shape=vac,fill=vac))+ geom_errorbar(aes(essai,ymin = nsujetAvAptotIm, ymax = nsujetAvAptotIp, color=vac),width = 0.02)+scale_color_manual(name="Post-vs-Pre",values=palette1[4:1])+scale_shape_manual(values=c(21,24,22,25),guide="none")+scale_fill_manual(values=palette1[4:1],guide="none")+ guides(color=guide_legend(override.aes=list(shape=c(21,24,22,25),fill=palette1[4:1])))

#Figure S4 -10 years after vaccine introduction
JeuCompletT10<-JeuCompletMIC[which(JeuCompletMIC$year==10 & JeuCompletMIC$partner=="all"),]
JeuCompletT10<-cbind(JeuCompletT10,essai=JeuCompletT10$Interaction+JeuCompletT10$vac*0.1-0.05)
#JeuCompletT10$Interaction<-as.factor(JeuCompletT10$Interaction)
JeuCompletT10$vac<-as.factor(JeuCompletT10$vac*100)

JeuCompletT10$dif_vac_unvacM_abs<-JeuCompletT10$dif_vac_unvacM_abs*100
JeuCompletT10$dif_vac_unvacIm_abs<-JeuCompletT10$dif_vac_unvacIm_abs*100
JeuCompletT10$dif_vac_unvacIp_abs<-JeuCompletT10$dif_vac_unvacIp_abs*100
JeuCompletT10$dif_AvApM_abs<-JeuCompletT10$dif_AvApM_abs*100
JeuCompletT10$dif_AvApIm_abs<-JeuCompletT10$dif_AvApIm_abs*100
JeuCompletT10$dif_AvApIp_abs<-JeuCompletT10$dif_AvApIp_abs*100

ggplot(JeuCompletT10)+geom_point(aes(essai,dif_vac_unvacM_abs,color=vac,shape=vac))+geom_hline(yintercept=seq(2,8,2),color='grey',size = 0.2)+ geom_vline(xintercept=seq(0.45,1.45,0.10),color='grey',size = 0.2)+  theme_classic() +geom_errorbar(aes(essai,ymin = dif_vac_unvacIm_abs, ymax = dif_vac_unvacIp_abs, color=vac),width = 0.02)+scale_shape_manual(values=c(1,2,0,6),guide="none")+ scale_color_manual(name="Vaccine  \n coverage \n \n Vac-vs-Unvac",values=palette2,guide=guide_legend(order=1))+guides(color=guide_legend(override.aes=list(shape=c(1,2,0,6),color=palette2)))+ geom_hline(yintercept=c(0),color='darkgrey', linetype="dashed")+labs(title= "10 years after vaccine introduction",x = "Strength of interaction", y = "Absolute NV-prevalence difference (%)")+ scale_y_continuous(limits=c(-0.0, 9),breaks=seq(0,8,2))+ new_scale_color() +new_scale("shape")+ 
  geom_point(aes(essai,dif_AvApM_abs,color=vac,shape=vac,fill=vac))+ geom_errorbar(aes(essai,ymin = dif_AvApIm_abs, ymax = dif_AvApIp_abs, color=vac),width = 0.02)+scale_shape_manual(values=c(21,24,22,25),guide="none")+ scale_color_manual(name="Post-vs-Pre",values=palette1[4:1],guide=guide_legend(order=2))+scale_fill_manual(values=palette1[4:1],guide="none")+ guides(color=guide_legend(override.aes=list(shape=c(21,24,22,25),fill=palette1[4:1])))+ geom_vline(xintercept=seq(0.45,1.45,0.10),color='grey',size = 0.2)



ggplot(JeuCompletT10) + 
  geom_vline(xintercept=seq(0.45,1.45,0.10),color='grey',size = 0.2)+geom_hline(yintercept=c(10^3,10^5,10^7,10^9),color='grey',size = 0.2)+theme_classic() +labs(x = "Strength of interaction", y = "")+labs(title= "10 years after vaccine introduction") + scale_x_continuous(limits=c(0.45, 1.55),breaks=seq(0.5,1.5,0.1))+ geom_point( aes(essai,nsujetvacUnvacM,color=vac,shape=vac,fill=vac))+ geom_errorbar(aes(essai,ymin = nsujetvacUnvacIm, ymax =nsujetvacUnvacIp, color=vac),width = 0.02)+ scale_color_manual(name="Vaccine  \n coverage (%) \n \n Vac-vs-Unvac",values=palette2,guide=guide_legend(override.aes=list(shape=c(1,2,0,6)),order=1))+scale_shape_manual(values=c(1,2,0,6),guide="none")+scale_y_log10(limits=c(5*10^2,5*10^9),breaks=c(10^3,10^5,10^7,10^9))+
  guides(color=guide_legend(override.aes=list(shape=c(1,2,0,65),color=palette2)))+new_scale_color() +new_scale("shape")+ 
  geom_point( aes(essai,sujetAvAptotM,color=vac,shape=vac,fill=vac))+ geom_errorbar(aes(essai,ymin = nsujetAvAptotIm, ymax = nsujetAvAptotIp, color=vac),width = 0.02)+scale_color_manual(name="Post-vs-Pre",values=palette1[4:1])+scale_shape_manual(values=c(21,24,22,25),guide="none")+scale_fill_manual(values=palette1[4:1],guide="none")+ guides(color=guide_legend(override.aes=list(shape=c(21,24,22,25),fill=palette1[4:1])))

#Figure S4 -20 years after vaccine introduction

JeuCompletT10<-JeuCompletMIC[which(JeuCompletMIC$year==20 & JeuCompletMIC$partner=="all"),]
JeuCompletT10<-cbind(JeuCompletT10,essai=JeuCompletT10$Interaction+JeuCompletT10$vac*0.1-0.05)
#JeuCompletT10$Interaction<-as.factor(JeuCompletT10$Interaction)
JeuCompletT10$vac<-as.factor(JeuCompletT10$vac*100)

JeuCompletT10$dif_vac_unvacM_abs<-JeuCompletT10$dif_vac_unvacM_abs*100
JeuCompletT10$dif_vac_unvacIm_abs<-JeuCompletT10$dif_vac_unvacIm_abs*100
JeuCompletT10$dif_vac_unvacIp_abs<-JeuCompletT10$dif_vac_unvacIp_abs*100
JeuCompletT10$dif_AvApM_abs<-JeuCompletT10$dif_AvApM_abs*100
JeuCompletT10$dif_AvApIm_abs<-JeuCompletT10$dif_AvApIm_abs*100
JeuCompletT10$dif_AvApIp_abs<-JeuCompletT10$dif_AvApIp_abs*100

ggplot(JeuCompletT10)+geom_point(aes(essai,dif_vac_unvacM_abs,color=vac,shape=vac))+geom_hline(yintercept=seq(2,8,2),color='grey',size = 0.2)+ geom_vline(xintercept=seq(0.45,1.45,0.10),color='grey',size = 0.2)+  theme_classic() +geom_errorbar(aes(essai,ymin = dif_vac_unvacIm_abs, ymax = dif_vac_unvacIp_abs, color=vac),width = 0.02)+scale_shape_manual(values=c(1,2,0,6),guide="none")+ scale_color_manual(name="Vaccine  \n coverage \n \n Vac-vs-Unvac",values=palette2,guide=guide_legend(order=1))+guides(color=guide_legend(override.aes=list(shape=c(1,2,0,6),color=palette2)))+ geom_hline(yintercept=c(0),color='darkgrey', linetype="dashed")+labs(title= "20 years after vaccine introduction",x = "Strength of interaction", y = "Absolute NV-prevalence difference (%)")+ scale_y_continuous(limits=c(-0.0, 9),breaks=seq(0,8,2))+ new_scale_color() +new_scale("shape")+ 
  geom_point(aes(essai,dif_AvApM_abs,color=vac,shape=vac,fill=vac))+ geom_errorbar(aes(essai,ymin = dif_AvApIm_abs, ymax = dif_AvApIp_abs, color=vac),width = 0.02)+scale_shape_manual(values=c(21,24,22,25),guide="none")+ scale_color_manual(name="Post-vs-Pre",values=palette1[4:1],guide=guide_legend(order=2))+scale_fill_manual(values=palette1[4:1],guide="none")+ guides(color=guide_legend(override.aes=list(shape=c(21,24,22,25),fill=palette1[4:1])))+ geom_vline(xintercept=seq(0.45,1.45,0.10),color='grey',size = 0.2)



ggplot(JeuCompletT10) + 
  geom_vline(xintercept=seq(0.45,1.45,0.10),color='grey',size = 0.2)+geom_hline(yintercept=c(10^3,10^5,10^7,10^9),color='grey',size = 0.2)+theme_classic() +labs(x = "Strength of interaction", y = "")+labs(title= "20 years after vaccine introduction") + scale_x_continuous(limits=c(0.45, 1.55),breaks=seq(0.5,1.5,0.1))+ geom_point( aes(essai,nsujetvacUnvacM,color=vac,shape=vac,fill=vac))+ geom_errorbar(aes(essai,ymin = nsujetvacUnvacIm, ymax =nsujetvacUnvacIp, color=vac),width = 0.02)+ scale_color_manual(name="Vaccine  \n coverage (%) \n \n Vac-vs-Unvac",values=palette2,guide=guide_legend(override.aes=list(shape=c(1,2,0,6)),order=1))+scale_shape_manual(values=c(1,2,0,6),guide="none")+scale_y_log10(limits=c(5*10^2,5*10^9),breaks=c(10^3,10^5,10^7,10^9))+
  guides(color=guide_legend(override.aes=list(shape=c(1,2,0,65),color=palette2)))+new_scale_color() +new_scale("shape")+ 
  geom_point( aes(essai,sujetAvAptotM,color=vac,shape=vac,fill=vac))+ geom_errorbar(aes(essai,ymin = nsujetAvAptotIm, ymax = nsujetAvAptotIp, color=vac),width = 0.02)+scale_color_manual(name="Post-vs-Pre",values=palette1[4:1])+scale_shape_manual(values=c(21,24,22,25),guide="none")+scale_fill_manual(values=palette1[4:1],guide="none")+ guides(color=guide_legend(override.aes=list(shape=c(21,24,22,25),fill=palette1[4:1])))

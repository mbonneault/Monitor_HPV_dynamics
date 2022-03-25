rm(list = ls())

library(ggplot2)
library(ggnewscale)


#### Data formatting ----

#setwd("Repositeries")
JeuComplet<-read.csv2("JeuComplet_tot_Prev_NV_V_part.csv")
JeuComplet<-cbind(JeuComplet,dif_vac_unvac=JeuComplet$val_prevV_Vac-JeuComplet$val_prevV_unVac)
JeuComplet<-cbind(JeuComplet,dif_AvAp=JeuComplet$val_prevVtot-JeuComplet$val_prevV_5)


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
  nsujetvac=c(nsujetvac,(kvac[pos]+1)/kvac[pos]*2.802^2/(2*(asin(sqrt(JeuComplet$val_prevV_Vac[int]))-asin(sqrt(JeuComplet$val_prevV_unVac[int]))))^2)
  nsujetunvac=c(nsujetunvac,kvac[pos]*nsujetvac[int])
  nsujetvactot=c(nsujetvactot,4*2.802^2/(2*(asin(sqrt(JeuComplet$val_prevV_Vac[int]))-asin(sqrt(JeuComplet$val_prevV_unVac[int]))))^2)
  nsujetAvAptot=c(nsujetAvAptot,4*2.802^2/(2*(asin(sqrt(JeuComplet$val_prevVtot[int]))-asin(sqrt(JeuComplet$val_prevV_5[int]))))^2)
  nsujetAvAp_min=c(nsujetAvAp_min,2.802^2/4/(asin(sqrt(JeuComplet$val_prevVtot[int]))-asin(sqrt(JeuComplet$val_prevV_5[int])))^2)
}

epi.ssxsectn(JeuComplet$val_prevNVtot[int], JeuComplet$val_prevNV_5[int], n=NA, power=0.8, r = 1, design = 0.5, sided.test = 2,conf.level = 0.95)

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
JeuCompletMIC<-JeuCompletMIC[which(JeuCompletMIC$Interaction>0.4 & (JeuCompletMIC$vac==0.2 |JeuCompletMIC$vac==0.4| JeuCompletMIC$vac==0.6|JeuCompletMIC$vac==0.8)),]



#selection of neutral interaction
JeuCompletT10<-JeuCompletMIC[which((JeuCompletMIC$Interaction==1) & (JeuCompletMIC$partner=="all")),]
JeuCompletT10$vac<-as.factor(JeuCompletT10$vac*100)
JeuCompletT10$dif_vac_unvacM<-JeuCompletT10$dif_vac_unvacM*100
JeuCompletT10$dif_vac_unvacIm<-JeuCompletT10$dif_vac_unvacIm*100
JeuCompletT10$dif_vac_unvacIp<-JeuCompletT10$dif_vac_unvacIp*100
JeuCompletT10$dif_AvApM<-JeuCompletT10$dif_AvApM*100
JeuCompletT10$dif_AvApIm<-JeuCompletT10$dif_AvApIm*100
JeuCompletT10$dif_AvApIp<-JeuCompletT10$dif_AvApIp*100

#Figure 4A
ggplot(JeuCompletT10)+geom_point(aes(year,dif_vac_unvacM,color=vac,shape=vac))+geom_line(aes(year,dif_vac_unvacM,color=vac))+geom_hline(yintercept=seq(-10,-2,2),color='grey',size = 0.2) +geom_vline(xintercept=15,color='grey',size = 0.2, linetype="dashed") +labs(title= "Vac-vs-Unvac comparison",x = "Years after vaccine introduction", y = "V-prevalence difference (%)") +geom_errorbar(aes(year,ymin = dif_vac_unvacIm, ymax = dif_vac_unvacIp, color=vac),width = 0.1)+ scale_color_manual(name="Vaccine  \n coverage (%) \n \n Vac-vs-Unvac",values=palette2,guide=guide_legend(override.aes=list(shape=c(1,2,0,6)),order=1))+ geom_hline(yintercept=c(0),color='red', linetype="dashed")+ scale_y_continuous(limits=c(-11.0,0.1),breaks=seq(-10,0,2))+scale_shape_manual(values=c(1,2,0,6),guide=F)+ theme_bw()

#Figure 4B
ggplot(JeuCompletT10)+ geom_hline(yintercept=c(0),color='red', linetype="dashed")+geom_hline(yintercept=seq(-10,-2,2),color='grey',size = 0.2) +labs(title= "Post-vs-Pre comparison",x = "Years after vaccine introduction", y = "V-prevalence difference (%)")+geom_vline(xintercept=15,color='grey',size = 0.2, linetype="dashed") + geom_point(aes(year,dif_AvApM,color=vac,shape=vac,fill=vac))+geom_line(aes(year,dif_AvApM,color=vac))+ geom_errorbar(aes(year,ymin = dif_AvApIm, ymax = dif_AvApIp, color=vac),width = 0.1)+scale_shape_manual(values=c(21,24,22,25),guide=FALSE)+ scale_color_manual(name="Vaccine  \n coverage (%) \n \n Post-vs-Pre",values=palette1[4:1])+scale_fill_manual(values=palette1[4:1],guide=FALSE)+guides(color=guide_legend(override.aes=list(shape=c(21,24,22,25),fill=palette1[4:1])))+ theme_bw()+ scale_y_continuous(limits=c(-11.0,0.1),breaks=seq(-10,0,2))


#Figure 4C
ggplot(JeuCompletT10)+labs(title= "Vac-vs-Unvac comparison",x = "Years after vaccine introduction", y = "Sample size required") +geom_vline(xintercept=15,color='grey',size = 0.2, linetype="dashed") + geom_hline(yintercept=c(10^2,5*10^2,10^3,5*10^3,10^4),color='grey',size = 0.2)+geom_point(aes(year,nsujetvacUnvacM,color=vac,shape=vac))+geom_line(aes(year,nsujetvacUnvacM,color=vac))+geom_errorbar(aes(year,ymin = nsujetvacUnvacIm, ymax = nsujetvacUnvacIp, color=vac),width = 0.02)+ scale_color_manual(name="Vaccine  \n coverage(%) \n \n Vac-vs-Unvac",values=palette2,guide=guide_legend(override.aes=list(shape=c(1,2,0,6)),order=1))+ scale_y_log10(limits=c(80,2*10^5),breaks=c(10^2,5*10^2,10^3,5*10^3,10^4,5*10^4,10^5))+scale_shape_manual(values=c(1,2,0,6),guide=F)+ theme_bw() 

#Figure 4D  
ggplot(JeuCompletT10)+labs(title= "Post-vs-Pre comparison",x = "Years after vaccine introduction", y = "Sample size required")+geom_vline(xintercept=15,color='grey',size = 0.2, linetype="dashed") + geom_hline(yintercept=c(10^2,5*10^2,10^3,5*10^3,10^4),color='grey',size = 0.2)+ geom_point( aes(year,sujetAvAptotM,color=vac,shape=vac,fill=vac))+geom_line( aes(year,sujetAvAptotM,color=vac))+ geom_errorbar(aes(year,ymin = nsujetAvAptotIm, ymax = nsujetAvAptotIp, color=vac),width = 0.02)+scale_color_manual(name="Post-vs-Pre",values=palette1[4:1])+scale_shape_manual(values=c(21,24,22,25),guide=FALSE)+scale_fill_manual(values=palette1[4:1],guide=FALSE)+ scale_y_log10(limits=c(80,2*10^5),breaks=c(10^2,5*10^2,10^3,5*10^3,10^4,5*10^4,10^5))+ theme_bw()

  
 

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

#figure 6 V genotypes
ggplot(JeuCompletT10)+geom_point(aes(year,dif_vac_unvacM,color=vac,shape=vac))+geom_line(aes(year,dif_vac_unvacM,color=vac))+geom_hline(yintercept=seq(-100,0,20),color='grey',size = 0.2) +geom_vline(xintercept=15,color='grey',size = 0.2, linetype="dashed")+    theme_bw() +labs(x = "Years after vaccine introduction", y = "V-prevalence differences (%)") +geom_errorbar(aes(year,ymin = dif_vac_unvacIm, ymax = dif_vac_unvacIp, color=vac),width = 0.1)+ scale_color_manual(name="Vaccine  \n coverage (%) \n \n Vac-vs-Unvac",values=palette2,guide=guide_legend(override.aes=list(shape=c(1,2,0,6)),order=1))+ geom_hline(yintercept=c(0),color='red', linetype="dashed")+labs(title= " V genotypes")+ scale_y_continuous(limits=c(-100.0, 05),breaks=seq(-100,0,20))+scale_shape_manual(values=c(1,2,0,6),guide=F)+ new_scale_color() +new_scale("shape")+
  geom_point(aes(year,dif_AvApM,color=vac,shape=vac,fill=vac))+geom_line(aes(year,dif_AvApM,color=vac))+ geom_errorbar(aes(year,ymin = dif_AvApIm, ymax = dif_AvApIp, color=vac),width = 0.1)+scale_shape_manual(values=c(21,24,22,25),guide=FALSE)+ scale_color_manual(name="Post-vs-Pre",values=palette1[4:1],guide=guide_legend(order=2))+scale_fill_manual(values=palette1[4:1],guide=FALSE)+
  guides(color=guide_legend(override.aes=list(shape=c(21,24,22,25),fill=palette1[4:1])))+ facet_grid(. ~ Part)


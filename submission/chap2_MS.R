library(readr)
library(dplyr)
library(metafor)  
library(ggplot2)
library(stringr)
library(tidyverse) 
library(ncmeta)
library(viridis)
library(ggthemes)
library(LSD)
library(yardstick)
library(ggplot2)
library(RColorBrewer)
library(dplyr)
library(gplots)
library(tidyselect)
library(extrafont)
library(cowplot)
library(ncdf4)
library(scales)
library(ggpubr)
library(MAd)
library(Hmisc)
library("missMDA")
library(FactoMineR)
library(plotrix)


#read csv
gcme <- read.csv("~/data/gcme/MS_data/plot_data.csv")

#output vcmax prediction under eCO2
gcme_co2 <- subset(gcme,condition=="co2")

#create a new column named Prediction to insert to boxplot
gcme_co2$prediction <- "Prediction"

vcmax_co2_fig <- gcme_co2 %>%
  ggplot( aes(x=ecosystem, y=vcmax)) +
  geom_boxplot(alpha = 0.6, outlier.shape = NA)+
  geom_boxplot(data=gcme_co2,aes(x=prediction,y=pred_vcmax), alpha = 0.6, width = 0.5,color="red",outlier.shape = NA) +
  geom_jitter(size=2) +
  geom_hline( linetype = 'dotted',yintercept=0.0, size=0.5)+ ylim(-1,1)+
  theme_classic()+coord_flip()+theme(axis.text=element_text(size=12))+
  labs(x="", y=~paste(V[cmax]," response to ",eCO[2])) 

jmax_co2_fig <- gcme_co2 %>%
  ggplot( aes(x=ecosystem, y=jmax)) +
  geom_boxplot(alpha = 0.6,outlier.shape = NA)+
  geom_boxplot(data=gcme_co2,aes(x=prediction,y=pred_jmax), alpha = 0.6, width = 0.5,color="red",outlier.shape = NA) +
  geom_jitter(size=2) +
  geom_hline( linetype = 'dotted',yintercept=0.0, size=0.5)+ ylim(-1,1)+
  theme_classic()+coord_flip()+theme(axis.text=element_text(size=12))+
  labs(x="", y=~paste(J[max]," response to ",eCO[2])) 

jv_co2_fig <- gcme_co2 %>%
  ggplot( aes(x=ecosystem, y=jmax_vcmax)) +
  geom_boxplot(alpha = 0.6,outlier.shape = NA)+
  geom_boxplot(data=gcme_co2,aes(x=prediction,y=pred_jmax_vcmax), alpha = 0.6, width = 0.5,color="red",outlier.shape = NA) +
  geom_jitter(size=2) +
  geom_hline( linetype = 'dotted',yintercept=0.0, size=0.5)+ ylim(-1,1)+
  theme_classic()+coord_flip()+theme(axis.text=element_text(size=12))+
  labs(x="", y=~paste(J[max],"/",V[cmax]," response to ",eCO[2])) 

plot_grid(vcmax_co2_fig,jmax_co2_fig,jv_co2_fig,nrow=1,label_size = 15)+theme(plot.background=element_rect(fill="white", color="white"))

ggsave(paste("./output/chap2_co2.jpg",sep=""),width = 15, height = 5)

#output vcmax prediction under light
gcme_light <- subset(gcme,condition=="light")
gcme_light$type_name[gcme_light$type_name=="shade_to_sun"] <- "Shade to sun"
gcme_light$type_name[gcme_light$type_name=="low_to_high_light"] <- "Low to high light"

vcmax_light_fig <- gcme_light %>% ggplot( aes(x=type_name, y=vcmax)) +
  geom_jitter(size=2) +geom_jitter(aes(x=type_name, y=pred_vcmax),color="red",size=2)+
  geom_hline( linetype = 'dotted',yintercept=0.0, size=0.5)+ ylim(-0.5,1)+
  geom_hline(yintercept=mean(gcme_light$pred_vcmax),color="red", size=0.5)+
  theme_classic()+coord_flip()+theme(axis.text=element_text(size=12))+
  labs(x="", y=~paste(V[cmax]," response to light")) 

jmax_light_fig <- gcme_light %>% ggplot( aes(x=type_name, y=jmax)) +
  geom_jitter(size=2) +geom_jitter(aes(x=type_name, y=pred_jmax),color="red",size=2)+
  geom_hline( linetype = 'dotted',yintercept=0.0, size=0.5)+ ylim(-0.5,1)+
  geom_hline(yintercept=mean(gcme_light$pred_jmax),color="red", size=0.5)+
  theme_classic()+coord_flip()+theme(axis.text=element_text(size=12))+
  labs(x="", y=~paste(J[max]," response to light")) 

jv_light_fig <- gcme_light %>% ggplot( aes(x=type_name, y=jmax_vcmax)) +
  geom_jitter(size=2) +geom_jitter(aes(x=type_name, y=pred_jmax_vcmax),color="red",size=2)+
  geom_hline( linetype = 'dotted',yintercept=0.0, size=0.5)+ ylim(-0.5,1)+
  geom_hline(yintercept=mean(gcme_light$pred_jmax_vcmax),color="red", size=0.5)+
  theme_classic()+coord_flip()+theme(axis.text=element_text(size=12))+
  labs(x="", y=~paste(J[max],"/",V[cmax]," response to light")) 

plot_grid(vcmax_light_fig,jmax_light_fig,jv_light_fig,nrow=1,label_size = 15)+theme(plot.background=element_rect(fill="white", color="white"))

ggsave(paste("./output/chap2_light.jpg",sep=""),width = 15, height = 5)

#output vcmax25 prediction under warming
gcme_warming <- subset(gcme,condition=="warming")
#create a new type name for below figure
gcme_warming$prediction <- "Prediction"
gcme_warming$observation <- "Observation"

w1 <- gcme_warming %>%
  ggplot( aes(x=observation, y=vcmax)) +
  geom_jitter(size=2,color="black") +
  geom_jitter(aes(x=prediction, y=pred_vcmax),size=2,color="red") +
  geom_hline(  linetype = 'dotted',yintercept=0.0, size=0.5)+ ylim(-0.15,0.05)+
  labs(x=" ", y=~paste(V[cmax25]," response to warming")) +
  theme_classic()+coord_flip()+theme(axis.text=element_text(size=12))

w2 <- gcme_warming %>%
  ggplot( aes(x=observation, y=jmax)) +
  geom_jitter(size=2,color="black") +
  geom_jitter(aes(x=prediction, y=pred_jmax),size=2,color="red") +
  geom_hline(  linetype = 'dotted',yintercept=0.0, size=0.5)+ ylim(-0.15,0.05)+
  labs(x=" ", y=~paste(J[max25]," response to warming")) +
  theme_classic()+coord_flip()+theme(axis.text=element_text(size=12))

w3 <- gcme_warming %>%
  ggplot( aes(x=observation, y=jmax_vcmax)) +
  geom_jitter(size=2,color="black") +
  geom_jitter(aes(x=prediction, y=pred_jmax_vcmax),size=2,color="red") +
  geom_hline( linetype = 'dotted', yintercept=0.0, size=0.5)+ ylim(-0.15,0.05)+
  labs(x=" ", y=~paste(J[max25],"/",V[cmax25],," response to warming")) +
  theme_classic()+coord_flip()+theme(axis.text=element_text(size=12))

plot_grid(w1,w2,w3,nrow=1,label_size = 15)+theme(plot.background=element_rect(fill="white", color="white"))

ggsave(paste("./output/chap2_warming.jpg",sep=""),width = 15, height = 5)

#N fertilization effect on eCO2-vcmax
gcme_nfer <-subset(gcme,type_name=="highN"|type_name=="lowN"|type_name=="No_fertilization"|type_name=="Fertilization")

#7 sites are low vs. high N; 7 sites are without vs. with N. Needs to divide them in colors
gcme_nfer %>% group_by(type_name)  %>% summarise(number = n())

#creat box name
gcme_nfer$type_box[gcme_nfer$type_name=="highN"|gcme_nfer$type_name=="Fertilization"] <- "N fertilization"
gcme_nfer$type_box[gcme_nfer$type_name=="lowN"|gcme_nfer$type_name=="No_fertilization"] <- "Non N fertilization"

f1 <- gcme_nfer %>%
  ggplot( aes(x=type_box, y=vcmax)) +
  geom_boxplot(alpha=0.5,width = 0.7, outlier.shape = NA)+
  geom_jitter(data=subset(gcme_nfer,type_name=="Fertilization"|type_name=="No_fertilization"),aes(x=type_box, y=vcmax),size=2,color="black")+
  geom_jitter(data=subset(gcme_nfer,type_name=="highN"|type_name=="lowN"),aes(x=type_box, y=vcmax),size=2,color="blue")+
  geom_hline(  linetype = 'dotted',yintercept=0.0, size=0.5)+ ylim(-1.5,0.5)+
  labs(x="", y=~paste(V[cmax]," response to eCO2")) +
  theme_classic()+coord_flip()+theme(axis.text=element_text(size=12))

f2 <- gcme_nfer %>%
  ggplot( aes(x=type_box, y=jmax)) +
  geom_boxplot(alpha=0.5,width = 0.7, outlier.shape = NA)+
  geom_jitter(data=subset(gcme_nfer,type_name=="Fertilization"|type_name=="No_fertilization"),aes(x=type_box, y=jmax),size=2,color="black")+
  geom_jitter(data=subset(gcme_nfer,type_name=="highN"|type_name=="lowN"),aes(x=type_box, y=jmax),size=2,color="blue")+
  geom_hline( linetype = 'dotted', yintercept=0.0, size=0.5)+ ylim(-1.5,1)+
  labs(x="", y=~paste(J[max]," response to eCO2")) +
  theme_classic()+coord_flip()+theme(axis.text=element_text(size=12))

f3 <- gcme_nfer %>%
  ggplot( aes(x=type_box, y=jmax_vcmax)) +
  geom_boxplot(alpha=0.5,width = 0.7, outlier.shape = NA)+
  geom_jitter(data=subset(gcme_nfer,type_name=="Fertilization"|type_name=="No_fertilization"),aes(x=type_box, y=jmax_vcmax),size=2,color="black")+
  geom_jitter(data=subset(gcme_nfer,type_name=="highN"|type_name=="lowN"),aes(x=type_box, y=jmax_vcmax),size=2,color="blue")+
  geom_hline( linetype = 'dotted', yintercept=0.0, size=0.5)+ ylim(-0.5,1.2)+
  labs(x="", y=~paste(J[max],"/",V[cmax]," response to eCO2")) +
  theme_classic()+coord_flip()+theme(axis.text=element_text(size=12))

plot_grid(f1,f2,f3,nrow=1,label_size = 15)+theme(plot.background=element_rect(fill="white", color="white"))

ggsave(paste("./output/chap2_fertilization.jpg",sep=""),width = 15, height = 5)

#anova test
summary(aov(vcmax ~ type_box, data = gcme_nfer))
summary(aov(jmax ~ type_box, data = gcme_nfer))
summary(aov(jmax_vcmax ~ type_box, data = gcme_nfer))

#now, meta-analysis of eCO2 responses 
#only keep co2 + fertilization (where only evaluate co2 effect) plots 
gcme_meta <- subset(gcme,condition=="co2")

#make points labelled well
gcme_meta$type_name[gcme_meta$condition=="co2"] <- "others"
gcme_meta$type_name[gcme_meta$ecm_type=="Nfix"] <- "N-fixing"

gcme_meta %>% group_by(type_name)  %>% summarise(number = n())

t1 <- ggplot(gcme_meta,aes_string(x="jmax", y="vcmax")) +geom_hline(yintercept=0)+geom_vline(xintercept=0)+geom_point(aes(color=type_name),size=3)+
  stat_cor(label.y.npc="top", label.x.npc = "left",aes(label = paste(..rr.label.., ..p.label.., sep = "~`,`~")),size=7)+
  geom_smooth(color="black",method="lm",se=T)+
  theme(axis.text=element_text(size=25),axis.title=element_text(size=25,face="bold"),legend.position="none")+
  labs(y=~paste(V[cmax]))+labs(x=~paste(J[max]))

t1a <- ggplot(gcme_meta,aes_string(x="vcmax", y="jmax")) +geom_hline(yintercept=0)+geom_vline(xintercept=0)+geom_point(aes(color=type_name),size=3)+stat_cor(label.y.npc="top", label.x.npc = "left",aes(label = paste(..rr.label.., ..p.label.., sep = "~`,`~")),size=7)+
  geom_smooth(color="black",method="lm",se=T)+
  theme(axis.text=element_text(size=25),axis.title=element_text(size=25,face="bold"),legend.position="none")+
  labs(y=~paste(J[max]))+labs(x=~paste(V[cmax]))

t2 <- ggplot(gcme_meta,aes_string(x="nmass", y="vcmax")) +geom_hline(yintercept=0)+geom_vline(xintercept=0)+geom_point(aes(color=type_name),size=3)+stat_cor(label.y.npc="top", label.x.npc = "left",aes(label = paste(..rr.label.., ..p.label.., sep = "~`,`~")),size=7)+
  geom_smooth(color="black",method="lm",se=T)+
  theme(axis.text=element_text(size=25),axis.title=element_text(size=25,face="bold"),legend.position="none")+
  labs(y=~paste(V[cmax]))+labs(x=~paste(N[mass]))

t3 <- ggplot(gcme_meta,aes_string(x="narea", y="vcmax")) +geom_hline(yintercept=0)+geom_vline(xintercept=0)+geom_point(aes(color=type_name),size=3)+stat_cor(label.y.npc="top", label.x.npc = "left",aes(label = paste(..rr.label.., ..p.label.., sep = "~`,`~")),size=7)+
  geom_smooth(color="black",method="lm",se=T)+
  theme(axis.text=element_text(size=25),axis.title=element_text(size=25,face="bold"),legend.position="none")+
  labs(y=~paste(V[cmax]))+labs(x=~paste(N[area]))

t4 <- ggplot(gcme_meta,aes_string(x="LMA", y="vcmax")) +geom_hline(yintercept=0)+geom_vline(xintercept=0)+geom_point(aes(color=type_name),size=3)+stat_cor(label.y.npc="top", label.x.npc = "left",aes(label = paste(..rr.label.., ..p.label.., sep = "~`,`~")),size=7)+
  theme(axis.text=element_text(size=25),axis.title=element_text(size=25,face="bold"),legend.position="none")+
  labs(y=~paste(V[cmax]))+labs(x=~paste(LMA))

p1 <- ggplot(gcme_meta,aes_string(x="nmass", y="jmax")) +geom_hline(yintercept=0)+geom_vline(xintercept=0)+geom_point(aes(color=type_name),size=3)+stat_cor(label.y.npc="top", label.x.npc = "left",aes(label = paste(..rr.label.., ..p.label.., sep = "~`,`~")),size=7)+
  geom_smooth(color="black",method="lm",se=T)+
  theme(axis.text=element_text(size=25),axis.title=element_text(size=25,face="bold"),legend.position="none")+
  labs(y=~paste(J[max]))+labs(x=~paste(N[mass]))

p2 <- ggplot(gcme_meta,aes_string(x="narea", y="jmax")) +geom_hline(yintercept=0)+geom_vline(xintercept=0)+geom_point(aes(color=type_name),size=3)+stat_cor(label.y.npc="top", label.x.npc = "left",aes(label = paste(..rr.label.., ..p.label.., sep = "~`,`~")),size=7)+
  geom_smooth(color="black",method="lm",se=T)+
  theme(axis.text=element_text(size=25),axis.title=element_text(size=25,face="bold"),legend.position="none")+
  labs(y=~paste(J[max]))+labs(x=~paste(N[area]))

p3 <- ggplot(gcme_meta,aes_string(x="LMA", y="jmax")) +geom_hline(yintercept=0)+geom_vline(xintercept=0)+geom_point(aes(color=type_name),size=3)+stat_cor(label.y.npc="top", label.x.npc = "left",aes(label = paste(..rr.label.., ..p.label.., sep = "~`,`~")),size=7)+
  theme(axis.text=element_text(size=25),axis.title=element_text(size=25,face="bold"),legend.position="none")+
  labs(y=~paste(J[max]))+labs(x=~paste(LMA))

#now - for bi-variate relationship
b1 <- ggplot(gcme_meta,aes_string(x="anpp", y="vcmax")) +geom_hline(yintercept=0)+geom_vline(xintercept=0)+geom_point(aes(color=type_name),size=3)+stat_cor(label.y.npc="top", label.x.npc = "left",aes(label = paste(..rr.label.., ..p.label.., sep = "~`,`~")),size=7)+
  geom_smooth(color="black",method="lm",se=F,linetype = "dashed")+
  theme(axis.text=element_text(size=25),axis.title=element_text(size=25,face="bold"),legend.position="none")+
  labs(y=~paste(V[cmax]))+labs(x=~paste(ANPP))

b2 <- ggplot(gcme_meta,aes_string(x="bnpp", y="vcmax")) +geom_hline(yintercept=0)+geom_vline(xintercept=0)+geom_point(aes(color=type_name),size=3)+stat_cor(label.y.npc="top", label.x.npc = "left",aes(label = paste(..rr.label.., ..p.label.., sep = "~`,`~")),size=7)+
  geom_smooth(color="black",method="lm",se=T)+
  theme(axis.text=element_text(size=25),axis.title=element_text(size=25,face="bold"),legend.position="none")+
  labs(y=~paste(V[cmax]))+labs(x=~paste(BNPP))

b3 <- ggplot(gcme_meta,aes_string(x="root_shoot_ratio", y="vcmax")) +geom_hline(yintercept=0)+geom_vline(xintercept=0)+geom_point(aes(color=type_name),size=3)+stat_cor(label.y.npc="top", label.x.npc = "left",aes(label = paste(..rr.label.., ..p.label.., sep = "~`,`~")),size=7)+
  geom_smooth(color="black",method="lm",se=F,linetype = "dashed")+
  theme(axis.text=element_text(size=25),axis.title=element_text(size=25,face="bold"),legend.position="none")+
  labs(y=~paste(V[cmax]))+labs(x=~paste(Root/Shoot))

b4 <- ggplot(gcme_meta,aes_string(x="lai", y="vcmax")) +geom_hline(yintercept=0)+geom_vline(xintercept=0)+geom_point(aes(color=type_name),size=3)+stat_cor(label.y.npc="top", label.x.npc = "left",aes(label = paste(..rr.label.., ..p.label.., sep = "~`,`~")),size=7)+
  geom_smooth(color="black",method="lm",se=T)+
  theme(axis.text=element_text(size=25),axis.title=element_text(size=25,face="bold"),legend.position="none")+
  labs(y=~paste(V[cmax]))+labs(x=~paste("LAI"))

b5 <- ggplot(gcme_meta,aes_string(x="soil_mineral_N", y="vcmax")) +geom_hline(yintercept=0)+geom_vline(xintercept=0)+geom_point(aes(color=type_name),size=3)+stat_cor(label.y.npc="top", label.x.npc = "left",aes(label = paste(..rr.label.., ..p.label.., sep = "~`,`~")),size=7)+
  theme(axis.text=element_text(size=25),axis.title=element_text(size=25,face="bold"),legend.position="none")+
  labs(y=~paste(V[cmax]))+labs(x=~paste("Soil inorganic N"))

c1 <- ggplot(gcme_meta,aes_string(x="anpp", y="jmax")) +geom_hline(yintercept=0)+geom_vline(xintercept=0)+geom_point(aes(color=type_name),size=3)+stat_cor(label.y.npc="top", label.x.npc = "left",aes(label = paste(..rr.label.., ..p.label.., sep = "~`,`~")),size=7)+
  geom_smooth(color="black",method="lm",se=T)+
  theme(axis.text=element_text(size=25),axis.title=element_text(size=25,face="bold"),legend.position="none")+
  labs(y=~paste(J[max]))+labs(x=~paste(ANPP))

c2 <- ggplot(gcme_meta,aes_string(x="bnpp", y="jmax")) +geom_hline(yintercept=0)+geom_vline(xintercept=0)+geom_point(aes(color=type_name),size=3)+stat_cor(label.y.npc="top", label.x.npc = "left",aes(label = paste(..rr.label.., ..p.label.., sep = "~`,`~")),size=7)+
  geom_smooth(color="black",method="lm",se=T)+
  theme(axis.text=element_text(size=25),axis.title=element_text(size=25,face="bold"),legend.position="none")+
  labs(y=~paste(J[max]))+labs(x=~paste(BNPP))

c3 <- ggplot(gcme_meta,aes_string(x="root_shoot_ratio", y="jmax")) +geom_hline(yintercept=0)+geom_vline(xintercept=0)+geom_point(aes(color=type_name),size=3)+stat_cor(label.y.npc="top", label.x.npc = "left",aes(label = paste(..rr.label.., ..p.label.., sep = "~`,`~")),size=7)+
  theme(axis.text=element_text(size=25),axis.title=element_text(size=25,face="bold"),legend.position="none")+
  labs(y=~paste(J[max]))+labs(x=~paste(Root/Shoot))

c4 <- ggplot(gcme_meta,aes_string(x="lai", y="jmax")) +geom_hline(yintercept=0)+geom_vline(xintercept=0)+geom_point(aes(color=type_name),size=3)+stat_cor(label.y.npc="top", label.x.npc = "left",aes(label = paste(..rr.label.., ..p.label.., sep = "~`,`~")),size=7)+
  theme(axis.text=element_text(size=25),axis.title=element_text(size=25,face="bold"),legend.position="none")+
  labs(y=~paste(J[max]))+labs(x=~paste("LAI"))

c5 <- ggplot(gcme_meta,aes_string(x="soil_mineral_N", y="jmax")) +geom_hline(yintercept=0)+geom_vline(xintercept=0)+geom_point(aes(color=type_name),size=3)+stat_cor(label.y.npc="top", label.x.npc = "left",aes(label = paste(..rr.label.., ..p.label.., sep = "~`,`~")),size=7)+
  theme(axis.text=element_text(size=25),axis.title=element_text(size=25,face="bold"),legend.position="none")+
  labs(y=~paste(J[max]))+labs(x=~paste("Soil inorganic N"))

#check p-value of leaf traits and output those figures with p < 0.05
gcme_meta_metrics <- gcme_meta[,c("vcmax","jmax","LMA","narea","nmass",
                                  "bnpp","root_shoot_ratio","soil_mineral_N",
                                  "anpp","lai")]
res2 <- rcorr(as.matrix(gcme_meta_metrics))

leaf_pvalue <- as.data.frame(res2$P[,c("narea","LMA","nmass")])
leaf_pvalue[leaf_pvalue<0.05] <- "yes"
leaf_pvalue

#additionally plot Narea vs. Nmass,  narea vs. anpp, LMA vs. root/shoot, LMA vs LAI, Nmass vs. bnpp

d1 <- ggplot(gcme_meta,aes_string(x="nmass", y="narea")) +geom_hline(yintercept=0)+geom_vline(xintercept=0)+
  geom_point(aes(color=type_name),size=3)+
  stat_cor(label.y.npc="top", label.x.npc = "left",aes(label = paste(..rr.label.., ..p.label.., sep = "~`,`~")),size=7)+
  geom_smooth(color="black",method="lm",se=T)+
  theme(axis.text=element_text(size=25),axis.title=element_text(size=25,face="bold"),legend.position="none")+
  labs(x=~paste(N[mass]))+labs(y=~paste(N[area]))

d2 <- ggplot(gcme_meta,aes_string(x="anpp", y="narea")) +geom_hline(yintercept=0)+
  geom_vline(xintercept=0)+geom_point(aes(color=type_name),size=3)+
  stat_cor(label.y.npc="top", label.x.npc = "left",aes(label = paste(..rr.label.., ..p.label.., sep = "~`,`~")),size=7)+
  geom_smooth(color="black",method="lm",se=T)+xlim(-0.1,0.7)+
  theme(axis.text=element_text(size=25),axis.title=element_text(size=25,face="bold"),legend.position="none")+
  labs(x=~paste("ANPP"))+labs(y=~paste(N[area]))

d3 <- ggplot(gcme_meta,aes_string(x="root_shoot_ratio", y="LMA")) +geom_hline(yintercept=0)+
  geom_vline(xintercept=0)+geom_point(aes(color=type_name),size=3)+
  stat_cor(label.y.npc="top", label.x.npc = "left",aes(label = paste(..rr.label.., ..p.label.., sep = "~`,`~")),size=7)+
  geom_smooth(color="black",method="lm",se=T)+
  theme(axis.text=element_text(size=25),axis.title=element_text(size=25,face="bold"),legend.position="none")+
  labs(x=~paste("Root / shoot"))+labs(y=~paste(LMA))

d4 <- ggplot(gcme_meta,aes_string(x="lai", y="LMA")) +geom_hline(yintercept=0)+geom_vline(xintercept=0)+geom_point(aes(color=type_name),size=3)+stat_cor(label.y.npc="top", label.x.npc = "left",aes(label = paste(..rr.label.., ..p.label.., sep = "~`,`~")),size=7)+
  geom_smooth(color="black",method="lm",se=T)+
  theme(axis.text=element_text(size=25),axis.title=element_text(size=25,face="bold"),legend.position="none")+
  labs(x=~paste(LAI))+labs(y=~paste(LMA))

d5 <- ggplot(gcme_meta,aes_string(x="bnpp", y="nmass")) +geom_hline(yintercept=0)+geom_vline(xintercept=0)+geom_point(aes(color=type_name),size=3)+stat_cor(label.y.npc="top", label.x.npc = "left",aes(label = paste(..rr.label.., ..p.label.., sep = "~`,`~")),size=7)+
  geom_smooth(color="black",method="lm",se=T)+
  theme(axis.text=element_text(size=25),axis.title=element_text(size=25,face="bold"),legend.position="none")+
  labs(x=~paste("BNPP"))+labs(y=~paste(N[mass]))


final1_legend <- ggplot(gcme_meta,aes_string(x="soil_mineral_N", y="vcmax")) +geom_hline(yintercept=0)+geom_vline(xintercept=0)+geom_point(aes(color=type_name),size=3)+stat_cor(label.y.npc="top", label.x.npc = "left",aes(label = paste(..rr.label.., ..p.label.., sep = "~`,`~")),size=7)+geom_smooth(color="black",method="lm",se=F)+theme(text = element_text(size=30),axis.text=element_text(size=25),axis.title=element_text(size=25,face="bold"))+
  labs(y=~paste(V[cmax]))+labs(x=~paste("Soil inorganic N"))+scale_colour_discrete(" ")

legend_info <- as_ggplot(get_legend(final1_legend))

white <- theme(plot.background=element_rect(fill="white", color="white"))

plot_grid(t1,t2,t3,t4,legend_info,
          t1a,p1,p2,p3,white,
          b1,b2,b3,b4,b5,
          c1,c2,c3,c4,c5,
          d1,d2,d3,d4,d5,
          nrow=5,
          labels = c('(a)','(b)','(c)','(d)',' ',
                     '(e)','(f)','(g)','(h)',' ',
                     '(i)','(j)','(k)','(l)','(m)',
                     '(n)','(o)','(p)','(q)','(r)',
                     '(s)','(t)','(u)','(v)','(w)'), label_size = 23)+theme(plot.background=element_rect(fill="white", color="white"))
ggsave(paste("./output/chap2_meta.jpg",sep=""),width = 25, height = 25)

#PCA figure

#
test <- gcme_meta[,c("lai","vcmax","jmax","narea","LMA","nmass","bnpp","anpp","root_shoot_ratio")]
ncomp <- estim_ncpPCA(test)
res.imp <- imputePCA(test, ncp = ncomp$ncp) #filled by pca function

res.pca <- PCA(res.imp$completeObs,axes=c(1,2))
ggsave(paste("./output/pca1.jpg",sep=""),width = 5, height = 5)
res.pca2 <- PCA(res.imp$completeObs,axes=c(3,4))
ggsave(paste("./output/pca2.jpg",sep=""),width = 5, height = 5)

#create a table for values and SE
mean(gcme_meta$vcmax,na.rm=TRUE);std.error(gcme_meta$vcmax)
mean(gcme_meta$jmax,na.rm=TRUE);std.error(gcme_meta$jmax)
mean(gcme_meta$narea,na.rm=TRUE);std.error(gcme_meta$narea)
mean(gcme_meta$LMA,na.rm=TRUE);std.error(gcme_meta$LMA)
mean(gcme_meta$nmass,na.rm=TRUE);std.error(gcme_meta$nmass)
mean(gcme_meta$bnpp,na.rm=TRUE);std.error(gcme_meta$bnpp)
mean(gcme_meta$anpp,na.rm=TRUE);std.error(gcme_meta$anpp)
mean(gcme_meta$root_shoot_ratio,na.rm=TRUE);std.error(gcme_meta$root_shoot_ratio)
mean(gcme_meta$soil_mineral_N,na.rm=TRUE);std.error(gcme_meta$soil_mineral_N)
mean(gcme_meta$lai,na.rm=TRUE);std.error(gcme_meta$lai)

t.test(gcme_meta$vcmax, mu = 0)
t.test(gcme_meta$jmax, mu = 0)
t.test(gcme_meta$narea, mu = 0)
t.test(gcme_meta$LMA, mu = 0)
t.test(gcme_meta$nmass, mu = 0)
t.test(gcme_meta$bnpp, mu = 0)
t.test(gcme_meta$anpp, mu = 0)
t.test(gcme_meta$root_shoot_ratio, mu = 0)
t.test(gcme_meta$soil_mineral_N, mu = 0)
t.test(gcme_meta$lai, mu = 0)

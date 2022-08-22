light <- read.csv("/Users/yunpeng/data/LxN_Greenhouse/data_sheets/LxN_data_withPhys.csv")
white <- theme(plot.background=element_rect(fill="white", color="white"))
#read ppfd
ppfd_smith <- mean(read.csv("/Users/yunpeng/data/gcme/kevin/forcing/climate_large/smith_light1.csv")$ppfd,na.rm=TRUE)*1000000 # in umol/m2/s
light$ppfd <- ppfd_smith

#convert
light$light <- light$ppfd*(1-light$shade.cover/100)
light$ln_vcmax <- log(light$Vcmax)
light$ln_light <- log(light$light)

#interactive regression
summary(lm(ln_vcmax~light*n.ppm,data=light))

light$n.ppm2 <- NA
light$n.ppm2[light$n.ppm==0] <- "Nfer=0"
light$n.ppm2[light$n.ppm==70] <- "Nfer=1 (70ppm)"
light$n.ppm2[light$n.ppm==210] <- "Nfer=2 (210ppm)"
light$n.ppm2[light$n.ppm==630] <- "Nfer=3 (630ppm)"

light$shade.cover2 <- NA
light$shade.cover2[light$shade.cover==0] <- "Shading=0%"
light$shade.cover2[light$shade.cover==30] <- "Shading=30%"
light$shade.cover2[light$shade.cover==50] <- "Shading=50%"
light$shade.cover2[light$shade.cover==80] <- "Shading=80%"

l0_n0 <-subset(light,n.ppm==0 &shade.cover==0)
l30_n0 <-subset(light,n.ppm==0 &shade.cover==30)
l50_n0 <-subset(light,n.ppm==0 &shade.cover==50)
l80_n0 <-subset(light,n.ppm==0 &shade.cover==80)

l0_n70 <-subset(light,n.ppm==70 &shade.cover==0)
l30_n70 <-subset(light,n.ppm==70 &shade.cover==30)
l50_n70 <-subset(light,n.ppm==70 &shade.cover==50)
l80_n70 <-subset(light,n.ppm==70 &shade.cover==80)

l0_n210 <-subset(light,n.ppm==210 &shade.cover==0)
l30_n210 <-subset(light,n.ppm==210 &shade.cover==30)
l50_n210 <-subset(light,n.ppm==210 &shade.cover==50)
l80_n210 <-subset(light,n.ppm==210 &shade.cover==80)

l0_n630 <-subset(light,n.ppm==630 &shade.cover==0)
l30_n630 <-subset(light,n.ppm==630 &shade.cover==30)
l50_n630 <-subset(light,n.ppm==630 &shade.cover==50)
l80_n630 <-subset(light,n.ppm==630 &shade.cover==80)


#have a look at light - vcmax, jmax reasonable

lall_n0 <-subset(light,n.ppm==0)
lall_n70 <-subset(light,n.ppm==70)
lall_n210 <-subset(light,n.ppm==210)
lall_n630 <-subset(light,n.ppm==630)

l1 <- lall_n0 %>% ggplot( aes(x=as.character(shade.cover2), y=Vcmax)) +geom_boxplot()+
  geom_point(alpha = 0.6, width = 0.5)+xlab("Nfer =0")

l2 <- lall_n70 %>% ggplot( aes(x=as.character(shade.cover2), y=Vcmax)) +geom_boxplot()+
  geom_point(alpha = 0.6, width = 0.5)+xlab("Nfer =70ppm")

l3 <- lall_n210 %>% ggplot( aes(x=as.character(shade.cover2), y=Vcmax)) +geom_boxplot()+
  geom_point(alpha = 0.6, width = 0.5)+xlab("Nfer =210ppm")

l4 <- lall_n630 %>% ggplot( aes(x=as.character(shade.cover2), y=Vcmax)) +geom_boxplot()+
  geom_point(alpha = 0.6, width = 0.5)+xlab("Nfer =530ppm")

plot_grid(l1,l2,l3,l4,nrow=1,label_size = 15)+
  theme(plot.background=element_rect(fill="white", color="white"))

ggsave(paste("~/data/output_gcme/colin/smith1.jpg",sep=""),width = 18, height = 5)

#calculate sensitivity coefficient
summary(lm(ln_vcmax~ln_light,data=lall_n0))$coef[2,1];summary(lm(ln_vcmax~ln_light,data=lall_n0))$coef[2,2]
summary(lm(ln_vcmax~ln_light,data=lall_n70))$coef[2,1];summary(lm(ln_vcmax~ln_light,data=lall_n70))$coef[2,2]
summary(lm(ln_vcmax~ln_light,data=lall_n210))$coef[2,1];summary(lm(ln_vcmax~ln_light,data=lall_n210))$coef[2,2]
summary(lm(ln_vcmax~ln_light,data=lall_n630))$coef[2,1];summary(lm(ln_vcmax~ln_light,data=lall_n630))$coef[2,2]

light %>% 
  ggplot(aes(x = ln_light, y = ln_vcmax)) +
  geom_smooth(aes(color=n.ppm2),alpha=0.1,method=lm, formula = y ~ x)+theme_classic()+theme(text = element_text(size=20))


#have a look at Nertilization
#have a look at light - vcmax, jmax reasonable

l0_nall <-subset(light,shade.cover==0)
l30_nall <-subset(light,shade.cover==30)
l50_nall <-subset(light,shade.cover==50)
l80_nall <-subset(light,shade.cover==80)

summary(lm(ln_vcmax~n.ppm,data=l0_nall))$coef[2,1];summary(lm(ln_vcmax~n.ppm,data=l0_nall))$coef[2,2]
summary(lm(ln_vcmax~n.ppm,data=l30_nall))$coef[2,1];summary(lm(ln_vcmax~n.ppm,data=l30_nall))$coef[2,2]
summary(lm(ln_vcmax~n.ppm,data=l50_nall))$coef[2,1];summary(lm(ln_vcmax~n.ppm,data=l50_nall))$coef[2,2]
summary(lm(ln_vcmax~n.ppm,data=l80_nall))$coef[2,1];summary(lm(ln_vcmax~n.ppm,data=l80_nall))$coef[2,2]

light %>% 
  ggplot(aes(x = n.ppm, y = ln_vcmax)) +
  geom_smooth(aes(color=shade.cover2),alpha=0.1,method=lm, formula = y ~ x)+theme_classic()+theme(text = element_text(size=20))

n1 <- l0_nall %>% ggplot( aes(x=as.character(n.ppm2), y=Vcmax)) +geom_boxplot()+
  geom_point(alpha = 0.6, width = 0.5)+xlab("0% shading")

n2 <- l30_nall %>% ggplot( aes(x=as.character(n.ppm2), y=Vcmax)) +geom_boxplot()+
  geom_point(alpha = 0.6, width = 0.5)+xlab("30% shading")

n3 <- l50_nall %>% ggplot( aes(x=as.character(n.ppm2), y=Vcmax)) +geom_boxplot()+
  geom_point(alpha = 0.6, width = 0.5)+xlab("50% shading")

n4 <- l80_nall %>% ggplot( aes(x=as.character(n.ppm2), y=Vcmax)) +geom_boxplot()+
  geom_point(alpha = 0.6, width = 0.5)+xlab("80% shading")

plot_grid(n1,n2,n3,n4,nrow=1,label_size = 15)+
  theme(plot.background=element_rect(fill="white", color="white"))

ggsave(paste("~/data/output_gcme/colin/smith2.jpg",sep=""),width = 20, height = 5)

#read csv for N fer
Nfer <- read.csv("/Users/yunpeng/data/NxI_soybean_phys/data/2021NxI_trait_data.csv")
#yi: inoculated
Nfer_yi <- subset(Nfer,inoc=="yi")
mylist <- list()

#high N / low N
for (i in 1:((nrow(Nfer_yi)/2))){
  mylist[[i]] <- log(Nfer_yi[2*i-1,7:56]/Nfer_yi[2*i,7:56])
}

Nfer_yi_logr <- do.call("rbind",mylist)
Nfer_yi_logr[sapply(Nfer_yi_logr, is.infinite)] <- NA

res2 <- rcorr(as.matrix(Nfer_yi_logr))
p_value <- as.data.frame(res2$P[,c("vcmax","jmax")])

vcmax_sig <- subset(p_value[order(p_value[,c("vcmax")] ),],vcmax<= 0.3)
vcmax_select <- rownames(vcmax_sig)
p <- list()
s <- list()


for(i in c(1:length(vcmax_select))){
  p[[i]] <- ggplot(Nfer_yi_logr,aes_string(x=vcmax_select[i], y="vcmax")) +
    geom_hline(yintercept=0)+geom_vline(xintercept=0)+
    geom_point(size=3)+
    stat_cor(aes(label = paste(..rr.label.., ..p.label.., sep = "~`,`~")))+
    geom_smooth(color="black",method="lm",se=F)+labs(y=~paste(V[cmax]))+
    theme(axis.text=element_text(size=20),axis.title=element_text(size=20,face="bold"),legend.position="none")}


plot_grid(p[[2]],p[[7]],p[[8]],p[[11]],p[[13]],p[[16]],p[[20]],p[[21]],p[[22]],p[[25]],p[[27]],p[[29]],nrow=3)+white
ggsave(paste("~/data/output_gcme/colin/smith3.jpg",sep=""),width = 20, height = 10)

#check weather Ci/Ca decreses with N fertilizations
boxplot(Nfer_yi_logr$ci.ca)
mean(Nfer_yi_logr$ci.ca,na.rm=TRUE);std.error(Nfer_yi_logr$ci.ca)
t.test(Nfer_yi_logr$ci.ca, mu = 0)


#without yi

Nfer <- read.csv("/Users/yunpeng/data/NxI_soybean_phys/data/2021NxI_trait_data.csv")
#ni: non-inoculated
Nfer_yi <- subset(Nfer,inoc=="ni")
mylist <- list()

#high N / low N
for (i in 1:((nrow(Nfer_yi)/2))){
  mylist[[i]] <- log(Nfer_yi[2*i-1,7:56]/Nfer_yi[2*i,7:56])
}

Nfer_yi_logr <- do.call("rbind",mylist)
Nfer_yi_logr[sapply(Nfer_yi_logr, is.infinite)] <- NA

res2 <- rcorr(as.matrix(Nfer_yi_logr))
p_value <- as.data.frame(res2$P[,c("vcmax","jmax")])

vcmax_sig <- subset(p_value[order(p_value[,c("vcmax")] ),],vcmax<= 0.3)
vcmax_select <- rownames(vcmax_sig)

p <- list()
s <- list()


for(i in c(1:length(vcmax_select))){
  p[[i]] <- ggplot(Nfer_yi_logr,aes_string(x=vcmax_select[i], y="vcmax")) +
    geom_hline(yintercept=0)+geom_vline(xintercept=0)+
    geom_point(size=3)+
    stat_cor(aes(label = paste(..rr.label.., ..p.label.., sep = "~`,`~")))+
    geom_smooth(color="black",method="lm",se=F)+labs(y=~paste(V[cmax]))+
    theme(axis.text=element_text(size=20),axis.title=element_text(size=20,face="bold"),legend.position="none")}

#check weather Ci/Ca decreses with N fertilizations
boxplot(Nfer_yi_logr$ci.ca)
mean(Nfer_yi_logr$ci.ca,na.rm=TRUE);std.error(Nfer_yi_logr$ci.ca)
t.test(Nfer_yi_logr$ci.ca, mu = 0)

plot_grid(p[[2]],p[[5]],p[[10]],p[[12]],p[[15]])+white
ggsave(paste("~/data/output_gcme/colin/smith4.jpg",sep=""),width = 15, height = 10)

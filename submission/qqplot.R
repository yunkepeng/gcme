library(readr) 
library(dplyr)
library(ggplot2)
library(cowplot)
library(ggpubr)
library(Hmisc)
library(missMDA)
library(FactoMineR)
library(plotrix)

# reset validation metrics info
white <- theme(plot.background=element_rect(fill="white", color="white"))

#read csv
gcme <- read.csv("~/data/gcme/MS_data/plot_data.csv")

#now, meta-analysis of eCO2 responses 
#only keep co2 + fertilization (where only evaluate co2 effect) plots 
gcme_meta <- subset(gcme,condition=="co2"|condition=="Fertilization"|
                      condition=="highN"|condition=="lowN")

#make points labelled well
gcme_meta$type_name[gcme_meta$condition=="highN"] <- "N-fertilization"
gcme_meta$type_name[gcme_meta$condition=="lowN"] <- "N-fertilization"
gcme_meta$type_name[gcme_meta$condition=="Fertilization"] <- "N-fertilization"
gcme_meta$type_name[gcme_meta$condition=="co2"] <- "others"
gcme_meta$type_name[gcme_meta$ecm_type=="Nfix"] <- "N-fixing"

#qq plot

# a for loop
qq.out <- qqplot(x=data1[,c("jmax")], y=data1[,c("vcmax")], plot.it=FALSE)

#create a function
df <- gcme_meta;y_name<-"vcmax";x_name<-"jmax";type_name<-"type_name"

data1 <- na.omit(df[,c(x_name,y_name,type_name)])
lm_output <- (lm(data1[,c(y_name)]~data1[,c(x_name)]))
data1$residual <- as.numeric(resid(lm_output))

ggplot(data1)+
  stat_qq(aes(sample = residual,color=type_name))+
  stat_qq_line(aes(sample = residual))+
  theme(axis.text=element_text(size=25),title=element_text(size=25),
        axis.title=element_text(size=25,face="bold"),legend.position="none")+
  xlab("Theoritical")+ylab("Residuals")+ ggtitle("Jmax vs. Vcmax")


qqplot_output <- function(df,y_name,x_name,type_name,lab_title){
  data1 <- na.omit(df[,c(x_name,y_name,type_name)])
  lm_output <- (lm(data1[,c(y_name)]~data1[,c(x_name)]))
  data1$residual <- as.numeric(resid(lm_output))
  
  output_fig <- ggplot(data1)+
    stat_qq(aes(sample = residual,color=type_name),size=3)+
    stat_qq_line(aes(sample = residual))+
    theme(axis.text=element_text(size=20),title=element_text(size=20),
          axis.title=element_text(size=20),legend.position="none")+
    xlab("Theoritical")+ylab("Residuals")+lab_title

  return(output_fig)
}

q1<- qqplot_output(gcme_meta,"vcmax","jmax","type_name",ggtitle(~paste(V[cmax]," vs. ",J[max],sep=" ")))
q2<- qqplot_output(gcme_meta,"vcmax","nmass","type_name",ggtitle(~paste(V[cmax]," vs. ",N[mass],sep=" ")))
q3<- qqplot_output(gcme_meta,"vcmax","narea","type_name",ggtitle(~paste(V[cmax]," vs. ",N[area],sep=" ")))
q4<- qqplot_output(gcme_meta,"vcmax","LMA","type_name",ggtitle(~paste(V[cmax]," vs. ",LMA,sep=" ")))
q5<- qqplot_output(gcme_meta,"jmax","vcmax","type_name",ggtitle(~paste(J[max]," vs. ",V[cmax],sep=" ")))
q6<- qqplot_output(gcme_meta,"jmax","nmass","type_name",ggtitle(~paste(J[max]," vs. ",N[mass],sep=" ")))
q7<- qqplot_output(gcme_meta,"jmax","narea","type_name",ggtitle(~paste(J[max]," vs. ",N[area],sep=" ")))
q8<- qqplot_output(gcme_meta,"jmax","LMA","type_name",ggtitle(~paste(J[max]," vs. ",LMA,sep=" ")))
q9<- qqplot_output(gcme_meta,"vcmax","anpp","type_name",ggtitle(~paste(V[cmax]," vs. ",ANPP,sep=" ")))
q10<- qqplot_output(gcme_meta,"vcmax","bnpp","type_name",ggtitle(~paste(V[cmax]," vs. ",BNPP,sep=" ")))
q11<- qqplot_output(gcme_meta,"vcmax","root_shoot_ratio","type_name",ggtitle(~paste(V[cmax]," vs. ",Root/Shoot,sep=" ")))
q12<- qqplot_output(gcme_meta,"vcmax","lai","type_name",ggtitle(~paste(V[cmax]," vs. ",LAI,sep=" ")))
q13<- qqplot_output(gcme_meta,"vcmax","soil_mineral_N","type_name",ggtitle(~paste(V[cmax]," vs. ","Soil N",sep=" ")))
q14<- qqplot_output(gcme_meta,"jmax","anpp","type_name",ggtitle(~paste(J[max]," vs. ",ANPP,sep=" ")))
q15<- qqplot_output(gcme_meta,"jmax","bnpp","type_name",ggtitle(~paste(J[max]," vs. ",BNPP,sep=" ")))
q16<- qqplot_output(gcme_meta,"jmax","root_shoot_ratio","type_name",ggtitle(~paste(J[max]," vs. ",Root/Shoot,sep=" ")))
q17<- qqplot_output(gcme_meta,"jmax","lai","type_name",ggtitle(~paste(J[max]," vs. ",LAI,sep=" ")))
q18<- qqplot_output(gcme_meta,"jmax","soil_mineral_N","type_name",ggtitle(~paste(J[max]," vs. ","Soil N",sep=" ")))
q19<- qqplot_output(gcme_meta,"LMA","root_shoot_ratio","type_name",ggtitle(~paste(LMA," vs. ",Root/Shoot,sep=" ")))
q20<- qqplot_output(gcme_meta,"LMA","lai","type_name",ggtitle(~paste(LMA," vs. ",LAI,sep=" ")))
q21<- qqplot_output(gcme_meta,"nmass","bnpp","type_name",ggtitle(~paste(N[mass]," vs. ",BNPP,sep=" ")))


#get legend
legend_fig<- ggplot(gcme_meta, aes( x= "jmax", y = "vcmax")) + 
  geom_point(aes(color=type_name),size=2)+
  theme(text = element_text(size=30),axis.text=element_text(size=25),axis.title=element_text(size=25,face="bold"))+scale_colour_discrete(" ")

legend_info <- as_ggplot(get_legend(legend_fig))

white <- theme(plot.background=element_rect(fill="white", color="white"))

plot_grid(q1,q2,q3,q4,legend_info,
          q5,q6,q7,q8,white,
          q9,q10,q11,q12,q13,
          q14,q15,q16,q17,q18,
          q19,q20,q21,white,white,
          nrow=5,
          labels = c('(a)','(b)','(c)','(d)',' ',
                     '(e)','(f)','(g)','(h)',' ',
                     '(i)','(j)','(k)','(l)','(m)',
                     '(n)','(o)','(p)','(q)','(r)',
                     '(s)','(t)','(u)', ' ' , ' '), label_size = 23)+
  theme(plot.background=element_rect(fill="white", color="white"))
ggsave(paste("/Users/yunpeng/yunkepeng/gcme_MS/output/chap2_meta_qqplot.jpg",sep=""),width = 25, height = 25)

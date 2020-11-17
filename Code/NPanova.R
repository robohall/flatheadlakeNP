


library(dplyr)
library(lme4)
library(ggplot2)
library(ggpubr)
library(wesanderson)



####Now for the data
npdata<-read.csv("Data/CompiledBioassay.csv")

#anova by themselves

m1<-lm(log(chl)~N*P, data=npdata[npdata$Experiment==1,])
summary(m1)
summary(aov(m1))

m2<-lm(log(chl)~N*P, data=npdata[npdata$Experiment==2,])
summary(m2)
summary(aov(m2))

m3<-lm(log(chl)~N*P, data=npdata[npdata$Experiment==3,])
summary(m3)
summary(aov(m3))

m4<-lm(log(chl)~N*P, data=npdata[npdata$Experiment==4,])
summary(m4)
summary(aov(m4))

m5<-lm(log(chl)~N*P, data=npdata[npdata$Experiment==5,])
summary(m5)
summary(aov(m5))

m6<-lm(log(chl)~N*P, data=npdata[npdata$Experiment==6,])
summary(m6)
summary(aov(m6))


mmfit<- lmer(log(chl)~N + P + N*P + (1|Experiment), data=npdata)
summary(mmfit)


coef(mmfit)$Experiment

npdata$NP<-npdata$N*npdata$P

#mmfit_re<- lmer(log(chl)~(1 | Experiment) + (1 | N:Experiment) + (1 | P:Experiment) +(1 | (NP):Experiment), data=npdata)

###plot

pdf(file="Figures/npanova.pdf", height=4.5, width = 6)

p <- ggplot(data = npdata, aes(x = N, y = chl, color = as.factor(P))) +  scale_color_manual( values=c("0"="#FF0000", "1"="#00A08A"), labels=c("0", "+P"), name="") + geom_smooth(method="lm", aes(group = P), fill=NA) + geom_point() + scale_y_continuous(trans='log10',limits=c(0.5,8))+ 
theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
panel.background = element_blank(), axis.line = element_line(colour = "black"), legend.position = "null") +
scale_x_continuous( name="", limits=c(-0.25,1.25), breaks=c(0,1), labels=c("0", "+N")) + ylab("Chlorophyll (µg/L)")

p<-p + facet_wrap(~Experiment)
p

sumNP<- data.frame(N=c(0,1,0,1), P=c(0,0,1,1), logchl=c(0.43,(0.43-0.02),(0.43 + 0.485), (0.41+0.485+0.384) ), sd=c(0.17,0.0677,0.0677,0.09574) )

propsdN<-(0.17^2+0.0677^2)^0.5
propsdNP<-(0.17^2+0.0677^2+0.0957^2)^0.5
sumNP$propsd<-c(0.17,propsdN,propsdN, propsdNP)
sumNP$chl<-exp(sumNP$logchl)
sumNP$ebup<-exp(sumNP$logchl+ sumNP$propsd)
sumNP$ebdown<-exp(sumNP$logchl-sumNP$propsd)

  
sp <- ggplot(data = sumNP, aes(x = N, y = chl, color = as.factor(P))) +  scale_color_manual( values=c("0"="#FF0000", "1"="#00A08A"), labels=c("0", "+P"), name="") + geom_line(size=1.2) + geom_pointrange(aes(ymin=ebdown, ymax=ebup))+
  scale_y_continuous(trans='log10', limits=c(0.5,8))+ theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
panel.background = element_blank(), axis.line = element_line(colour = "black"), axis.title.y = element_blank(), legend.position = c(0.5,0.8) ) +
scale_x_continuous( name="", limits=c(-0.25,1.25), breaks=c(0,1), labels=c("0", "+N")) 
sp


ggarrange(p, sp,  labels = c("A", "B") )

dev.off()

#common.legend = TRUE, legend = "bottom",








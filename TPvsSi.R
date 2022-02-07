
library(dplyr)
library(ggplot2)
library(lme4)
library (brms)

TPSi<-read.csv("Data/TPSi_Hall.csv")
TPSi$sitefac<- as.integer(as.factor(TPSi$Site))
TPSi$SiO2<-as.numeric(TPSi$SiO2) # data imported as characters
TPSi$TP<-as.numeric(TPSi$TP)
TPSi<-TPSi[TPSi$SiO2<23 & TPSi$TP<80 & !is.na(TPSi$sitefac),]
TPSi$logTP<- log(TPSi$TP)
TPSi<-TPSi[!is.na(TPSi$logTP), ]
TPSi<-TPSi[TPSi$logTP> -5, ]

TPSi<-TPSi[TPSi$sitefac !=3 , ]

plot(TPSi$SiO2,TPSi$TP)
plot(TPSi$SiO2,TPSi$TP, log="y")

ggplot(TPSi, aes(SiO2, TP)) + 
  geom_point(aes(color = Si_corrected)) +
  facet_grid(. ~ sitefac)

ggplot(TPSi, aes(SiO2, TP)) + scale_y_log10()+
  geom_point(aes(color = Si_corrected)) +
  facet_grid(. ~ sitefac)

ggplot(TPSi, aes(Si_corrected, TP)) + scale_y_log10()+
  geom_point() +
  stat_smooth(method = "lm", col = "red")+
  facet_grid(. ~ sitefac)

ggplot(TPSi, aes(Si_corrected, TP)) + 
  geom_point() +
  stat_smooth(method = "lm", col = "red")+
  facet_grid(. ~ sitefac)

##the below are not the way to go, too complicated
lmfit<- lmer(TP~ 1 + Si_corrected + (Si_corrected | sitefac), data=TPSi )
ranef(lmfit)

lmfit<- lmer(log(TP)~ 1 + Si_corrected + (Si_corrected | sitefac), data=TPSi )
ranef(lmfit)

bayesfit<- brm(formula = log(TP)~ 0 + Intercept + Si_corrected + (Si_corrected | sitefac), data = TPSi, family="gaussian", 
                     prior = c(prior(normal(0, 1), class = b),
                               prior(normal(0, 10), class = b, coef=Intercept) ))

print(bayesfit, digits=4)
ranef(bayesfit)
exp(1.91)
exp(-0.44)
exp(1.91)-exp(1.91-0.44)#2.4 is average correction
z<-as.data.frame(ranef(lmfit))
zslope<-z[z$term=="Si_corrected", ]
zint<-z[z$term=="(Intercept)", ]
predranef<-exp(1.91+zint$condval)-exp((1.91+zint$condval)+(-0.44+zslope$condval))

bayesgamfit<- brm(formula = TP~ 0 + Intercept + Si_corrected + (Si_corrected | sitefac), data = TPSi, family=Gamma(link="identity"),
               prior = c(prior(normal(0, 1), class = b),
                         prior(normal(0, 10), class = b, coef=Intercept) ), cores=4)


bayesgamlogfit<- brm(formula = TP~ 0 + Intercept + Si_corrected + (Si_corrected | sitefac), data = TPSi, family=Gamma(link="log"),
                  prior = c(prior(normal(0, 1), class = b),
                            prior(normal(0, 10), class = b, coef=Intercept) ), cores=4)




print(bayesgamfit, digits=4)
ranef(bayesgamfit)

TPSi_FL<-TPSi[TPSi$Site=="Flathead Lake, Midlake Deep",]
#TPSi_FL<-TPSi_FL[TPSi$TP<20,]  ##remove clear outliers, only 2
#TPSi_FL<-TPSi_FL[TPSi$SiO2<20,]


plot(TPSi_FL$SiO2[TPSi_FL$Si_corrected==0], TPSi_FL$TP[TPSi_FL$Si_corrected==0], 
     cex=0.8, pch=16, col="blue", ylim = c(0,12),xlim = c(2,7),  xlab="SiO2", ylab="TP", main="MLD")
points(TPSi_FL$SiO2[TPSi_FL$Si_corrected==1], TPSi_FL$TP[TPSi_FL$Si_corrected==1], 
       cex=0.8, pch=16, col="red")

TPSi_RC<-TPSi[TPSi$Site=="Roys Creek @ Springhead",]

#TPSi_FL<-TPSi[TPSi$TP<20,]  ##remove clear outliers, only 2
#TPSi_FL<-TPSi[TPSi$SiO2<20,]



plot(TPSi_RC$SiO2[TPSi_RC$Si_corrected==0], TPSi_RC$TP[TPSi_RC$Si_corrected==0], 
     cex=0.8, pch=16, col="blue", ylim = c(0,22), xlim=c(10,20), xlab="SiO2", ylab="TP", main="Roy's")
points(TPSi_RC$SiO2[TPSi_RC$Si_corrected==1], TPSi_RC$TP[TPSi_RC$Si_corrected==1], 
       cex=0.8, pch=16, col="red")



TPSi_FR<-TPSi[TPSi$Site=="Flathead River, Mainstem, @ Holt (Sportsman Bri...",]
plot(TPSi_FR$SiO2[TPSi_FR$Si_corrected==0], TPSi_FR$TP[TPSi_FR$Si_corrected==0], 
     cex=0.8, pch=16, col="blue", ylim = c(0,40), xlim=c(0,15), xlab="SiO2", ylab="TP", main="Flathead River")
points(TPSi_FR$SiO2[TPSi_FR$Si_corrected==1], TPSi_FR$TP[TPSi_FR$Si_corrected==1], 
       cex=0.8, pch=16, col="red")

TPSi_SW<-TPSi[TPSi$Site=="Swan River @ Bigfork Steel Bridge (Electric Ave...",]
plot(TPSi_SW$SiO2[TPSi_SW$Si_corrected==0], TPSi_SW$TP[TPSi_SW$Si_corrected==0], 
     cex=0.8, pch=16, col="blue", ylim = c(0,20), xlim=c(0,10), xlab="SiO2", ylab="TP", main="Swan River")
points(TPSi_SW$SiO2[TPSi_SW$Si_corrected==1], TPSi_SW$TP[TPSi_SW$Si_corrected==1], 
       cex=0.8, pch=16, col="red")

TPSi_PO<-TPSi[TPSi$Site=="Flathead River, Lower, @ Polson",]
plot(TPSi_PO$SiO2[TPSi_PO$Si_corrected==0], TPSi_PO$TP[TPSi_PO$Si_corrected==0],
     cex=0.8, pch=16, col="blue", ylim = c(1,20), xlim=c(0,10), xlab="SiO2", ylab="TP", main="Flathead, Polson")
points(TPSi_PO$SiO2[TPSi_PO$Si_corrected==1], TPSi_PO$TP[TPSi_PO$Si_corrected==1], 
       cex=0.8, pch=16, col="red")



TPSi_FL %>%
  filter( TP<20 ) %>%
  ggplot( aes(x=as.numeric(TP))) +
   geom_density(aes(color = Si_corrected))
  
  geom_density(fill="#69b3a2", color="#e9ecef")

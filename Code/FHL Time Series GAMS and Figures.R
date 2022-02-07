##### Flathead Lake Time Series Analysis- Oct 2021 #####
# GAMs used to assess trends and variability in N, P, and N:P at MLD in Flathead and Annual Nutrient loading data ##
library(tidyverse)
library(lubridate)
library(mgcv)
library(gratia)
library(readxl)
library(gridExtra)
library(ggpubr)
library(viridis)
library(gtsummary)
library(ggsci)
library(gt)
library(ggridges)
library(ggthemes)
#####Libraries#####
####Data Call ####

EpiChems<-read_csv("./Data/FL_EPI_Chems.csv")

####hOUSEKEEPING#####

EpiChems1<- EpiChems  %>%
  group_by(CollectionDate) %>%
  summarise_all(funs(first(.[!is.na(.)]))) %>%
  mutate(WYear = smwrBase::waterYear(CollectionDate, numeric=TRUE),
         Year = year(CollectionDate),
         Month = month(CollectionDate)) %>%
  filter(WYear > 1982) %>%
  filter(!is.na(CollectionDate))
   


###Flathead Lake  MLD GAM and Figure####
####TP GAM & Fig####
corr_TPGAM<-gam(corr_TPmol ~  s(NDate), method ="REML", data=EpiChems1)
summary(corr_TPGAM)
gam.check(corr_TPGAM, rep = 500)
appraise(corr_TPGAM, point_col = "steelblue", point_alpha = 0.4, method = "simulate")

maxdev<- derivatives(corr_TPGAM)%>% select(derivative, data, lower,  upper) 

maxdevDate <- maxdev %>%
  filter(derivative == min(derivative))%>%
  select(data)

EpiChems1 <- EpiChems1 %>%
  mutate(corr_TP_GAM = predict(corr_TPGAM),
         corr_TP_GAM_SE  = predict(corr_TPGAM, se.fit =T)[[2]]) 

MLD_TP<-ggplot(EpiChems1, aes(x = CollectionDate, y =corr_TPmol)) +
  geom_line(size=1.1, alpha=.9) +
  geom_ribbon(aes(ymin = corr_TP_GAM-corr_TP_GAM_SE, ymax = corr_TP_GAM+corr_TP_GAM_SE, x = CollectionDate), alpha = 0.75, inherit.aes = FALSE, fill = "gray70") +
  geom_line(aes(y = corr_TP_GAM, x=CollectionDate), inherit.aes = FALSE, size = 1.2, color = "black") +
  theme_bw() +
  theme(axis.line=element_line(size=1, lineend="square"),
        axis.text=element_text(size=12, color="black"),
        legend.position = "none",
        axis.title = element_text(size = 12, color="black"),
        plot.title=element_text(size=12, hjust=0.5, face="bold", vjust=-0.5),
        panel.background = element_rect(fill ="white", color = "black", size=1),
        panel.grid.major.y = element_blank(),
        panel.grid.minor.y = element_blank(),
        panel.grid.major.x = element_blank(),
        panel.grid.minor.x = element_blank()) +
  xlab("") +
  ylim(0, 0.4) +
  ylab(expression(paste("Total P", "  ",  "(µmol", " ", " ",  L^-1, "",")", sep=""))) 
  #ggsave("./FHL MLD and Loading Elser 2021/Plots/corr TP_GAM_Year_ms.png", device = "png", width=8, height=6, units = "in")
ggsave("FL_MLD_TP_final.pdf", device = "pdf", width=6, height=4.5, units = "in")


#### TN GAM and Fig####

TNGAM<-gam(TotalNmol ~  s(NDate), method ="REML", data=EpiChems1)
summary(TNGAM)
gam.check(TNGAM, rep = 500)
appraise(TNGAM, point_col = "steelblue", point_alpha = 0.4, method = "simulate")

maxdev<- derivatives(TNGAM)%>% select(derivative, data, lower,  upper) 

maxdevDate <- maxdev %>%
  filter(derivative == min(derivative))%>%
  select(data)

EpiChems1 <- EpiChems1 %>%
  mutate(TN_GAM = predict(TNGAM),
         TN_GAM_SE = predict(TNGAM, se.fit =T)[[2]])

MLD_TN<-ggplot(EpiChems1, aes(x = CollectionDate, y = TotalNmol)) +
  geom_line(size=1.1, alpha=.85) +
  geom_ribbon(aes(ymin = TN_GAM-TN_GAM_SE, ymax = TN_GAM+TN_GAM_SE, x = CollectionDate), alpha = 0.75, inherit.aes = FALSE, fill = "gray70") +
  geom_line(aes(y = TN_GAM, x=CollectionDate), inherit.aes = FALSE, size = 1.2, color = "black") +
  theme_bw() +
  ylim(0,12.5) +
  theme(axis.line=element_line(size=1, lineend="square"),
        axis.text=element_text(size=12, color="black"),
        legend.position = "none",
        axis.title = element_text(size = 12, color="black"),
        plot.title=element_text(size=12, hjust=0.5, face="bold", vjust=-0.5),
        panel.background = element_rect(fill ="white", color = "black", size=1),
        panel.grid.major.y = element_blank(),
        panel.grid.minor.y = element_blank(),
        panel.grid.major.x = element_blank(),
        panel.grid.minor.x = element_blank()) +
  xlab("") +
  ylab(expression(paste("Total N", "  ",  "(µmol", " ", " ",  L^-1, "",")", sep="")))

  #ggsave("FL_MLD_TN_final.png", device = "png", width=8, height=6, units = "in")
  
ggsave("./Figures/FL_MLD_TN_final.pdf", device = "pdf", width=6, height=4.5, units = "in")

#### MLD N:P Ratio GAM & Fig)####

corrNtoPGAM<-gam(corr_NtoP ~  s(NDate), method ="REML", data=EpiChems1)
summary(corrNtoPGAM)
gam.check(corrNtoPGAM, rep = 500)

appraise(corrNtoPGAM, point_col = "steelblue", point_alpha = 0.4, method = "simulate")
draw(corrNtoPGAM)
draw(corr_derivatives(corrNtoPGAM))

maxdev<- derivatives(corrNtoPGAM) %>% select(derivative, data, lower,  upper) 

maxdevDate <- maxdev %>%
  filter(derivative == min(derivative))%>%
  select(data)

EpiChems1 <- EpiChems1 %>%
  mutate(corNtoP_GAM = predict(corrNtoPGAM),
         corNtoP_GAM_SE = predict(corrNtoPGAM, se.fit =T)[[2]])

tail(EpiChems1)
Nto__predict<- predict(corrNtoPGAM)

MLD_NtoP<-ggplot(EpiChems1, aes(x = CollectionDate, y = log10(corr_NtoP))) +
  geom_line(size=1.1) +
  geom_ribbon(aes(ymin = log10(corNtoP_GAM-corNtoP_GAM_SE), ymax = log10(corNtoP_GAM + corNtoP_GAM_SE), x = CollectionDate), alpha = 0.8, inherit.aes = FALSE, fill = "gray70") +
  geom_line(aes(y = log10(corNtoP_GAM), x=CollectionDate), inherit.aes = FALSE, size = 1.2, color = "black") +
  theme_bw() +
  scale_y_continuous( breaks = c(1.5, 2, 2.5), label = c("30", "100", "325")) +
  theme(axis.line=element_line(size=1, lineend="square"),
        axis.text=element_text(size=12, color="black"),
        legend.position = "none",
        axis.title = element_text(size = 12, color="black"),
        plot.title=element_text(size=12, hjust=0.5, face="bold", vjust=-0.5),
        panel.background = element_rect(fill ="white", color = "black", size=1),
        panel.grid.major.y = element_blank(),
        panel.grid.minor.y = element_blank(),
        panel.grid.major.x = element_blank(),
        panel.grid.minor.x = element_blank()) +
  xlab("Year") +
  ylab("N:P (molar)") 
  # ggsave("./FHL MLD and Loading Elser 2021/Plots/NtoP_GAM_Yearlogged_ms.png", device = "png", width=8, height=6, units = "in")
  ggsave("FL_MLD_NtoP_final.png", device = "png", width=8, height=6, units = "in")

lay <- rbind(1, 2, 3)
MLDFig<-grid.arrange(grobs = list(MLD_TN, MLD_TP, MLD_NtoP), layout_matrix = lay)
ggsave("./Figures/MLD Figure.pdf", plot=MLDFig, device = "pdf", width=6, height= 7.5, units = "in")





#####
#### Flathead River Loading GAMS####
##### Datacall for Loading####
 
FL_DailyLoading <- read_csv("Data/FL Daily Loading via loadFlex.csv") 
Holt_Q <- read_csv("Data/LOADEST Holt Q.csv") %>%
  rename(Q = q) %>%
  mutate(date = mdy(Date)) %>%
  select(-Date)

VWC<-left_join(FL_DailyLoading, Holt_Q)



Load <-  FL_DailyLoading %>%
        mutate(NtoP_D = (TN_kg_D/14.007)/(TP_kg_D/30.973)) %>%
        mutate(WYear = smwrBase::waterYear(date, numeric=TRUE)) %>%
        left_join(Holt_Q)


YLoad<-Load %>%
  group_by(WYear) %>%
  summarize(AnnualQ = sum(Q),
  TP_MT_y= sum(TP_kg_D)/1000,
  TN_MT_y = sum(TN_kg_D)/1000) %>% 
  filter(WYear > 1982)


####TP Loading GAM ####

TPGAM<- gam(TP_MT_y ~  s(WYear) , method ="REML", data=YLoad)
summary(TPGAM)
gam.check(TPGAM, rep = 500)

appraise(TPGAM, point_col = "steelblue", point_alpha = 0.4, method = "simulate")
draw(TPGAM)
draw(derivatives(TPGAM))

YLoad <- YLoad %>%
  mutate(TP_GAM = predict(TPGAM),
         TP_GAM_SE= predict(TPGAM, se.fit =T)[[2]])
####TP Loading Figure####

LoadingTP<-ggplot(YLoad, aes(x = WYear, y = TP_MT_y)) +
  geom_line(size=1.1, color="black") +
  geom_point(size=2) +
  geom_ribbon(data = YLoad, aes(ymin = (TP_GAM-TP_GAM_SE), ymax = (TP_GAM+TP_GAM_SE), x = WYear), alpha = 0.5, inherit.aes = FALSE, fill = "gray70") +
  geom_line(data = YLoad, aes(y = TP_GAM, x = WYear), inherit.aes = FALSE, size = 1.2, color = "black") +
  theme_few() +
  xlab("") +
  ylab("Total P (Mg/y)") 
  ggsave("./FL TP loadFlex.png", device = "png", width=8, height=8, units = "in")

#### TN Loading GAM####

TNGAM<-gam(TN_MT_y ~  s(WYear), method ="REML", data=YLoad)
summary(TNGAM)
gam.check(TNGAM, rep = 500)

appraise(TNGAM, point_col = "steelblue", point_alpha = 0.4, method = "simulate")
draw(TNGAM)
draw(derivatives(TNGAM))

YLoad <- YLoad %>%
  mutate(TN_GAM = predict(TNGAM),
         TN_GAM_SE = predict(TNGAM, se.fit =T)[[2]])


LoadingTN<-ggplot(YLoad, aes(x = WYear, y = TN_MT_y)) +
  geom_line(size=1.1, color="black") +
  geom_point(size=2) +
  geom_ribbon(aes(ymin = TN_GAM-TN_GAM_SE, ymax = TN_GAM+TN_GAM_SE, x = WYear), alpha = 0.5, inherit.aes = FALSE, fill = "gray70") +
  geom_line(aes(y = TN_GAM, x = WYear), inherit.aes = FALSE, size = 1.2, color = "black") +
  theme_few() +
  xlab("") +
  ylab("Total N (Mg/y)") +
  scale_y_continuous( limits = c(0,4000)) 
  ggsave("./FL TN loadFlex.png", device = "png", width=8, height=8, units = "in")


#### TN:TP Loading GAM####  
LoadNtoPGAM<-gam((TN_MT_y/14.007)/(TP_MT_y/30.973) ~  s(WYear), method ="REML", data=YLoad)
summary(LoadNtoPGAM)
gam.check(LoadNtoPGAM, rep = 500)

appraise(LoadNtoPGAM, point_col = "steelblue", point_alpha = 0.4, method = "simulate")
draw(LoadNtoPGAM)
draw(derivatives(LoadNtoPGAM))

YLoad <- YLoad %>%
  mutate(NtoP_GAM = predict(LoadNtoPGAM),
         NtoP_GAM_SE = predict(LoadNtoPGAM, se.fit =T)[[2]])

LoadingNtoP<- ggplot(YLoad, aes(x = WYear, y =((TN_MT_y/14.007)/(TP_MT_y/30.973)))) +
  geom_line(color="black", size=1.1) +
  geom_point(size=2) +
  geom_ribbon(aes(ymin = NtoP_GAM-NtoP_GAM_SE, ymax = NtoP_GAM+NtoP_GAM_SE, x = WYear), alpha = 0.5, inherit.aes = FALSE, fill = "gray70") +
  geom_line(aes(y = NtoP_GAM, x = WYear), inherit.aes = FALSE, size = 1.2, color = "black", linetype="dashed") +
  theme_few() +
  scale_y_continuous(trans=("log10"), limits = c(3,70)) +
  xlab("Water Year") +
  ylab("Molar N:P") 

#ggsave("./Figures/NtoPloadFlex.png", device = "png", width=10, height=8, units = "in")

ggsave("./Figures/NtoPloadFlex.pdf", device = "pdf", width=6, height=4.8, units = "in")

#### Compilation Figure of Flathead River Loading ####
lay <- rbind(c(1,2),
            c(3,3))

LoadingFig<-grid.arrange(grobs = list(LoadingTN, LoadingTP, LoadingNtoP), layout_matrix = lay)
#ggsave("Loading Figure.png", plot=LoadingFig, device = "png", width=8, height= 10, units = "in")
ggsave("./Figures/Loading Figure.pdf", plot=LoadingFig, device = "pdf", width=6, height= 6, units = "in")


##mean of TN:TP in loading

YLoad$NP <- (YLoad$TN_MT_y/14.007)/(YLoad$TP_MT_y/30.973)

mean(  YLoad$NP )  

exp (mean(  log( YLoad$NP) ) ) #geometric mean


#Correlation between loading NP and annual Q
plot( YLoad$AnnualQ  ,  YLoad$NP, log='y' )


cor(YLoad$AnnualQ , log(YLoad$NP)  )
cor.test(YLoad$AnnualQ , log(YLoad$NP)  )


#### Volume weighted Concentrations###

ggplot(VWC, aes(x=date, y=TP_kg_D/(Q*86400))) +
geom_line() +
scale_y_continuous(trans=("log10")) 

mean(VWC$TP_kg_D /(VWC$Q*86400))*1e6

ggplot(VWC, aes(x=date, y=TN_kg_D/(Q*86400))) +
  geom_line() +
  scale_y_continuous(trans=("log10")) 

ggplot(VWC, aes(x=date, y=TN_kg_D/TP_kg_D) )+
  geom_line() +
  scale_y_continuous(trans=("log10")) 


ggplot(YLoad, aes(x=WYear, y=TP_MT_y/(AnnualQ))) +
  geom_line() +
  scale_y_continuous(trans=("log10")) 

ggplot(YLoad, aes(x=WYear, y=(TN_MT_y/TP_MT_y))) +
  geom_line() +
  scale_y_continuous(trans=("log10")) 



ggplot(YLoad, aes(x=WYear, y=TN_MT_y/AnnualQ)) +
  geom_line() +
  scale_y_continuous(trans=("log10"))

ggplot(YLoad, aes(x=AnnualQ, y=TN_kg_D/TP_kg_D)) +
  geom_point() +
  scale_y_continuous(trans=("log10"))

acf(YLoad$TP_MT_y/(YLoad$AnnualQ))

acf(YLoad$TN_MT_y/(YLoad$AnnualQ))

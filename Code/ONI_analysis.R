

#ONI analysis

oni<-read_csv("./Data/oni.csv")

oniDJF <- oni %>% filter (grepl("DJF", SEAS)) 


oniMAM <- oni %>% filter (grepl("MAM", SEAS)) 
oniWIN <- oniDJF
oniWIN$ANOM_MAM<-oniMAM$ANOM

oniWIN <- oniWIN %>% 
  rename(
    WYear = YR)

YLoad <-left_join(YLoad, oniWIN, by = "WYear")
YLoad$ANOMTOT <- YLoad$ANOM + YLoad$ANOM_MAM

Yload_oni <- YLoad %>% select(WYear, AnnualQ, TP_MT_y, TN_MT_y, ANOM, ANOM_MAM, ANOMTOT)

pairs(Yload_oni)


plot(Yload_oni$ANOMTOT, Yload_oni$TP_MT_y, log="y",  xlab="ONI", ylab="Tp Load")

plot(Yload_oni$ANOMTOT, Yload_oni$TN_MT_y, log="y", xlab="ONI", ylab="TN Load")

plot(Yload_oni$ANOMTOT, Yload_oni$TP_MT_y/Yload_oni$AnnualQ, log="y")
plot(Yload_oni$ANOMTOT, Yload_oni$TN_MT_y/Yload_oni$AnnualQ, log="y", xlab="ONI", ylab="TN Conc")
plot(Yload_oni$ANOMTOT, Yload_oni$AnnualQ, log="y")

summary(lm(log(TN_MT_y/AnnualQ) ~ ANOMTOT, data=Yload_oni))
summary(lm(log(TP_MT_y/AnnualQ) ~ ANOMTOT, data=Yload_oni))
summary(lm(log(TN_MT_y) ~ ANOMTOT, data=Yload_oni))
summary(lm(log(TP_MT_y) ~ ANOMTOT, data=Yload_oni))
summary(lm(log(AnnualQ) ~ ANOMTOT, data=Yload_oni))

acf(log(Yload_oni$TN_MT_y))

plot(Yload_oni$AnnualQ,(Yload_oni$TN_MT_y/Yload_oni$TP_MT_y), log="y")

Yload_oni$NP <- (Yload_oni$TN_MT_y/14.007)/(Yload_oni$TP_MT_y/30.973)
plot(Yload_oni$ANOMTOT, Yload_oni$NP, log="y")

##analyses for paper



cor.test(Yload_oni$ANOMTOT, Yload_oni$TP_MT_y)
cor.test(Yload_oni$ANOMTOT, Yload_oni$TN_MT_y)


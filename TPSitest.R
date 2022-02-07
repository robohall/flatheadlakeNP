
library(ggplot2)
library(brms)
sitest<-read.csv("Data/TPSi_test.csv")


TP_vs_Si<- ggplot(sitest, aes(SI_spike, P_conc) ) + xlab(expression(paste("SiO"[2], " spike (mg/L)" )) )  + ylab("Measured TP (µg/L)") +
  geom_point(aes(color = factor(Si_corr))) +
  geom_smooth(method = "lm", alpha = .15, aes(fill =  factor(Si_corr)))+
  facet_grid(. ~ P_spike)

pdf("/Users/bob.hall/FlatheadLake/flatheadlakeNP/Figures/TP_vs_Si.pdf", height=6, width=6)
print(TP_vs_Si)
dev.off()
#Lets remove the P spike as I see no interaction
sitest$P_corr<-sitest$P_conc-sitest$P_spike

ggplot(sitest, aes(SI_spike, P_corr)) + 
  geom_point(aes(color = factor(Si_corr))) +
  geom_smooth(method = "lm",  aes(fill =  factor(Si_corr)))



fit1 <- brm(formula = P_corr ~ SI_spike * Si_corr,
            data = sitest, 
            prior = c(set_prior("normal(0,5)", class = "b")
                      ),
                   warmup = 1000, iter = 2000, chains = 4)

summary(fit1)

fit_extract<- rstan::extract(fit1)

steps<-as.mcmc(
  fit1,
  fixed = FALSE,
  combine_chains = TRUE,
  inc_warmup = FALSE
)

corr<-NA
for ( i in 1:1000){
  
  corr[i]<- steps[i,1]+steps[i,2]* 6
}

plot(density(corr))

steps[1,2]

summary(lm(formula = P_corr ~ SI_spike * Si_corr,
           data = sitest))

fit2 <- brm(formula = bf(P_corr ~ SI_spike * Si_corr, sigma~Si_corr),
            data = sitest, 
            prior = c(set_prior("normal(0,5)", class = "b")
            ),
            warmup = 1000, iter = 2000, chains = 4)

summary(fit2)

sd_uncorr<-exp(0.41)
sdcorr<- exp(0.41-0.88)

### next fit is to just do the uncorrcted with Si as a predictor

sitest_uncorr<- sitest[sitest$Si_corr==0,]

fit3 <- brm(formula = bf(P_corr ~ SI_spike + Si_conc),
            data = sitest_uncorr, 
            prior = c(set_prior("normal(0,5)", class = "b")
            ),
            warmup = 1000, iter = 2000, chains = 4)

summary(fit3)

fit4 <- brm(formula = bf(P_corr ~ SI_spike),
            data = sitest_uncorr, 
            prior = c(set_prior("normal(0,5)", class = "b")
            ),
            warmup = 1000, iter = 2000, chains = 4)

summary(fit4)

residspike<-resid(fit4)

plot(sitest_uncorr$Si_conc, residspike[,1])
summary(lm(residspike[,1]~sitest_uncorr$Si_conc) )

fit5 <- brm(formula = bf(P_corr ~ Si_conc),
            family="student",
            data = sitest_uncorr, 
            prior = c(set_prior("normal(0,5)", class = "b")
            ),
            warmup = 1000, iter = 2000, chains = 4)
print(fit5)





pfit<-lm(sitest_uncorr$P_corr~sitest_uncorr$Si_conc) 
summary(pfit)

pdf(file="/Users/bob.hall/FlatheadLake/flatheadlakeNP/Figures/TP_corr.pdf", height=5, width=6)
plot(sitest_uncorr$Si_conc, sitest_uncorr$P_corr, xlab=expression(paste("Measured SiO"[2], " (mg/L)" ))  , ylab="P anomaly (µg P/L)", pch=16)
lines(sitest_uncorr$Si_conc, -11.25 + 0.33*sitest_uncorr$Si_conc)
dev.off()

sifit<-lm(sitest_uncorr$Si_conc~sitest_uncorr$SI_spike)
summary(sifit)

plot(sitest_uncorr$SI_spike,sitest_uncorr$Si_conc, xlab="Si spike", ylab="Measured Si")
abline(sifit)
lines(c(0,20), c(30,50), col='red')



##get correction error

steps<-as.mcmc( fit5,
  pars = c("Intercept", "Si_conc"),
  fixed = FALSE,
  combine_chains = TRUE,
  inc_warmup = FALSE
)

MLD_Si<- TPSi$SiO2[TPSi$sitefac==1]
mean(MLD_Si)
phat<- steps[1:1000,2]*(sample(MLD_Si, 1000, replace=T))

hist(phat)
quantile(phat, c(0.025,0.5,0.975))

roys_Si<- TPSi$SiO2[TPSi$sitefac==4]

phat_roys<- steps[1:1000,2]*(sample(roys_Si, 1000, replace=T))

hist(phat_roys)
quantile(phat_roys, c(0.025,0.5,0.975))

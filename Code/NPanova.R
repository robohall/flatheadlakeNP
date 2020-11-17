y ~ x1 + x2 + x1*x2 + (1|siteid)

setwd("/Users/bobhall/misc")
library(dplyr)
library(rstan)
library(lme4)
library(ggplot2)
library(wesanderson)

#some fakedata

N<-c(0,0,0,0,0,0,0,0,0,0,1,1,1,1,1,1,1,1,1,1)
P<-c(0,0,0,0,0,1,1,1,1,1,0,0,0,0,0,1,1,1,1,1)

chl<-c(5,5,5,5,5,10,10,10,10,10,5,5,5,5,5,40,40,40,40,40)

fdf<-data.frame(N=N,P=P,chl=log(chl)+rnorm(20,0,0.05))
plot(fdf)

summary(aov(chl~N*P, data=fdf))

summary(lm(chl~N*P, data=fdf))
summary(lm(chl~N+P, data=fdf))
#summary(glm(chl~N*P, data=fdf))


#aov(lm(chl~N*P, data=fdf))


group_by(fdf, N, P) %>%
  summarise(
    count = n(),
    mean = mean(chl, na.rm = TRUE),
    sd = sd(chl, na.rm = TRUE)
  )


####Now for the data
npdata<-read.csv("/Users/bobhall/misc/CompiledBioassay.csv")

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

#common.legend = TRUE, legend = "bottom",






#mmfit_re1<- lmer(log(chl)~1+ (1 | N:Experiment) + (1 | P:Experiment) +(1 | NP:Experiment), data=npdata)
#coef(mmfit_re1)








##
#setwd("/Users/bobhall/ecomodels/")

sink("2wayreg.stan")

cat("
    
    data {
    int<lower = 1> N;
    vector[N] f1;
  vector[N] f2;
    vector[N] y;
    }
    parameters {
    real b0;               
    real b1;                
     real   b2;
     real   b3;
    real<lower = 0> sigma; //  standard deviation
    }
    model {
    vector[N] yhat;
      yhat=b0 + b1 * f1 + b2 * f2 + b3 * f1 .* f2;
    y ~ normal(yhat, sigma); // likelihood
    
    b0~normal(0,5);
    b1~normal(0,1);
    b2~normal(0,1);
    b3~normal(0,1);
    
    
    }
    
   
    
    "
    ,fill=TRUE)
sink()


fdf_data<-list(N=length(fdf$chl), f1=fdf$N, f2=fdf$P, y=fdf$chl)

fit <- stan("2wayreg.stan", data = fdf_data,  iter = 2000, chains = 4)


print(fit)




sink("2wayreg_H.stan")

cat("
    
    data {
    int<lower = 1> N;
    int<lower = 1> Ngroups;
    vector[N] f1;
    vector[N] f2;
    vector[N] y;
     int <lower=0>  groupid[N];
    }

    parameters {
    vector[Ngroups] b0;               
    real b1;                
    real b2;
    real  b3;
    real<lower = 0> sigma; //  standard deviation
    real muintercept;
     real<lower = 0> intercept_sigma; //  standard deviation
    }
    model {
    vector[N] yhat;
for (i in 1:N){
    yhat[i]=b0[groupid[i]] + b1 * f1[i] + b2 * f2[i] + b3 * f1[i]* f2[i];
}
    y ~ normal(yhat, sigma); // likelihood
    
    b0~normal(muintercept,intercept_sigma);
    b1~normal(0,1);
    b2~normal(0,1);
    b3~normal(0,1);
    muintercept~normal(0,5);
  intercept_sigma~normal(0,5);
    
    }
    
    
    
    "
    ,fill=TRUE)
sink()


np_data<-list(N=length(npdata$chl), f1=npdata$N, f2=npdata$P, y=log(npdata$chl, groupid=npdata$Experiment, Ngroups=6)

fit <- stan("2wayreg_H.stan", data = np_data,  iter = 2000, chains = 4)


print(fit)



data {
  int N_obs; // number of observations
  int N_pts; // number of participants
  int K; // number of predictors + intercept
  int pid[N_obs]; // participant id vector
  matrix[N_obs, K] x; // matrix of predictors
  real y[N_obs]; // y vector
}
parameters {
  vector[K] beta_p[N_pts]; // ind participant intercept and slope coefficients by group
  vector<lower=0>[K] sigma_p; // sd for intercept and slope
  vector[K] beta; // intercept and slope hyper-priors
  corr_matrix[K] Omega; // correlation matrix
  real<lower=0> sigma; // population sigma
}
model {
  vector[N_obs] mu;
  
  // priors
  beta ~ normal(0, 1);
  Omega ~ lkj_corr(2);
  sigma_p ~ exponential(1);
  sigma ~ exponential(1);
  beta_p ~ multi_normal(beta, quad_form_diag(Omega, sigma_p));
  
  // likelihood
  for(i in 1:N_obs) {
    mu[i] = x[i] * (beta_p[pid[i]]); // * is matrix multiplication in this context
  }
  
  y ~ normal(mu, sigma);
}</lower=0></lower=0>
  
  
  stan_dat2 <- list(
    N_obs = nrow(dat),
    N_pts = max(as.numeric(dat$pid)),
    K = 2, # intercept + slope
    pid = as.numeric(dat$pid),
    x = matrix(c(rep(1, nrow(dat)), (dat$x - mean(dat$x)) / sd(dat$x)), ncol = 2), # z-score for x
    y = (dat$y - mean(dat$y)) / sd(dat$y) # z-score for y
  )



######
#Loading

  GP_dep<-read.csv("NTN-MT05-wy.csv")
  
  
  plot(GP_dep$yr,  GP_dep$totalN, type="l" )
acf(GP_dep$totalN)  

fit<-arima(GP_dep$totalN, order=c(1,0,0), xreg=GP_dep$yr)
  
  summary(fit)
  fit2<-lm(GP_dep$totalN~GP_dep$yr )
  summary(fit2)
  fit2<-gls()
  

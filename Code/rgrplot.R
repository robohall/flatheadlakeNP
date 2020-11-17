


fq<- read.csv("Data/food_quality_and_quantity.csv")
head(fq)


#model with interaction  term.  Wek evidence for interatcion
fqfit<- lm(rgr~log10(conc)*P, data=fq)
summary(fqfit)


#model with no interaction term, we went with this one.  This model assumes constan slope and variable intercept
fqfit_noint<- lm(rgr~log10(conc)+P, data=fq)
summary(fqfit_noint)

params<-coefficients(fqfit_noint)

##what if we ran an anova with these data?
aovfit<-aov(rgr~log10(conc)*P, data=fq)
summary(aovfit)


##make the figure
pdf(file="Figures/rgrplot.pdf", height=3, width = 3.5)

par(mai=c(0.7,0.7,0.1,0.1), mgp=c(2,0.8,0))
plot((fq$conc[fq$P==0]), fq$rgr[fq$P==0], ylim=c(0,0.4), xlim=c(5,60), log='x', 
     xlab=expression(paste("Seston concentration (", mu,"mol C L"^{-1},")")), ylab=expression(paste("Relative growth rate (","d"^{-1},")"))   )
points((fq$conc[fq$P==1]), fq$rgr[fq$P==1], pch=16)
lines((fq$conc[fq$P==0]), log10(10^(params[1]) * fq$conc[fq$P==0]^params[2]) )
lines((fq$conc[fq$P==1]), log10(10^(params[1]+params[3]) * fq$conc[fq$P==1]^params[2]) )

legend(x=25, y=0.1, legend=c("+P", "Control"), pch=c(16,1), bty="n")

dev.off()



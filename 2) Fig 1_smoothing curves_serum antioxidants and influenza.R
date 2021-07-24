#########################
#### Smoothing curve ####
#########################

library(foreign)
library(Hmisc)
library(survey)
library(splines)
options(survey.lonely.psu="adjust")

setwd("C:\\Users\\great\\OneDrive\\01 on going\\2020- antioxidant\\antioxidant and influenza mortality\\Github")
load("serum.rda")

attach(serum)
bpdsn45<-svydesign(id=~psu, strata=~strata, weights=~wt_mec, data=serum, nest=T)

par(mfrow=c(2,2))

#####
#vitC
svygam.vitC<-svycoxph(Surv(HSAITMOR, age_die, m_inf)~ns(s_vitC,df=4)+factor(sex)+factor(raceth)+factor(phase)
                        +factor(educ)+chol+bmi+factor(smk), bpdsn45)
summary(svygam.vitC)

fit.svygam.vitC<-predict(svygam.vitC, type="terms", se.fit=T)
fit.svygam.vitC$fit[1:10,]
fit.svygam.vitC$fit[,1] <- fit.svygam.vitC$fit[,1] +1

lo.s_vitC<-fit.svygam.vitC$fit[,1]-1.96*fit.svygam.vitC$se.fit[,1]
hi.s_vitC<-fit.svygam.vitC$fit[,1]+1.96*fit.svygam.vitC$se.fit[,1]
zero.s_vitC<-fit.svygam.vitC$fit/fit.svygam.vitC$fit
summary(lo.s_vitC)
summary(hi.s_vitC)

par(mfrow=c(2,2), new=F)
quantile(s_vitC, .995)
hist(s_vitC,  nclass=20, col = "#E6E6E6", border=F, xlim=c(0,2.56595),
     main="", axes = F, ylab="", xlab="")
o<-order(s_vitC)
par(new=T)
plot(s_vitC[o], fit.svygam.vitC$fit[,1][o], type="l",  lwd="2",
      ylim=c(0,3), xlim=c(0,2.56595),
     main="A. Vitamin C", ylab="Hazard ratio", xlab="Serum vitamin C (mg/dL)")
lines(s_vitC[o], lo.s_vitC[o], lty=2, col=1)
lines(s_vitC[o], hi.s_vitC[o], lty=2, col=1)
lines(s_vitC[o], zero.s_vitC[o], lty=1, col=8)
par(new=F)

#####
#carotene
svygam.car<-svycoxph(Surv(HSAITMOR, age_die, m_inf)~ns(I(s_car),df=4)+factor(sex)+factor(raceth)+factor(phase)
                      +factor(educ)+chol+bmi+factor(smk), bpdsn45)
summary(svygam.car)

fit.svygam.car<-predict(svygam.car, type="terms", se.fit=T)
fit.svygam.car$fit[1:10,]
fit.svygam.car$fit[,1] <- fit.svygam.car$fit[,1] +1

lo.s_car<-fit.svygam.car$fit[,1]-1.96*fit.svygam.car$se.fit[,1]
hi.s_car<-fit.svygam.car$fit[,1]+1.96*fit.svygam.car$se.fit[,1]
zero.s_car<-fit.svygam.car$fit/fit.svygam.car$fit
summary(lo.s_car)
summary(hi.s_car)

quantile(s_car, .995)
hist(s_car,  nclass=60, col = "#E6E6E6", border=F, xlim=c(0,186.865),
     main="", axes = F, ylab="", xlab="")
o<-order(s_car)
par(new=T)
plot(s_car[o], fit.svygam.car$fit[,1][o], type="l", lwd="2",
     ylim=c(0,5), xlim=c(0,186.865),
     main="B. Sum of serum メ- and モ-carotene", ylab="Hazard ratio", xlab="Sum of serum メ- and モ-carotene (レg/dL)")
lines(s_car[o], lo.s_car[o], lty=2, col=1)
lines(s_car[o], hi.s_car[o], lty=2, col=1)
lines(s_car[o], zero.s_car[o], lty=1, col=8)
par(new=F)

#####
#lycopene
svygam.lyp<-svycoxph(Surv(HSAITMOR, age_die, m_inf)~ns(I(s_lyp),df=4)+factor(sex)+factor(raceth)+factor(phase)
                     +factor(educ)+chol+bmi+factor(smk), bpdsn45)
summary(svygam.lyp)

fit.svygam.lyp<-predict(svygam.lyp, type="terms", se.fit=T)
fit.svygam.lyp$fit[1:10,]
fit.svygam.lyp$fit[,1] <- fit.svygam.lyp$fit[,1] +1

lo.s_lyp<-fit.svygam.lyp$fit[,1]-1.96*fit.svygam.lyp$se.fit[,1]
hi.s_lyp<-fit.svygam.lyp$fit[,1]+1.96*fit.svygam.lyp$se.fit[,1]
zero.s_lyp<-fit.svygam.lyp$fit/fit.svygam.lyp$fit
summary(lo.s_lyp)
summary(hi.s_lyp)

quantile(s_lyp, .995)
hist(s_lyp,  nclass=30, col = "#E6E6E6", border=F, xlim=c(0,57),
     main="", axes = F, ylab="", xlab="")
o<-order(s_lyp)
par(new=T)
plot(s_lyp[o], fit.svygam.lyp$fit[,1][o], type="l", lwd="2",
     ylim=c(-3,4), xlim=c(0,57),
     main="C. Lycopene", ylab="Hazard ratio", xlab="Serum lycopene (レg/dL)")
lines(s_lyp[o], lo.s_lyp[o], lty=2, col=1)
lines(s_lyp[o], hi.s_lyp[o], lty=2, col=1)
lines(s_lyp[o], zero.s_lyp[o], lty=1, col=8)
par(new=F)

#####
#TAC
svygam.TAC<-svycoxph(Surv(HSAITMOR, age_die, m_inf)~ns(I(s_TAC),df=4)+factor(sex)+factor(raceth)+factor(phase)
                     +factor(educ)+chol+bmi+factor(smk), bpdsn45)
summary(svygam.TAC)

fit.svygam.TAC<-predict(svygam.TAC, type="terms", se.fit=T)
fit.svygam.TAC$fit[1:10,]
fit.svygam.TAC$fit[,1] <- fit.svygam.TAC$fit[,1] +1

lo.s_TAC<-fit.svygam.TAC$fit[,1]-1.96*fit.svygam.TAC$se.fit[,1]
hi.s_TAC<-fit.svygam.TAC$fit[,1]+1.96*fit.svygam.TAC$se.fit[,1]
zero.s_TAC<-fit.svygam.TAC$fit/fit.svygam.TAC$fit
summary(lo.s_TAC)
summary(hi.s_TAC)

quantile(s_TAC, .995)
hist(s_TAC,  nclass=20, col = "#E6E6E6", border=F, xlim=c(0,3.39527),
     main="", axes = F, ylab="", xlab="")
o<-order(s_TAC)
par(new=T)
plot(s_TAC[o], fit.svygam.TAC$fit[,1][o], type="l", lwd="2",
     ylim=c(-1,4), xlim=c(0,3.39527),
     main="D. Total antioxidant capacity", ylab="Hazard ratio", xlab="Serum TAC (mg VCE/dL)")
lines(s_TAC[o], lo.s_TAC[o], lty=2, col=1)
lines(s_TAC[o], hi.s_TAC[o], lty=2, col=1)
lines(s_TAC[o], zero.s_TAC[o], lty=1, col=8)
par(new=F)


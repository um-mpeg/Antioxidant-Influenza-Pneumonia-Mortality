library(haven)
library(survival)
library(epiDisplay)
library(survey)
library(dplyr)
options(survey.lonely.psu="adjust")

setwd("C:\\Users\\great\\OneDrive\\01 on going\\2020- antioxidant\\antioxidant and influenza mortality\\Github")
load("serum.rda")
names(serum)
dim(serum)

bpdsn45<-svydesign(id=~psu, strata=~strata, weights=~wt_mec, data=serum, nest=T)

###############################
### TABLE 1 Characteristics ###
###############################

#quartile group
serum.q1 <- serum %>% filter(s_TAC.gr4=='[0.152,0.802)')
serum.q2 <- serum %>% filter(s_TAC.gr4=='[0.802,1.215)')
serum.q3 <- serum %>% filter(s_TAC.gr4=='[1.215,1.588)')
serum.q4 <- serum %>% filter(s_TAC.gr4=='[1.588,5.008]')
dim(serum.q1); dim(serum.q2); dim(serum.q3); dim(serum.q4)
bpdsn45.q1<-svydesign(id=~psu, strata=~strata, weights=~wt_mec, data=serum.q1, nest=T)
bpdsn45.q2<-svydesign(id=~psu, strata=~strata, weights=~wt_mec, data=serum.q2, nest=T)
bpdsn45.q3<-svydesign(id=~psu, strata=~strata, weights=~wt_mec, data=serum.q3, nest=T)
bpdsn45.q4<-svydesign(id=~psu, strata=~strata, weights=~wt_mec, data=serum.q4, nest=T)

# incident rate
table(serum$m_inf);table(serum.q1$m_inf);table(serum.q2$m_inf);table(serum.q3$m_inf);table(serum.q4$m_inf)

summary(serum$PERMTH_EXM/12);summary(serum.q1$PERMTH_EXM/12);summary(serum.q2$PERMTH_EXM/12);summary(serum.q3$PERMTH_EXM/12);summary(serum.q4$PERMTH_EXM/12)

ci.poisson(sum(serum$m_inf), sum(serum$PERMTH_EXM)/(1000*12), alpha=.05)
ci.poisson(sum(serum.q1$m_inf), sum(serum.q1$PERMTH_EXM)/(1000*12), alpha=.05)
ci.poisson(sum(serum.q2$m_inf), sum(serum.q2$PERMTH_EXM)/(1000*12), alpha=.05)
ci.poisson(sum(serum.q3$m_inf), sum(serum.q3$PERMTH_EXM)/(1000*12), alpha=.05)
ci.poisson(sum(serum.q4$m_inf), sum(serum.q4$PERMTH_EXM)/(1000*12), alpha=.05)

ci.poisson(sum(serum$m_inf*serum$wt_mec), sum(serum$PERMTH_EXM*serum$wt_mec)/(1000*12), alpha=.05)
ci.poisson(sum(serum.q1$m_inf*serum.q1$wt_mec), sum(serum.q1$PERMTH_EXM*serum.q1$wt_mec)/(1000*12), alpha=.05)
ci.poisson(sum(serum.q2$m_inf*serum.q2$wt_mec), sum(serum.q2$PERMTH_EXM*serum.q2$wt_mec)/(1000*12), alpha=.05)
ci.poisson(sum(serum.q3$m_inf*serum.q3$wt_mec), sum(serum.q3$PERMTH_EXM*serum.q3$wt_mec)/(1000*12), alpha=.05)
ci.poisson(sum(serum.q4$m_inf*serum.q4$wt_mec), sum(serum.q4$PERMTH_EXM*serum.q4$wt_mec)/(1000*12), alpha=.05)
names(serum)

##continuous variables
round(svymean(~age+bmi+chol,design=bpdsn45),1)
round(confint(svymean(~age+bmi+chol,design=bpdsn45)),1)
round(svyby(~age+bmi+chol, by=s_TAC.gr4, design=bpdsn45, svymean),1)
round(confint(svyby(~age+bmi+chol, by=s_TAC.gr4, design=bpdsn45, svymean)),1)

#cotinine
serum.cot <- subset(serum, is.na(cotinine)==F)
serum.q1.cot <- subset(serum.q1, is.na(cotinine)==F)
serum.q2.cot <- subset(serum.q2, is.na(cotinine)==F)
serum.q3.cot <- subset(serum.q3, is.na(cotinine)==F)
serum.q4.cot <- subset(serum.q4, is.na(cotinine)==F)
bpdsn45.cot<-svydesign(id=~psu, strata=~strata, weights=~wt_mec, data=serum.cot, nest=T)
bpdsn45.q1.cot<-svydesign(id=~psu, strata=~strata, weights=~wt_mec, data=serum.q1.cot, nest=T)
bpdsn45.q2.cot<-svydesign(id=~psu, strata=~strata, weights=~wt_mec, data=serum.q2.cot, nest=T)
bpdsn45.q3.cot<-svydesign(id=~psu, strata=~strata, weights=~wt_mec, data=serum.q3.cot, nest=T)
bpdsn45.q4.cot<-svydesign(id=~psu, strata=~strata, weights=~wt_mec, data=serum.q4.cot, nest=T)
dim(serum.cot);dim(serum.q1.cot);dim(serum.q2.cot);dim(serum.q3.cot);dim(serum.q4.cot)

cot.v<-svymean(log(serum.cot$cotinine), bpdsn45.cot)
round(c(exp(cot.v[1]), exp(cot.v[1]-1.96*0.1011), exp(cot.v[1]+1.96*0.1011)),2)
q1.cot.v<-svymean(log(serum.q1.cot$cotinine), bpdsn45.q1.cot)
round(c(exp(q1.cot.v[1]), exp(q1.cot.v[1]-1.96*0.1011), exp(q1.cot.v[1]+1.96*0.1011)),2)
q2.cot.v<-svymean(log(serum.q2.cot$cotinine), bpdsn45.q2.cot)
round(c(exp(q2.cot.v[1]), exp(q2.cot.v[1]-1.96*0.1011), exp(q2.cot.v[1]+1.96*0.1011)),2)
q3.cot.v<-svymean(log(serum.q3.cot$cotinine), bpdsn45.q3.cot)
round(c(exp(q3.cot.v[1]), exp(q3.cot.v[1]-1.96*0.1011), exp(q3.cot.v[1]+1.96*0.1011)),2)
q4.cot.v<-svymean(log(serum.q4.cot$cotinine), bpdsn45.q4.cot)
round(c(exp(q4.cot.v[1]), exp(q4.cot.v[1]-1.96*0.1011), exp(q4.cot.v[1]+1.96*0.1011)),2)

##categorical variables
svytable(~sex, Ntotal=100,design=bpdsn45)
svytable(~raceth, Ntotal=100,design=bpdsn45)
svytable(~educ, Ntotal=100,design=bpdsn45)
svytable(~smk, Ntotal=100,design=bpdsn45)

N2 <- matrix(nrow=2,ncol=4)
N2[1,] <- N2[2,] <- svytable(~s_TAC.gr4,design=bpdsn45)
N3 <- matrix(nrow=3,ncol=4)
N3[3,] <- N3[2,] <- N3[1,] <- svytable(~s_TAC.gr4,design=bpdsn45)
N4 <- matrix(nrow=4,ncol=4)
N4[3:4,] <- N4[1:2,] <- N2

round(svytable(~sex+s_TAC.gr4,design=bpdsn45)/N2*100,1)
round(svytable(~raceth+s_TAC.gr4,design=bpdsn45)/N4*100,1)
round(svytable(~educ+s_TAC.gr4,design=bpdsn45)/N3*100,1)
round(svytable(~smk+s_TAC.gr4,design=bpdsn45)/N3*100,1)

#pir
serum.pir <- serum %>% filter(pir1!=3)
dim(serum)-dim(serum.pir)
bpdsn45.pir<-svydesign(id=~psu, strata=~strata, weights=~wt_mec, data=serum.pir, nest=T)

N2 <- matrix(nrow=2,ncol=4)
N2[1,] <- N2[2,] <- svytable(~s_TAC.gr4,design=bpdsn45.pir)

round(svytable(~pir1,Ntotal=100,design=bpdsn45.pir),1)
round(svytable(~pir1+s_TAC.gr4,design=bpdsn45.pir)/N2*100,1)

#alcohol consumption
serum.alc <- serum %>% filter(is.na(alc4)==F)
dim(serum)-dim(serum.alc)
bpdsn45.alc<-svydesign(id=~psu, strata=~strata, weights=~wt_mec, data=serum.alc, nest=T)

N2 <- matrix(nrow=2,ncol=4)
N2[1,] <- N2[2,] <- svytable(~s_TAC.gr4,design=bpdsn45.alc)

round(svytable(~alc4,Ntotal=100,design=bpdsn45.alc),1)
round(svytable(~alc4+s_TAC.gr4,design=bpdsn45.alc)/N2*100,1)

#supplement use
serum.suppl <- serum %>% filter(is.na(suppl)==F)
dim(serum)-dim(serum.suppl)
bpdsn45.suppl<-svydesign(id=~psu, strata=~strata, weights=~wt_mec, data=serum.suppl, nest=T)

N2 <- matrix(nrow=2,ncol=4)
N2[1,] <- N2[2,] <- svytable(~s_TAC.gr4,design=bpdsn45.suppl)

round(svytable(~suppl,Ntotal=100,design=bpdsn45.suppl),1)
round(svytable(~suppl+s_TAC.gr4,design=bpdsn45.suppl)/N2*100,1)

#################################
### TABLE 2 serum antioxidant ###
#################################

quantile(serum$s_vitC, probs = seq(0, 1, 0.25))
quantile(serum$s_vitA, probs = seq(0, 1, 0.25))
quantile(serum$s_vitE, probs = seq(0, 1, 0.25))
quantile(serum$s_car, probs = seq(0, 1, 0.25))
quantile(serum$s_Bxp, probs = seq(0, 1, 0.25))
quantile(serum$s_lup, probs = seq(0, 1, 0.25))
quantile(serum$s_lyp, probs = seq(0, 1, 0.25))
quantile(serum$s_TAC, probs = seq(0, 1, 0.25))

#################################
### TABLE 3 Survival Analysis ###
#################################

#########################################
### use age_die instead of PERMTH_EXM ###
#########################################
#####
# deaths/total
attach(serum)
table(s_vitC.gr4, m_inf)
table(s_vitA.gr4, m_inf)
table(s_vitE.gr4, m_inf)
table(s_car.gr4, m_inf)
table(s_Bxp.gr4, m_inf)
table(s_lup.gr4, m_inf)
table(s_lyp.gr4, m_inf)
table(s_TAC.gr4, m_inf)

#####
### model 1

## vitC
# HRs and CIs
svycox.vitC.gr4<-svycoxph(Surv(HSAITMOR, age_die, m_inf)~factor(s_vitC.gr4)+factor(sex)+factor(raceth)+factor(phase)
                          , bpdsn45)
summary(svycox.vitC.gr4)
# p for trend
svycox.vitC<-svycoxph(Surv(HSAITMOR, age_die, m_inf)~as.numeric(s_vitC.gr4)+factor(sex)+factor(raceth)+factor(phase)
                     , bpdsn45)
summary(svycox.vitC)

## vitA
# HRs and CIs
svycox.vitA.gr4<-svycoxph(Surv(HSAITMOR, age_die, m_inf)~factor(s_vitA.gr4)+factor(sex)+factor(raceth)+factor(phase)
                          , bpdsn45)
summary(svycox.vitA.gr4)
# p for trend
svycox.vitA<-svycoxph(Surv(HSAITMOR, age_die, m_inf)~as.numeric(s_vitA.gr4)+factor(sex)+factor(raceth)+factor(phase)
                      , bpdsn45)
summary(svycox.vitA)

## vitE
# HRs and CIs
svycox.vitE.gr4<-svycoxph(Surv(HSAITMOR, age_die, m_inf)~factor(s_vitE.gr4)+factor(sex)+factor(raceth)+factor(phase)
                          , bpdsn45)
summary(svycox.vitE.gr4)
# p for trend
svycox.vitE<-svycoxph(Surv(HSAITMOR, age_die, m_inf)~as.numeric(s_vitE.gr4)+factor(sex)+factor(raceth)+factor(phase)
                     , bpdsn45)
summary(svycox.vitE)

## carotene
# HRs and CIs
svycox.car.gr4<-svycoxph(Surv(HSAITMOR, age_die, m_inf)~factor(s_car.gr4)+factor(sex)+factor(raceth)+factor(phase)
                          , bpdsn45)
summary(svycox.car.gr4)
# p for trend
svycox.Bcar<-svycoxph(Surv(HSAITMOR, age_die, m_inf)~as.numeric(s_car.gr4)+factor(sex)+factor(raceth)+factor(phase)
                     , bpdsn45)
summary(svycox.Bcar)

## beta cryptoxanthin
# HRs and CIs
svycox.Bxp.gr4<-svycoxph(Surv(HSAITMOR, age_die, m_inf)~factor(s_Bxp.gr4)+factor(sex)+factor(raceth)+factor(phase)
                         , bpdsn45)
summary(svycox.Bxp.gr4)
# p for trend
svycox.Bxp<-svycoxph(Surv(HSAITMOR, age_die, m_inf)~as.numeric(s_Bxp.gr4)+factor(sex)+factor(raceth)+factor(phase)
                     , bpdsn45)
summary(svycox.Bxp)

## lutein + zeaxanthin
# HRs and CIs
svycox.lup.gr4<-svycoxph(Surv(HSAITMOR, age_die, m_inf)~factor(s_lup.gr4)+factor(sex)+factor(raceth)+factor(phase)
                         , bpdsn45)
summary(svycox.lup.gr4)
# p for trend
svycox.lup<-svycoxph(Surv(HSAITMOR, age_die, m_inf)~as.numeric(s_lup.gr4)+factor(sex)+factor(raceth)+factor(phase)
                     , bpdsn45)
summary(svycox.lup)

## lycopene
# HRs and CIs
svycox.lyp.gr4<-svycoxph(Surv(HSAITMOR, age_die, m_inf)~factor(s_lyp.gr4)+factor(sex)+factor(raceth)+factor(phase)
                         , bpdsn45)
summary(svycox.lyp.gr4)
# p for trend
svycox.lyp<-svycoxph(Surv(HSAITMOR, age_die, m_inf)~as.numeric(s_lyp.gr4)+factor(sex)+factor(raceth)+factor(phase)
                     , bpdsn45)
summary(svycox.lyp)

## total antioxidant capacity
# HRs and CIs
svycox.TAC.gr4<-svycoxph(Surv(HSAITMOR, age_die, m_inf)~factor(s_TAC.gr4)+factor(sex)+factor(raceth)+factor(phase)
                         , bpdsn45)
summary(svycox.TAC.gr4)
# p for trend
svycox.TAC<-svycoxph(Surv(HSAITMOR, age_die, m_inf)~as.numeric(s_TAC.gr4)+factor(sex)+factor(raceth)+factor(phase)
                     , bpdsn45)
summary(svycox.TAC)

#####
### model 2

## vitC
# HRs and CIs
svycox.vitC.gr4.2<-svycoxph(Surv(HSAITMOR, age_die, m_inf)~factor(s_vitC.gr4)+factor(sex)+factor(raceth)+factor(phase)
                          +factor(educ)+chol+bmi+factor(smk), bpdsn45)
summary(svycox.vitC.gr4.2)
# p for trend
svycox.vitC.2<-svycoxph(Surv(HSAITMOR, age_die, m_inf)~as.numeric(s_vitC.gr4)+factor(sex)+factor(raceth)+factor(phase)
                      +factor(educ)+chol+bmi+factor(smk), bpdsn45)
summary(svycox.vitC.2)

## vitA
# HRs and CIs
svycox.vitA.gr4.2<-svycoxph(Surv(HSAITMOR, age_die, m_inf)~factor(s_vitA.gr4)+factor(sex)+factor(raceth)+factor(phase)
                          +factor(educ)+chol+bmi+factor(smk), bpdsn45)
summary(svycox.vitA.gr4.2)
# p for trend
svycox.vitA.2<-svycoxph(Surv(HSAITMOR, age_die, m_inf)~as.numeric(s_vitA.gr4)+factor(sex)+factor(raceth)+factor(phase)
                      +factor(educ)+chol+bmi+factor(smk), bpdsn45)
summary(svycox.vitA.2)

## vitE
# HRs and CIs
svycox.vitE.gr4.2<-svycoxph(Surv(HSAITMOR, age_die, m_inf)~factor(s_vitE.gr4)+factor(sex)+factor(raceth)+factor(phase)
                          +factor(educ)+chol+bmi+factor(smk), bpdsn45)
summary(svycox.vitE.gr4.2)
# p for trend
svycox.vitE.2<-svycoxph(Surv(HSAITMOR, age_die, m_inf)~as.numeric(s_vitE.gr4)+factor(sex)+factor(raceth)+factor(phase)
                      +factor(educ)+chol+bmi+factor(smk), bpdsn45)
summary(svycox.vitE.2)

## carotene
# HRs and CIs
svycox.car.gr4.2<-svycoxph(Surv(HSAITMOR, age_die, m_inf)~factor(s_car.gr4)+factor(sex)+factor(raceth)+factor(phase)
                          +factor(educ)+chol+bmi+factor(smk), bpdsn45)
summary(svycox.car.gr4.2)
# p for trend
svycox.car.2<-svycoxph(Surv(HSAITMOR, age_die, m_inf)~as.numeric(s_car.gr4)+factor(sex)+factor(raceth)+factor(phase)
                      +factor(educ)+chol+bmi+factor(smk), bpdsn45)
summary(svycox.car.2)

## beta cryptoxanthin
# HRs and CIs
svycox.Bxp.gr4.2<-svycoxph(Surv(HSAITMOR, age_die, m_inf)~factor(s_Bxp.gr4)+factor(sex)+factor(raceth)+factor(phase)
                         +factor(educ)+chol+bmi+factor(smk), bpdsn45)
summary(svycox.Bxp.gr4.2)
# p for trend
svycox.Bxp.2<-svycoxph(Surv(HSAITMOR, age_die, m_inf)~as.numeric(s_Bxp.gr4)+factor(sex)+factor(raceth)+factor(phase)
                     +factor(educ)+chol+bmi+factor(smk), bpdsn45)
summary(svycox.Bxp.2)

## lutein + zeaxanthin
# HRs and CIs
svycox.lup.gr4.2<-svycoxph(Surv(HSAITMOR, age_die, m_inf)~factor(s_lup.gr4)+factor(sex)+factor(raceth)+factor(phase)
                         +factor(educ)+chol+bmi+factor(smk), bpdsn45)
summary(svycox.lup.gr4.2)
# p for trend
svycox.lup.2<-svycoxph(Surv(HSAITMOR, age_die, m_inf)~as.numeric(s_lup.gr4)+factor(sex)+factor(raceth)+factor(phase)
                     +factor(educ)+chol+bmi+factor(smk), bpdsn45)
summary(svycox.lup.2)

## lycopene
# HRs and CIs
svycox.lyp.gr4.2<-svycoxph(Surv(HSAITMOR, age_die, m_inf)~factor(s_lyp.gr4)+factor(sex)+factor(raceth)+factor(phase)
                         +factor(educ)+chol+bmi+factor(smk), bpdsn45)
summary(svycox.lyp.gr4.2)
# p for trend
svycox.lyp.2<-svycoxph(Surv(HSAITMOR, age_die, m_inf)~as.numeric(s_lyp.gr4)+factor(sex)+factor(raceth)+factor(phase)
                     +factor(educ)+chol+bmi+factor(smk), bpdsn45)
summary(svycox.lyp.2)

## total antioxidant capacity
# HRs and CIs
svycox.TAC.gr4.2<-svycoxph(Surv(HSAITMOR, age_die, m_inf)~factor(s_TAC.gr4)+factor(sex)+factor(raceth)+factor(phase)
                         +factor(educ)+chol+bmi+factor(smk), bpdsn45)
summary(svycox.TAC.gr4.2)
# p for trend
svycox.TAC.2<-svycoxph(Surv(HSAITMOR, age_die, m_inf)~as.numeric(s_TAC.gr4)+factor(sex)+factor(raceth)+factor(phase)
                     +factor(educ)+chol+bmi+factor(smk), bpdsn45)
summary(svycox.TAC.2)


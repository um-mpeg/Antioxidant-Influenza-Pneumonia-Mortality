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

### incident rate
summ(PERMTH_EXM)
summary(PERMTH_EXM)
summary(PERMTH_EXM/12)
sum(PERMTH_EXM)
sum(wt_mec)
table(serum$m_inf)
svymean(~PERMTH_EXM, bpdsn45)
confint(svymean(~PERMTH_EXM, bpdsn45))
sum(m_inf)
1000*(sum(m_inf)/sum(PERMTH_EXM))*12
1000*(sum(m_inf*wt_mec)/sum(PERMTH_EXM*wt_mec))*12
1000*(svytotal(~m_inf, bpdsn45)/svytotal(~PERMTH_EXM, bpdsn45))*12

ci.poisson(sum(m_inf), sum(PERMTH_EXM)/(1000*12), alpha=.05)
ci.poisson(sum(m_inf*wt_mec), sum(PERMTH_EXM*wt_mec)/(1000*12), alpha=.05)

#continuous variables
svymean(~age+bmi+chol, bpdsn45)
confint(svymean(~age+bmi+chol, bpdsn45))

serum.cot <- subset(serum, is.na(cotinine)==F)
bpdsn45.cot<-svydesign(id=~psu, strata=~strata, weights=~wt_mec, data=serum.cot, nest=T)
cot.v<-svymean(log(serum.cot$cotinine), bpdsn45.cot)
exp(cot.v[1])
exp(cot.v[1]-1.96*0.1011)
exp(cot.v[1]+1.96*0.1011)

#categorical variables
svytable(~sex, Ntotal=100,design=bpdsn45)
svytable(~raceth, Ntotal=100,design=bpdsn45)
svytable(~educ, Ntotal=100,design=bpdsn45)
svytable(~pir1, Ntotal=100,design=bpdsn45)
svytable(~smk, Ntotal=100,design=bpdsn45)
htn_bp1 <- ifelse(is.na(htn_bp), 9,htn_bp)
bpdsn45<-svydesign(id=~psu, strata=~strata, weights=~wt_mec, data=serum, nest=T)
svytable(~htn_bp1, Ntotal=100,design=bpdsn45)

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


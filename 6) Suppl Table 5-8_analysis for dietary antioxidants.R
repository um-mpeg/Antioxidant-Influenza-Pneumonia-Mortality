library(haven)
library(survival)
library(epiDisplay)
library(survey)
library(dplyr)
options(survey.lonely.psu="adjust")

setwd("C:\\Users\\great\\OneDrive\\01 on going\\2020- antioxidant\\antioxidant and influenza mortality\\Github")
load("diet.rda")
names(diet)
dim(diet)

bpdsn45<-svydesign(id=~psu, strata=~strata, weights=~wt_mec, data=diet, nest=T)

############################
### Supplemental Table 5 ###
############################

### incident rate
summ(PERMTH_EXM)
summary(PERMTH_EXM)
summary(PERMTH_EXM/12)
sum(PERMTH_EXM)
sum(wt_mec)
table(diet$m_inf)
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
svymean(~d_vitC+d_reti+d_vitE+d_Bcar+d_TAC, bpdsn45)
confint(svymean(~d_vitC+d_reti+d_vitE+d_Bcar+d_TAC, bpdsn45))
svyquantile(~d_vitC+d_reti+d_vitE+d_Bcar+d_TAC, bpdsn45, quantiles=c(.25, .5, .75))

diet.cot <- subset(diet, cotinine > 0)
bpdsn45.cot<-svydesign(id=~psu, strata=~strata, weights=~wt_mec, data=diet.cot, nest=T)
cot.v<-svymean(log(diet.cot$cotinine), bpdsn45.cot)
exp(cot.v[1])
exp(cot.v[1]-1.96*0.1023)
exp(cot.v[1]+1.96*0.1023)


#categorical variables
svytable(~sex, Ntotal=100,design=bpdsn45)
svytable(~raceth, Ntotal=100,design=bpdsn45)
svytable(~educ, Ntotal=100,design=bpdsn45)
svytable(~pir1, Ntotal=100,design=bpdsn45)
svytable(~smk, Ntotal=100,design=bpdsn45)

############################
### Supplemental Table 6 ###
############################

quantile(diet$a_vitC, probs = seq(0, 1, 0.25))
quantile(diet$a_reti, probs = seq(0, 1, 0.25))
quantile(diet$a_vitE, probs = seq(0, 1, 0.25))
quantile(diet$a_Bcar, probs = seq(0, 1, 0.25))

############################
### Supplemental Table 7 ###
############################

#########################################
### use age_die instead of PERMTH_EXM ###
#########################################

## deaths/total
attach(diet)
table(a_vitC.gr4, m_inf)
table(a_reti.gr4, m_inf)
table(a_vitE.gr4, m_inf)
table(a_Bcar.gr4, m_inf)
table(a_vitC.gr4, m_inf)[,1]+table(a_vitC.gr4, m_inf)[,2]
table(a_reti.gr4, m_inf)[,1]+table(a_reti.gr4, m_inf)[,2]
table(a_vitE.gr4, m_inf)[,1]+table(a_vitE.gr4, m_inf)[,2]
table(a_Bcar.gr4, m_inf)[,1]+table(a_Bcar.gr4, m_inf)[,2]
names(diet)
#####
## model 1

## vitC
# HRs and CIs
svycox.vitC.gr4<-svycoxph(Surv(HSAITMOR, age_die, m_inf)~factor(a_vitC.gr4)+factor(sex)+factor(raceth)+factor(phase)+lnenergy
                          , bpdsn45)
summary(svycox.vitC.gr4)
# p for trend
svycox.vitC<-svycoxph(Surv(HSAITMOR, age_die, m_inf)~as.numeric(a_vitC.gr4)+factor(sex)+factor(raceth)+factor(phase)+lnenergy
                     , bpdsn45)
summary(svycox.vitC)

## reti
# HRs and CIs
svycox.reti.gr4<-svycoxph(Surv(HSAITMOR, age_die, m_inf)~factor(a_reti.gr4)+factor(sex)+factor(raceth)+factor(phase)+lnenergy
                          , bpdsn45)
summary(svycox.reti.gr4)
# p for trend
svycox.reti<-svycoxph(Surv(HSAITMOR, age_die, m_inf)~as.numeric(a_reti.gr4)+factor(sex)+factor(raceth)+factor(phase)+lnenergy
                      , bpdsn45)
summary(svycox.reti)

## vitE
# HRs and CIs
svycox.vitE.gr4<-svycoxph(Surv(HSAITMOR, age_die, m_inf)~factor(a_vitE.gr4)+factor(sex)+factor(raceth)+factor(phase)+lnenergy
                          , bpdsn45)
summary(svycox.vitE.gr4)
# p for trend
svycox.vitE<-svycoxph(Surv(HSAITMOR, age_die, m_inf)~as.numeric(a_vitE.gr4)+factor(sex)+factor(raceth)+factor(phase)+lnenergy
                     , bpdsn45)
summary(svycox.vitE)

## beta carotene
# HRs and CIs
svycox.Bcar.gr4<-svycoxph(Surv(HSAITMOR, age_die, m_inf)~factor(a_Bcar.gr4)+factor(sex)+factor(raceth)+factor(phase)+lnenergy
                          , bpdsn45)
summary(svycox.Bcar.gr4)
# p for trend
svycox.Bcar<-svycoxph(Surv(HSAITMOR, age_die, m_inf)~as.numeric(a_Bcar.gr4)+factor(sex)+factor(raceth)+factor(phase)+lnenergy
                     , bpdsn45)
summary(svycox.Bcar)


#####
## model 2

## vitC
# HRs and CIs
svycox.vitC.gr4.2<-svycoxph(Surv(HSAITMOR, age_die, m_inf)~factor(a_vitC.gr4)+factor(sex)+factor(raceth)+factor(phase)+lnenergy
                            +factor(educ)+chol+bmi+factor(smk), bpdsn45)
summary(svycox.vitC.gr4.2)
# p for trend
svycox.vitC.2<-svycoxph(Surv(HSAITMOR, age_die, m_inf)~as.numeric(a_vitC.gr4)+factor(sex)+factor(raceth)+factor(phase)+lnenergy
                        +factor(educ)+chol+bmi+factor(smk), bpdsn45)
summary(svycox.vitC.2)

## reti
# HRs and CIs
svycox.reti.gr4.2<-svycoxph(Surv(HSAITMOR, age_die, m_inf)~factor(a_reti.gr4)+factor(sex)+factor(raceth)+factor(phase)+lnenergy
                            +factor(educ)+chol+bmi+factor(smk), bpdsn45)
summary(svycox.reti.gr4.2)
# p for trend
svycox.reti.2<-svycoxph(Surv(HSAITMOR, age_die, m_inf)~as.numeric(a_reti.gr4)+factor(sex)+factor(raceth)+factor(phase)+lnenergy
                        +factor(educ)+chol+bmi+factor(smk), bpdsn45)
summary(svycox.reti.2)

## vitE
# HRs and CIs
svycox.vitE.gr4.2<-svycoxph(Surv(HSAITMOR, age_die, m_inf)~factor(a_vitE.gr4)+factor(sex)+factor(raceth)+factor(phase)+lnenergy
                            +factor(educ)+chol+bmi+factor(smk), bpdsn45)
summary(svycox.vitE.gr4.2)
# p for trend
svycox.vitE.2<-svycoxph(Surv(HSAITMOR, age_die, m_inf)~as.numeric(a_vitE.gr4)+factor(sex)+factor(raceth)+factor(phase)+lnenergy
                        +factor(educ)+chol+bmi+factor(smk), bpdsn45)
summary(svycox.vitE.2)

## beta carotene
# HRs and CIs
svycox.Bcar.gr4.2<-svycoxph(Surv(HSAITMOR, age_die, m_inf)~factor(a_Bcar.gr4)+factor(sex)+factor(raceth)+factor(phase)+lnenergy
                            +factor(educ)+chol+bmi+factor(smk), bpdsn45)
summary(svycox.Bcar.gr4.2)
# p for trend
svycox.Bcar.2<-svycoxph(Surv(HSAITMOR, age_die, m_inf)~as.numeric(a_Bcar.gr4)+factor(sex)+factor(raceth)+factor(phase)+lnenergy
                        +factor(educ)+chol+bmi+factor(smk), bpdsn45)
summary(svycox.Bcar.2)
############################
### Supplemental Table 8 ###
############################

#####
# +vitamin C

#reti
svycox.reti.gr4<-svycoxph(Surv(HSAITMOR, age_die, m_inf)~factor(a_reti.gr4)+factor(sex)+factor(raceth)+factor(phase)+lnenergy
                          +factor(educ)+chol+bmi+factor(smk)+factor(a_vitC.gr4), bpdsn45)
summary(svycox.reti.gr4)
svycox.reti.gr4<-svycoxph(Surv(HSAITMOR, age_die, m_inf)~as.numeric(a_reti.gr4)+factor(sex)+factor(raceth)+factor(phase)+lnenergy
                          +factor(educ)+chol+bmi+factor(smk)+factor(a_vitC.gr4), bpdsn45)
summary(svycox.reti.gr4)

#vitE
svycox.vitE.gr4<-svycoxph(Surv(HSAITMOR, age_die, m_inf)~factor(a_vitE.gr4)+factor(sex)+factor(raceth)+factor(phase)+lnenergy
                          +factor(educ)+chol+bmi+factor(smk)+factor(a_vitC.gr4), bpdsn45)
summary(svycox.vitE.gr4)
svycox.vitE.gr4<-svycoxph(Surv(HSAITMOR, age_die, m_inf)~as.numeric(a_vitE.gr4)+factor(sex)+factor(raceth)+factor(phase)+lnenergy
                          +factor(educ)+chol+bmi+factor(smk)+factor(a_vitC.gr4), bpdsn45)
summary(svycox.vitE.gr4)

#beta carotene
svycox.Bcar.gr4<-svycoxph(Surv(HSAITMOR, age_die, m_inf)~factor(a_Bcar.gr4)+factor(sex)+factor(raceth)+factor(phase)+lnenergy
                          +factor(educ)+chol+bmi+factor(smk)+factor(a_vitC.gr4), bpdsn45)
summary(svycox.Bcar.gr4)
svycox.Bcar.gr4<-svycoxph(Surv(HSAITMOR, age_die, m_inf)~as.numeric(a_Bcar.gr4)+factor(sex)+factor(raceth)+factor(phase)+lnenergy
                          +factor(educ)+chol+bmi+factor(smk)+factor(a_vitC.gr4), bpdsn45)
summary(svycox.Bcar.gr4)

#####
# +vitamin E

#vitC
svycox.vitC.gr4<-svycoxph(Surv(HSAITMOR, age_die, m_inf)~factor(a_vitC.gr4)+factor(sex)+factor(raceth)+factor(phase)+lnenergy
                          +factor(educ)+chol+bmi+factor(smk)+factor(a_vitE.gr4), bpdsn45)
summary(svycox.vitC.gr4)
svycox.vitC.gr4<-svycoxph(Surv(HSAITMOR, age_die, m_inf)~as.numeric(a_vitC.gr4)+factor(sex)+factor(raceth)+factor(phase)+lnenergy
                          +factor(educ)+chol+bmi+factor(smk)+factor(a_vitE.gr4), bpdsn45)
summary(svycox.vitC.gr4)

#reti
svycox.reti.gr4<-svycoxph(Surv(HSAITMOR, age_die, m_inf)~factor(a_reti.gr4)+factor(sex)+factor(raceth)+factor(phase)+lnenergy
                          +factor(educ)+chol+bmi+factor(smk)+factor(a_vitE.gr4), bpdsn45)
summary(svycox.reti.gr4)
svycox.reti.gr4<-svycoxph(Surv(HSAITMOR, age_die, m_inf)~as.numeric(a_reti.gr4)+factor(sex)+factor(raceth)+factor(phase)+lnenergy
                          +factor(educ)+chol+bmi+factor(smk)+factor(a_vitE.gr4), bpdsn45)
summary(svycox.reti.gr4)

#beta carotene
svycox.Bcar.gr4<-svycoxph(Surv(HSAITMOR, age_die, m_inf)~factor(a_Bcar.gr4)+factor(sex)+factor(raceth)+factor(phase)+lnenergy
                          +factor(educ)+chol+bmi+factor(smk)+factor(a_vitE.gr4), bpdsn45)
summary(svycox.Bcar.gr4)
svycox.Bcar.gr4<-svycoxph(Surv(HSAITMOR, age_die, m_inf)~as.numeric(a_Bcar.gr4)+factor(sex)+factor(raceth)+factor(phase)+lnenergy
                          +factor(educ)+chol+bmi+factor(smk)+factor(a_vitE.gr4), bpdsn45)
summary(svycox.Bcar.gr4)


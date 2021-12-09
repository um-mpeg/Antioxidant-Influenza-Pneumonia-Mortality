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
### Supplemental Table 7 ###
############################

#########################################
### use age_die instead of PERMTH_EXM ###
#########################################

## deaths/total
attach(diet)
table(a_vitC.gr4, m_inf)
table(a_reti.gr4, m_inf)
table(a_alpha.gr4, m_inf)
table(a_gamma.gr4, m_inf)
table(a_Bcar.gr4, m_inf)
table(a_vitC.gr4, m_inf)[,1]+table(a_vitC.gr4, m_inf)[,2]
table(a_reti.gr4, m_inf)[,1]+table(a_reti.gr4, m_inf)[,2]
table(a_alpha.gr4, m_inf)[,1]+table(a_alpha.gr4, m_inf)[,2]
table(a_gamma.gr4, m_inf)[,1]+table(a_gamma.gr4, m_inf)[,2]
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

## alpha
# HRs and CIs
svycox.alpha.gr4<-svycoxph(Surv(HSAITMOR, age_die, m_inf)~factor(a_alpha.gr4)+factor(sex)+factor(raceth)+factor(phase)+lnenergy
                          , bpdsn45)
summary(svycox.alpha.gr4)
# p for trend
svycox.alpha<-svycoxph(Surv(HSAITMOR, age_die, m_inf)~as.numeric(a_alpha.gr4)+factor(sex)+factor(raceth)+factor(phase)+lnenergy
                     , bpdsn45)
summary(svycox.alpha)

## gamma-tocopherol
# HRs and CIs
svycox.gamma.gr4<-svycoxph(Surv(HSAITMOR, age_die, m_inf)~factor(a_gamma.gr4)+factor(sex)+factor(raceth)+factor(phase)+lnenergy
                          , bpdsn45)
summary(svycox.gamma.gr4)
# p for trend
svycox.gamma<-svycoxph(Surv(HSAITMOR, age_die, m_inf)~as.numeric(a_gamma.gr4)+factor(sex)+factor(raceth)+factor(phase)+lnenergy
                      , bpdsn45)
summary(svycox.gamma)

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

## alpha-tocopherol
# HRs and CIs
svycox.alpha.gr4.2<-svycoxph(Surv(HSAITMOR, age_die, m_inf)~factor(a_alpha.gr4)+factor(sex)+factor(raceth)+factor(phase)+lnenergy
                            +factor(educ)+chol+bmi+factor(smk), bpdsn45)
summary(svycox.alpha.gr4.2)
# p for trend
svycox.alpha.2<-svycoxph(Surv(HSAITMOR, age_die, m_inf)~as.numeric(a_alpha.gr4)+factor(sex)+factor(raceth)+factor(phase)+lnenergy
                        +factor(educ)+chol+bmi+factor(smk), bpdsn45)
summary(svycox.alpha.2)

## gamma-tocopherol
# HRs and CIs
svycox.gamma.gr4.2<-svycoxph(Surv(HSAITMOR, age_die, m_inf)~factor(a_gamma.gr4)+factor(sex)+factor(raceth)+factor(phase)+lnenergy
                             +factor(educ)+chol+bmi+factor(smk), bpdsn45)
summary(svycox.gamma.gr4.2)
# p for trend
svycox.gamma.2<-svycoxph(Surv(HSAITMOR, age_die, m_inf)~as.numeric(a_gamma.gr4)+factor(sex)+factor(raceth)+factor(phase)+lnenergy
                         +factor(educ)+chol+bmi+factor(smk), bpdsn45)
summary(svycox.gamma.2)

## beta carotene
# HRs and CIs
svycox.Bcar.gr4.2<-svycoxph(Surv(HSAITMOR, age_die, m_inf)~factor(a_Bcar.gr4)+factor(sex)+factor(raceth)+factor(phase)+lnenergy
                            +factor(educ)+chol+bmi+factor(smk), bpdsn45)
summary(svycox.Bcar.gr4.2)
# p for trend
svycox.Bcar.2<-svycoxph(Surv(HSAITMOR, age_die, m_inf)~as.numeric(a_Bcar.gr4)+factor(sex)+factor(raceth)+factor(phase)+lnenergy
                        +factor(educ)+chol+bmi+factor(smk), bpdsn45)
summary(svycox.Bcar.2)


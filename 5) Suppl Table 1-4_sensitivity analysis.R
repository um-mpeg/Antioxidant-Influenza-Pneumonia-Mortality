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


############################
### Supplemental Table 1 ###
############################

## deaths/total
attach(serum)
table(s_vitC.gr, m_inf)
table(s_vitA.gr, m_inf)
table(s_vitE.gr, m_inf)

## model 1
svycox.vitC.gr<-svycoxph(Surv(HSAITMOR, age_die, m_inf)~factor(s_vitC.gr)+factor(sex)+factor(raceth)+factor(phase)
                         , bpdsn45)
summary(svycox.vitC.gr)

svycox.vitA.gr<-svycoxph(Surv(HSAITMOR, age_die, m_inf)~factor(s_vitA.gr)+factor(sex)+factor(raceth)+factor(phase)
                         , bpdsn45)
summary(svycox.vitA.gr)

svycox.vitE.gr<-svycoxph(Surv(HSAITMOR, age_die, m_inf)~factor(s_vitE.gr)+factor(sex)+factor(raceth)+factor(phase)
                         , bpdsn45)
summary(svycox.vitE.gr)

## model 2
svycox.vitC.gr.2<-svycoxph(Surv(HSAITMOR, age_die, m_inf)~factor(s_vitC.gr)+factor(sex)+factor(raceth)+factor(phase)
                           +factor(educ)+chol+bmi+factor(smk), bpdsn45)
summary(svycox.vitC.gr.2)

svycox.vitA.gr.2<-svycoxph(Surv(HSAITMOR, age_die, m_inf)~factor(s_vitA.gr)+factor(sex)+factor(raceth)+factor(phase)
                           +factor(educ)+chol+bmi+factor(smk), bpdsn45)
summary(svycox.vitA.gr.2)

svycox.vitE.gr.2<-svycoxph(Surv(HSAITMOR, age_die, m_inf)~factor(s_vitE.gr)+factor(sex)+factor(raceth)+factor(phase)
                           +factor(educ)+chol+bmi+factor(smk), bpdsn45)
summary(svycox.vitE.gr.2)



############################
### Supplemental Table 2 ###
############################
#####
## +hypertension
serum.hpt<-subset(serum, is.na(htn_bp)==F)
bpdsn.hpt.45<-svydesign(id=~psu, strata=~strata, weights=~wt_mec, data=serum.hpt, nest=T)

# vitC
svycox.vitC.gr4.2<-svycoxph(Surv(HSAITMOR, age_die, m_inf)~factor(s_vitC.gr4)+factor(sex)+factor(raceth)+factor(phase)
                            +factor(educ)+chol+bmi+factor(smk)+factor(htn_bp), bpdsn.hpt.45)
summary(svycox.vitC.gr4.2)

# vitA
svycox.vitA.gr4.2<-svycoxph(Surv(HSAITMOR, age_die, m_inf)~factor(s_vitA.gr4)+factor(sex)+factor(raceth)+factor(phase)
                            +factor(educ)+chol+bmi+factor(smk)+factor(htn_bp), bpdsn.hpt.45)
summary(svycox.vitA.gr4.2)

# vitE
svycox.vitE.gr4.2<-svycoxph(Surv(HSAITMOR, age_die, m_inf)~factor(s_vitE.gr4)+factor(sex)+factor(raceth)+factor(phase)
                            +factor(educ)+chol+bmi+factor(smk)+factor(htn_bp), bpdsn.hpt.45)
summary(svycox.vitE.gr4.2)

# carotene
svycox.car.gr4.2<-svycoxph(Surv(HSAITMOR, age_die, m_inf)~factor(s_car.gr4)+factor(sex)+factor(raceth)+factor(phase)
                           +factor(educ)+chol+bmi+factor(smk)+factor(htn_bp), bpdsn.hpt.45)
summary(svycox.car.gr4.2)

# beta cryptoxanthin
svycox.Bxp.gr4.2<-svycoxph(Surv(HSAITMOR, age_die, m_inf)~factor(s_Bxp.gr4)+factor(sex)+factor(raceth)+factor(phase)
                           +factor(educ)+chol+bmi+factor(smk)+factor(htn_bp), bpdsn.hpt.45)
summary(svycox.Bxp.gr4.2)

# lutein + zeaxanthin
svycox.lup.gr4.2<-svycoxph(Surv(HSAITMOR, age_die, m_inf)~factor(s_lup.gr4)+factor(sex)+factor(raceth)+factor(phase)
                           +factor(educ)+chol+bmi+factor(smk)+factor(htn_bp), bpdsn.hpt.45)
summary(svycox.lup.gr4.2)

# lycopene
svycox.lyp.gr4.2<-svycoxph(Surv(HSAITMOR, age_die, m_inf)~factor(s_lyp.gr4)+factor(sex)+factor(raceth)+factor(phase)
                           +factor(educ)+chol+bmi+factor(smk)+factor(htn_bp), bpdsn.hpt.45)
summary(svycox.lyp.gr4.2)

# total antioxidant capacity
svycox.TAC.gr4.2<-svycoxph(Surv(HSAITMOR, age_die, m_inf)~factor(s_TAC.gr4)+factor(sex)+factor(raceth)+factor(phase)
                           +factor(educ)+chol+bmi+factor(smk)+factor(htn_bp), bpdsn.hpt.45)
summary(svycox.TAC.gr4.2)

#####
## +serum cotinine
serum.cot<-subset(serum, is.na(cotinine)==FALSE)
bpdsn.cot.45<-svydesign(id=~psu, strata=~strata, weights=~wt_mec, data=serum.cot, nest=T)

# vitC
svycox.vitC.gr4.2<-svycoxph(Surv(HSAITMOR, age_die, m_inf)~factor(s_vitC.gr4)+factor(sex)+factor(raceth)+factor(phase)
                            +factor(educ)+chol+bmi+factor(smk)+I(log(cotinine)), bpdsn.cot.45)
summary(svycox.vitC.gr4.2)

# vitA
svycox.vitA.gr4.2<-svycoxph(Surv(HSAITMOR, age_die, m_inf)~factor(s_vitA.gr4)+factor(sex)+factor(raceth)+factor(phase)
                            +factor(educ)+chol+bmi+factor(smk)+I(log(cotinine)), bpdsn.cot.45)
summary(svycox.vitA.gr4.2)

# vitE
svycox.vitE.gr4.2<-svycoxph(Surv(HSAITMOR, age_die, m_inf)~factor(s_vitE.gr4)+factor(sex)+factor(raceth)+factor(phase)
                            +factor(educ)+chol+bmi+factor(smk)+I(log(cotinine)), bpdsn.cot.45)
summary(svycox.vitE.gr4.2)

# carotene
svycox.car.gr4.2<-svycoxph(Surv(HSAITMOR, age_die, m_inf)~factor(s_car.gr4)+factor(sex)+factor(raceth)+factor(phase)
                           +factor(educ)+chol+bmi+factor(smk)+I(log(cotinine)), bpdsn.cot.45)
summary(svycox.car.gr4.2)

# beta cryptoxanthin
svycox.Bxp.gr4.2<-svycoxph(Surv(HSAITMOR, age_die, m_inf)~factor(s_Bxp.gr4)+factor(sex)+factor(raceth)+factor(phase)
                           +factor(educ)+chol+bmi+factor(smk)+I(log(cotinine)), bpdsn.cot.45)
summary(svycox.Bxp.gr4.2)

# lutein + zeaxanthin
svycox.lup.gr4.2<-svycoxph(Surv(HSAITMOR, age_die, m_inf)~factor(s_lup.gr4)+factor(sex)+factor(raceth)+factor(phase)
                           +factor(educ)+chol+bmi+factor(smk)+I(log(cotinine)), bpdsn.cot.45)
summary(svycox.lup.gr4.2)

# lycopene
svycox.lyp.gr4.2<-svycoxph(Surv(HSAITMOR, age_die, m_inf)~factor(s_lyp.gr4)+factor(sex)+factor(raceth)+factor(phase)
                           +factor(educ)+chol+bmi+factor(smk)+I(log(cotinine)), bpdsn.cot.45)
summary(svycox.lyp.gr4.2)

# total antioxidant capacity
svycox.TAC.gr4.2<-svycoxph(Surv(HSAITMOR, age_die, m_inf)~factor(s_TAC.gr4)+factor(sex)+factor(raceth)+factor(phase)
                           +factor(educ)+chol+bmi+factor(smk)+I(log(cotinine)), bpdsn.cot.45)
summary(svycox.TAC.gr4.2)

#####
## +PIR
serum.pir<-subset(serum, is.na(pir1)==F)
bpdsn.pir.45<-svydesign(id=~psu, strata=~strata, weights=~wt_mec, data=serum.pir, nest=T)

# vitC
svycox.vitC.gr4.2<-svycoxph(Surv(HSAITMOR, age_die, m_inf)~factor(s_vitC.gr4)+factor(sex)+factor(raceth)+factor(phase)
                            +factor(educ)+chol+bmi+factor(smk)+pir, bpdsn.pir.45)
summary(svycox.vitC.gr4.2)

# vitA
svycox.vitA.gr4.2<-svycoxph(Surv(HSAITMOR, age_die, m_inf)~factor(s_vitA.gr4)+factor(sex)+factor(raceth)+factor(phase)
                            +factor(educ)+chol+bmi+factor(smk)+pir, bpdsn.pir.45)
summary(svycox.vitA.gr4.2)

# vitE
svycox.vitE.gr4.2<-svycoxph(Surv(HSAITMOR, age_die, m_inf)~factor(s_vitE.gr4)+factor(sex)+factor(raceth)+factor(phase)
                            +factor(educ)+chol+bmi+factor(smk)+pir, bpdsn.pir.45)
summary(svycox.vitE.gr4.2)

# carotene
svycox.car.gr4.2<-svycoxph(Surv(HSAITMOR, age_die, m_inf)~factor(s_car.gr4)+factor(sex)+factor(raceth)+factor(phase)
                           +factor(educ)+chol+bmi+factor(smk)+pir, bpdsn.pir.45)
summary(svycox.car.gr4.2)

# beta cryptoxanthin
svycox.Bxp.gr4.2<-svycoxph(Surv(HSAITMOR, age_die, m_inf)~factor(s_Bxp.gr4)+factor(sex)+factor(raceth)+factor(phase)
                           +factor(educ)+chol+bmi+factor(smk)+pir, bpdsn.pir.45)
summary(svycox.Bxp.gr4.2)

# lutein + zeaxanthin
svycox.lup.gr4.2<-svycoxph(Surv(HSAITMOR, age_die, m_inf)~factor(s_lup.gr4)+factor(sex)+factor(raceth)+factor(phase)
                           +factor(educ)+chol+bmi+factor(smk)+pir, bpdsn.pir.45)
summary(svycox.lup.gr4.2)

# lycopene
svycox.lyp.gr4.2<-svycoxph(Surv(HSAITMOR, age_die, m_inf)~factor(s_lyp.gr4)+factor(sex)+factor(raceth)+factor(phase)
                           +factor(educ)+chol+bmi+factor(smk)+pir, bpdsn.pir.45)
summary(svycox.lyp.gr4.2)

# total antioxidant capacity
svycox.TAC.gr4.2<-svycoxph(Surv(HSAITMOR, age_die, m_inf)~factor(s_TAC.gr4)+factor(sex)+factor(raceth)+factor(phase)
                           +factor(educ)+chol+bmi+factor(smk)+pir, bpdsn.pir.45)
summary(svycox.TAC.gr4.2)

############################
### Supplemental Table 3 ###
############################

#####
## +vitamin C

# vitA
svycox.vitA.gr4.2<-svycoxph(Surv(HSAITMOR, age_die, m_inf)~factor(s_vitA.gr4)+factor(sex)+factor(raceth)+factor(phase)
                            +factor(educ)+chol+bmi+factor(smk)+factor(s_vitC.gr4), bpdsn45)
summary(svycox.vitA.gr4.2)

# vitE
svycox.vitE.gr4.2<-svycoxph(Surv(HSAITMOR, age_die, m_inf)~factor(s_vitE.gr4)+factor(sex)+factor(raceth)+factor(phase)
                            +factor(educ)+chol+bmi+factor(smk)+factor(s_vitC.gr4), bpdsn45)
summary(svycox.vitE.gr4.2)

# carotene
svycox.car.gr4.2<-svycoxph(Surv(HSAITMOR, age_die, m_inf)~factor(s_car.gr4)+factor(sex)+factor(raceth)+factor(phase)
                           +factor(educ)+chol+bmi+factor(smk)+factor(s_vitC.gr4), bpdsn45)
summary(svycox.car.gr4.2)

# beta cryptoxanthin
svycox.Bxp.gr4.2<-svycoxph(Surv(HSAITMOR, age_die, m_inf)~factor(s_Bxp.gr4)+factor(sex)+factor(raceth)+factor(phase)
                           +factor(educ)+chol+bmi+factor(smk)+factor(s_vitC.gr4), bpdsn45)
summary(svycox.Bxp.gr4.2)

# lutein + zeaxanthin
svycox.lup.gr4.2<-svycoxph(Surv(HSAITMOR, age_die, m_inf)~factor(s_lup.gr4)+factor(sex)+factor(raceth)+factor(phase)
                           +factor(educ)+chol+bmi+factor(smk)+factor(s_vitC.gr4), bpdsn45)
summary(svycox.lup.gr4.2)

# lycopene
svycox.lyp.gr4.2<-svycoxph(Surv(HSAITMOR, age_die, m_inf)~factor(s_lyp.gr4)+factor(sex)+factor(raceth)+factor(phase)
                           +factor(educ)+chol+bmi+factor(smk)+factor(s_vitC.gr4), bpdsn45)
summary(svycox.lyp.gr4.2)

#####
## +carotene

# vitC
svycox.vitA.gr4.2<-svycoxph(Surv(HSAITMOR, age_die, m_inf)~factor(s_vitC.gr4)+factor(sex)+factor(raceth)+factor(phase)
                            +factor(educ)+chol+bmi+factor(smk)+factor(s_car.gr4), bpdsn45)
summary(svycox.vitA.gr4.2)

# vitA
svycox.vitA.gr4.2<-svycoxph(Surv(HSAITMOR, age_die, m_inf)~factor(s_vitA.gr4)+factor(sex)+factor(raceth)+factor(phase)
                            +factor(educ)+chol+bmi+factor(smk)+factor(s_car.gr4), bpdsn45)
summary(svycox.vitA.gr4.2)

# vitE
svycox.vitE.gr4.2<-svycoxph(Surv(HSAITMOR, age_die, m_inf)~factor(s_vitE.gr4)+factor(sex)+factor(raceth)+factor(phase)
                            +factor(educ)+chol+bmi+factor(smk)+factor(s_car.gr4), bpdsn45)
summary(svycox.vitE.gr4.2)

# beta cryptoxanthin
svycox.Bxp.gr4.2<-svycoxph(Surv(HSAITMOR, age_die, m_inf)~factor(s_Bxp.gr4)+factor(sex)+factor(raceth)+factor(phase)
                           +factor(educ)+chol+bmi+factor(smk)+factor(s_car.gr4), bpdsn45)
summary(svycox.Bxp.gr4.2)

#l utein + zeaxanthin
svycox.lup.gr4.2<-svycoxph(Surv(HSAITMOR, age_die, m_inf)~factor(s_lup.gr4)+factor(sex)+factor(raceth)+factor(phase)
                           +factor(educ)+chol+bmi+factor(smk)+factor(s_car.gr4), bpdsn45)
summary(svycox.lup.gr4.2)

# lycopene
svycox.lyp.gr4.2<-svycoxph(Surv(HSAITMOR, age_die, m_inf)~factor(s_lyp.gr4)+factor(sex)+factor(raceth)+factor(phase)
                           +factor(educ)+chol+bmi+factor(smk)+factor(s_car.gr4), bpdsn45)
summary(svycox.lyp.gr4.2)

#####
## +lycopene

# vitC
svycox.vitA.gr4.2<-svycoxph(Surv(HSAITMOR, age_die, m_inf)~factor(s_vitC.gr4)+factor(sex)+factor(raceth)+factor(phase)
                            +factor(educ)+chol+bmi+factor(smk)+factor(s_lyp.gr4), bpdsn45)
summary(svycox.vitA.gr4.2)

# vitA
svycox.vitA.gr4.2<-svycoxph(Surv(HSAITMOR, age_die, m_inf)~factor(s_vitA.gr4)+factor(sex)+factor(raceth)+factor(phase)
                            +factor(educ)+chol+bmi+factor(smk)+factor(s_lyp.gr4), bpdsn45)
summary(svycox.vitA.gr4.2)

# vitE
svycox.vitE.gr4.2<-svycoxph(Surv(HSAITMOR, age_die, m_inf)~factor(s_vitE.gr4)+factor(sex)+factor(raceth)+factor(phase)
                            +factor(educ)+chol+bmi+factor(smk)+factor(s_lyp.gr4), bpdsn45)
summary(svycox.vitE.gr4.2)

# carotene
svycox.vitE.gr4.2<-svycoxph(Surv(HSAITMOR, age_die, m_inf)~factor(s_car.gr4)+factor(sex)+factor(raceth)+factor(phase)
                            +factor(educ)+chol+bmi+factor(smk)+factor(s_lyp.gr4), bpdsn45)
summary(svycox.vitE.gr4.2)

# beta cryptoxanthin
svycox.Bxp.gr4.2<-svycoxph(Surv(HSAITMOR, age_die, m_inf)~factor(s_Bxp.gr4)+factor(sex)+factor(raceth)+factor(phase)
                           +factor(educ)+chol+bmi+factor(smk)+factor(s_lyp.gr4), bpdsn45)
summary(svycox.Bxp.gr4.2)

# lutein + zeaxanthin
svycox.lup.gr4.2<-svycoxph(Surv(HSAITMOR, age_die, m_inf)~factor(s_lup.gr4)+factor(sex)+factor(raceth)+factor(phase)
                           +factor(educ)+chol+bmi+factor(smk)+factor(s_lyp.gr4), bpdsn45)
summary(svycox.lup.gr4.2)

############################
### Supplemental Table 4 ###
############################

### Stratification by follow time < 10
serum.short<-subset(serum, PERMTH_EXM/12<10)
bpdsn.short.45<-svydesign(id=~psu, strata=~strata, weights=~wt_mec, data=serum.short, nest=T)

attach(serum.short)
table(m_inf)
table(s_vitC.gr4, m_inf)
table(s_car.gr4, m_inf)
table(s_lyp.gr4, m_inf)
table(s_TAC.gr4, m_inf)

#vitC
svycox.vitC.gr4.2<-svycoxph(Surv(HSAITMOR, age_die, m_inf)~factor(s_vitC.gr4)+factor(sex)+factor(raceth)+factor(phase)
                            +factor(educ)+chol+bmi+factor(smk), bpdsn.short.45)
summary(svycox.vitC.gr4.2)

#carotene
svycox.car.gr4.2<-svycoxph(Surv(HSAITMOR, age_die, m_inf)~factor(s_car.gr4)+factor(sex)+factor(raceth)+factor(phase)
                           +factor(educ)+chol+bmi+factor(smk), bpdsn.short.45)
summary(svycox.car.gr4.2)

#lycopene
svycox.lyp.gr4.2<-svycoxph(Surv(HSAITMOR, age_die, m_inf)~factor(s_lyp.gr4)+factor(sex)+factor(raceth)+factor(phase)
                           +factor(educ)+chol+bmi+factor(smk), bpdsn.short.45)
summary(svycox.lyp.gr4.2)

#total antioxidant capacity
svycox.TAC.gr4.2<-svycoxph(Surv(HSAITMOR, age_die, m_inf)~factor(s_TAC.gr4)+factor(sex)+factor(raceth)+factor(phase)
                           +factor(educ)+chol+bmi+factor(smk), bpdsn.short.45)
summary(svycox.TAC.gr4.2)

### Stratification by age <85
serum.young<-subset(serum, age<85)
bpdsn.young.45<-svydesign(id=~psu, strata=~strata, weights=~wt_mec, data=serum.young, nest=T)

attach(serum.young)
table(m_inf)
table(s_vitC.gr4, m_inf)
table(s_car.gr4, m_inf)
table(s_lyp.gr4, m_inf)
table(s_TAC.gr4, m_inf)

#vitC
svycox.vitC.gr4.2<-svycoxph(Surv(HSAITMOR, age_die, m_inf)~factor(s_vitC.gr4)+factor(sex)+factor(raceth)+factor(phase)
                            +factor(educ)+chol+bmi+factor(smk), bpdsn.young.45)
summary(svycox.vitC.gr4.2)

#carotene
svycox.car.gr4.2<-svycoxph(Surv(HSAITMOR, age_die, m_inf)~factor(s_car.gr4)+factor(sex)+factor(raceth)+factor(phase)
                           +factor(educ)+chol+bmi+factor(smk), bpdsn.young.45)
summary(svycox.car.gr4.2)

#lycopene
svycox.lyp.gr4.2<-svycoxph(Surv(HSAITMOR, age_die, m_inf)~factor(s_lyp.gr4)+factor(sex)+factor(raceth)+factor(phase)
                           +factor(educ)+chol+bmi+factor(smk), bpdsn.young.45)
summary(svycox.lyp.gr4.2)

#total antioxidant capacity
svycox.TAC.gr4.2<-svycoxph(Surv(HSAITMOR, age_die, m_inf)~factor(s_TAC.gr4)+factor(sex)+factor(raceth)+factor(phase)
                           +factor(educ)+chol+bmi+factor(smk), bpdsn.young.45)
summary(svycox.TAC.gr4.2)

### Stratification by follow time < 10 and age < 85
serum.short.old<-subset(serum, PERMTH_EXM/12<10 & age < 85)
bpdsn.short.old.45<-svydesign(id=~psu, strata=~strata, weights=~wt_mec, data=serum.short.old, nest=T)

attach(serum.short.old)
table(m_inf)
table(s_vitC.gr4, m_inf)
table(s_car.gr4, m_inf)
table(s_lyp.gr4, m_inf)
table(s_TAC.gr4, m_inf)

#vitC
svycox.vitC.gr4.2<-svycoxph(Surv(HSAITMOR, age_die, m_inf)~factor(s_vitC.gr4)+factor(sex)+factor(raceth)+factor(phase)
                            +factor(educ)+chol+bmi+factor(smk), bpdsn.short.old.45)
summary(svycox.vitC.gr4.2)

#carotene
svycox.car.gr4.2<-svycoxph(Surv(HSAITMOR, age_die, m_inf)~factor(s_car.gr4)+factor(sex)+factor(raceth)+factor(phase)
                           +factor(educ)+chol+bmi+factor(smk), bpdsn.short.old.45)
summary(svycox.car.gr4.2)

#lycopene
svycox.lyp.gr4.2<-svycoxph(Surv(HSAITMOR, age_die, m_inf)~factor(s_lyp.gr4)+factor(sex)+factor(raceth)+factor(phase)
                           +factor(educ)+chol+bmi+factor(smk), bpdsn.short.old.45)
summary(svycox.lyp.gr4.2)

#total antioxidant capacity
svycox.TAC.gr4.2<-svycoxph(Surv(HSAITMOR, age_die, m_inf)~factor(s_TAC.gr4)+factor(sex)+factor(raceth)+factor(phase)
                           +factor(educ)+chol+bmi+factor(smk), bpdsn.short.old.45)
summary(svycox.TAC.gr4.2)


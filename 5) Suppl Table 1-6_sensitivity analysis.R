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


#####
## +alcohol consumption
serum.alc<-subset(serum, is.na(alc4)==F)
bpdsn.alc.45<-svydesign(id=~psu, strata=~strata, weights=~wt_mec, data=serum.alc, nest=T)

# vitC
svycox.vitC.gr4.2<-svycoxph(Surv(HSAITMOR, age_die, m_inf)~factor(s_vitC.gr4)+factor(sex)+factor(raceth)+factor(phase)
                            +factor(educ)+chol+bmi+factor(smk)+factor(alc4), bpdsn.alc.45)
summary(svycox.vitC.gr4.2)

# vitA
svycox.vitA.gr4.2<-svycoxph(Surv(HSAITMOR, age_die, m_inf)~factor(s_vitA.gr4)+factor(sex)+factor(raceth)+factor(phase)
                            +factor(educ)+chol+bmi+factor(smk)+factor(alc4), bpdsn.alc.45)
summary(svycox.vitA.gr4.2)

# vitE
svycox.vitE.gr4.2<-svycoxph(Surv(HSAITMOR, age_die, m_inf)~factor(s_vitE.gr4)+factor(sex)+factor(raceth)+factor(phase)
                            +factor(educ)+chol+bmi+factor(smk)+factor(alc4), bpdsn.alc.45)
summary(svycox.vitE.gr4.2)

# carotene
svycox.car.gr4.2<-svycoxph(Surv(HSAITMOR, age_die, m_inf)~factor(s_car.gr4)+factor(sex)+factor(raceth)+factor(phase)
                           +factor(educ)+chol+bmi+factor(smk)+factor(alc4), bpdsn.alc.45)
summary(svycox.car.gr4.2)

# beta cryptoxanthin
svycox.Bxp.gr4.2<-svycoxph(Surv(HSAITMOR, age_die, m_inf)~factor(s_Bxp.gr4)+factor(sex)+factor(raceth)+factor(phase)
                           +factor(educ)+chol+bmi+factor(smk)+factor(alc4), bpdsn.alc.45)
summary(svycox.Bxp.gr4.2)

# lutein + zeaxanthin
svycox.lup.gr4.2<-svycoxph(Surv(HSAITMOR, age_die, m_inf)~factor(s_lup.gr4)+factor(sex)+factor(raceth)+factor(phase)
                           +factor(educ)+chol+bmi+factor(smk)+factor(alc4), bpdsn.alc.45)
summary(svycox.lup.gr4.2)

# lycopene
svycox.lyp.gr4.2<-svycoxph(Surv(HSAITMOR, age_die, m_inf)~factor(s_lyp.gr4)+factor(sex)+factor(raceth)+factor(phase)
                           +factor(educ)+chol+bmi+factor(smk)+factor(alc4), bpdsn.alc.45)
summary(svycox.lyp.gr4.2)

# total antioxidant capacity
svycox.TAC.gr4.2<-svycoxph(Surv(HSAITMOR, age_die, m_inf)~factor(s_TAC.gr4)+factor(sex)+factor(raceth)+factor(phase)
                           +factor(educ)+chol+bmi+factor(smk)+factor(alc4), bpdsn.alc.45)
summary(svycox.TAC.gr4.2)



#####
## +supplement use
serum.suppl<-subset(serum, is.na(suppl)==F)
bpdsn.suppl.45<-svydesign(id=~psu, strata=~strata, weights=~wt_mec, data=serum.suppl, nest=T)

# vitC
svycox.vitC.gr4.2<-svycoxph(Surv(HSAITMOR, age_die, m_inf)~factor(s_vitC.gr4)+factor(sex)+factor(raceth)+factor(phase)
                            +factor(educ)+chol+bmi+factor(smk)+factor(suppl), bpdsn.suppl.45)
summary(svycox.vitC.gr4.2)

# vitA
svycox.vitA.gr4.2<-svycoxph(Surv(HSAITMOR, age_die, m_inf)~factor(s_vitA.gr4)+factor(sex)+factor(raceth)+factor(phase)
                            +factor(educ)+chol+bmi+factor(smk)+factor(suppl), bpdsn.suppl.45)
summary(svycox.vitA.gr4.2)

# vitE
svycox.vitE.gr4.2<-svycoxph(Surv(HSAITMOR, age_die, m_inf)~factor(s_vitE.gr4)+factor(sex)+factor(raceth)+factor(phase)
                            +factor(educ)+chol+bmi+factor(smk)+factor(suppl), bpdsn.suppl.45)
summary(svycox.vitE.gr4.2)

# carotene
svycox.car.gr4.2<-svycoxph(Surv(HSAITMOR, age_die, m_inf)~factor(s_car.gr4)+factor(sex)+factor(raceth)+factor(phase)
                           +factor(educ)+chol+bmi+factor(smk)+factor(suppl), bpdsn.suppl.45)
summary(svycox.car.gr4.2)

# beta cryptoxanthin
svycox.Bxp.gr4.2<-svycoxph(Surv(HSAITMOR, age_die, m_inf)~factor(s_Bxp.gr4)+factor(sex)+factor(raceth)+factor(phase)
                           +factor(educ)+chol+bmi+factor(smk)+factor(suppl), bpdsn.suppl.45)
summary(svycox.Bxp.gr4.2)

# lutein + zeaxanthin
svycox.lup.gr4.2<-svycoxph(Surv(HSAITMOR, age_die, m_inf)~factor(s_lup.gr4)+factor(sex)+factor(raceth)+factor(phase)
                           +factor(educ)+chol+bmi+factor(smk)+factor(suppl), bpdsn.suppl.45)
summary(svycox.lup.gr4.2)

# lycopene
svycox.lyp.gr4.2<-svycoxph(Surv(HSAITMOR, age_die, m_inf)~factor(s_lyp.gr4)+factor(sex)+factor(raceth)+factor(phase)
                           +factor(educ)+chol+bmi+factor(smk)+factor(suppl), bpdsn.suppl.45)
summary(svycox.lyp.gr4.2)

# total antioxidant capacity
svycox.TAC.gr4.2<-svycoxph(Surv(HSAITMOR, age_die, m_inf)~factor(s_TAC.gr4)+factor(sex)+factor(raceth)+factor(phase)
                           +factor(educ)+chol+bmi+factor(smk)+factor(suppl), bpdsn.suppl.45)
summary(svycox.TAC.gr4.2)


#####
## +crp_g use
serum.crp_g<-subset(serum, is.na(crp_g)==F)
bpdsn.crp_g.45<-svydesign(id=~psu, strata=~strata, weights=~wt_mec, data=serum.crp_g, nest=T)

# vitC
svycox.vitC.gr4.2<-svycoxph(Surv(HSAITMOR, age_die, m_inf)~factor(s_vitC.gr4)+factor(sex)+factor(raceth)+factor(phase)
                            +factor(educ)+chol+bmi+factor(smk)+factor(crp_g), bpdsn.crp_g.45)
summary(svycox.vitC.gr4.2)

# vitA
svycox.vitA.gr4.2<-svycoxph(Surv(HSAITMOR, age_die, m_inf)~factor(s_vitA.gr4)+factor(sex)+factor(raceth)+factor(phase)
                            +factor(educ)+chol+bmi+factor(smk)+factor(crp_g), bpdsn.crp_g.45)
summary(svycox.vitA.gr4.2)

# vitE
svycox.vitE.gr4.2<-svycoxph(Surv(HSAITMOR, age_die, m_inf)~factor(s_vitE.gr4)+factor(sex)+factor(raceth)+factor(phase)
                            +factor(educ)+chol+bmi+factor(smk)+factor(crp_g), bpdsn.crp_g.45)
summary(svycox.vitE.gr4.2)

# carotene
svycox.car.gr4.2<-svycoxph(Surv(HSAITMOR, age_die, m_inf)~factor(s_car.gr4)+factor(sex)+factor(raceth)+factor(phase)
                           +factor(educ)+chol+bmi+factor(smk)+factor(crp_g), bpdsn.crp_g.45)
summary(svycox.car.gr4.2)

# beta cryptoxanthin
svycox.Bxp.gr4.2<-svycoxph(Surv(HSAITMOR, age_die, m_inf)~factor(s_Bxp.gr4)+factor(sex)+factor(raceth)+factor(phase)
                           +factor(educ)+chol+bmi+factor(smk)+factor(crp_g), bpdsn.crp_g.45)
summary(svycox.Bxp.gr4.2)

# lutein + zeaxanthin
svycox.lup.gr4.2<-svycoxph(Surv(HSAITMOR, age_die, m_inf)~factor(s_lup.gr4)+factor(sex)+factor(raceth)+factor(phase)
                           +factor(educ)+chol+bmi+factor(smk)+factor(crp_g), bpdsn.crp_g.45)
summary(svycox.lup.gr4.2)

# lycopene
svycox.lyp.gr4.2<-svycoxph(Surv(HSAITMOR, age_die, m_inf)~factor(s_lyp.gr4)+factor(sex)+factor(raceth)+factor(phase)
                           +factor(educ)+chol+bmi+factor(smk)+factor(crp_g), bpdsn.crp_g.45)
summary(svycox.lyp.gr4.2)

# total antioxidant capacity
svycox.TAC.gr4.2<-svycoxph(Surv(HSAITMOR, age_die, m_inf)~factor(s_TAC.gr4)+factor(sex)+factor(raceth)+factor(phase)
                           +factor(educ)+chol+bmi+factor(smk)+factor(crp_g), bpdsn.crp_g.45)
summary(svycox.TAC.gr4.2)


############################
### Supplemental Table 3 ###
############################
N2 <- matrix(nrow=2,ncol=4)
N2[1,] <- N2[2,] <- svytable(~s_TAC.gr4,design=bpdsn45)

serum.htn <- serum %>% filter(is.na(htn_bp)==F)
dim(serum)-dim(serum.htn)
bpdsn45.htn<-svydesign(id=~psu, strata=~strata, weights=~wt_mec, data=serum.htn, nest=T)
round(svytable(~htn_bp+s_TAC.gr4,design=bpdsn45.htn)/N2*100,1)

serum.dm <- serum %>% filter(is.na(diag_dm)==F)
dim(serum)-dim(serum.dm)
bpdsn45.dm<-svydesign(id=~psu, strata=~strata, weights=~wt_mec, data=serum.dm, nest=T)
round(svytable(~diag_dm+s_TAC.gr4,design=bpdsn45.dm)/N2*100,1)

serum.heart <- serum %>% filter(is.na(diag_heart_attack)==F)
dim(serum)-dim(serum.heart)
bpdsn45.heart<-svydesign(id=~psu, strata=~strata, weights=~wt_mec, data=serum.heart, nest=T)
round(svytable(~diag_heart_attack+s_TAC.gr4,design=bpdsn45.heart)/N2*100,1)

serum.COPD <- serum %>% filter(is.na(diag_COPD)==F)
dim(serum)-dim(serum.COPD)
bpdsn45.COPD<-svydesign(id=~psu, strata=~strata, weights=~wt_mec, data=serum.COPD, nest=T)
round(svytable(~diag_COPD+s_TAC.gr4,design=bpdsn45.COPD)/N2*100,1)

serum.ca <- serum %>% filter(is.na(diag_ca)==F)
dim(serum)-dim(serum.ca)
bpdsn45.ca<-svydesign(id=~psu, strata=~strata, weights=~wt_mec, data=serum.ca, nest=T)
round(svytable(~diag_ca+s_TAC.gr4,design=bpdsn45.ca)/N2*100,1)


############################
### Supplemental Table 4 ###
############################
#####
## +hypertension
# vitC
svycox.vitC.gr4.2<-svycoxph(Surv(HSAITMOR, age_die, m_inf)~factor(s_vitC.gr4)+factor(sex)+factor(raceth)+factor(phase)
                            +factor(educ)+chol+bmi+factor(smk)+factor(htn_bp), bpdsn45.htn)
summary(svycox.vitC.gr4.2)

# vitA
svycox.vitA.gr4.2<-svycoxph(Surv(HSAITMOR, age_die, m_inf)~factor(s_vitA.gr4)+factor(sex)+factor(raceth)+factor(phase)
                            +factor(educ)+chol+bmi+factor(smk)+factor(htn_bp), bpdsn45.htn)
summary(svycox.vitA.gr4.2)

# vitE
svycox.vitE.gr4.2<-svycoxph(Surv(HSAITMOR, age_die, m_inf)~factor(s_vitE.gr4)+factor(sex)+factor(raceth)+factor(phase)
                            +factor(educ)+chol+bmi+factor(smk)+factor(htn_bp), bpdsn45.htn)
summary(svycox.vitE.gr4.2)

# carotene
svycox.car.gr4.2<-svycoxph(Surv(HSAITMOR, age_die, m_inf)~factor(s_car.gr4)+factor(sex)+factor(raceth)+factor(phase)
                           +factor(educ)+chol+bmi+factor(smk)+factor(htn_bp), bpdsn45.htn)
summary(svycox.car.gr4.2)

# beta cryptoxanthin
svycox.Bxp.gr4.2<-svycoxph(Surv(HSAITMOR, age_die, m_inf)~factor(s_Bxp.gr4)+factor(sex)+factor(raceth)+factor(phase)
                           +factor(educ)+chol+bmi+factor(smk)+factor(htn_bp), bpdsn45.htn)
summary(svycox.Bxp.gr4.2)

# lutein + zeaxanthin
svycox.lup.gr4.2<-svycoxph(Surv(HSAITMOR, age_die, m_inf)~factor(s_lup.gr4)+factor(sex)+factor(raceth)+factor(phase)
                           +factor(educ)+chol+bmi+factor(smk)+factor(htn_bp), bpdsn45.htn)
summary(svycox.lup.gr4.2)

# lycopene
svycox.lyp.gr4.2<-svycoxph(Surv(HSAITMOR, age_die, m_inf)~factor(s_lyp.gr4)+factor(sex)+factor(raceth)+factor(phase)
                           +factor(educ)+chol+bmi+factor(smk)+factor(htn_bp), bpdsn45.htn)
summary(svycox.lyp.gr4.2)

# total antioxidant capacity
svycox.TAC.gr4.2<-svycoxph(Surv(HSAITMOR, age_die, m_inf)~factor(s_TAC.gr4)+factor(sex)+factor(raceth)+factor(phase)
                           +factor(educ)+chol+bmi+factor(smk)+factor(htn_bp), bpdsn45.htn)
summary(svycox.TAC.gr4.2)

#####
## +diabetes
# vitC
svycox.vitC.gr4.2<-svycoxph(Surv(HSAITMOR, age_die, m_inf)~factor(s_vitC.gr4)+factor(sex)+factor(raceth)+factor(phase)
                            +factor(educ)+chol+bmi+factor(smk)+factor(diag_dm), bpdsn45.dm)
summary(svycox.vitC.gr4.2)

# vitA
svycox.vitA.gr4.2<-svycoxph(Surv(HSAITMOR, age_die, m_inf)~factor(s_vitA.gr4)+factor(sex)+factor(raceth)+factor(phase)
                            +factor(educ)+chol+bmi+factor(smk)+factor(diag_dm), bpdsn45.dm)
summary(svycox.vitA.gr4.2)

# vitE
svycox.vitE.gr4.2<-svycoxph(Surv(HSAITMOR, age_die, m_inf)~factor(s_vitE.gr4)+factor(sex)+factor(raceth)+factor(phase)
                            +factor(educ)+chol+bmi+factor(smk)+factor(diag_dm), bpdsn45.dm)
summary(svycox.vitE.gr4.2)

# carotene
svycox.car.gr4.2<-svycoxph(Surv(HSAITMOR, age_die, m_inf)~factor(s_car.gr4)+factor(sex)+factor(raceth)+factor(phase)
                           +factor(educ)+chol+bmi+factor(smk)+factor(diag_dm), bpdsn45.dm)
summary(svycox.car.gr4.2)

# beta cryptoxanthin
svycox.Bxp.gr4.2<-svycoxph(Surv(HSAITMOR, age_die, m_inf)~factor(s_Bxp.gr4)+factor(sex)+factor(raceth)+factor(phase)
                           +factor(educ)+chol+bmi+factor(smk)+factor(diag_dm), bpdsn45.dm)
summary(svycox.Bxp.gr4.2)

# lutein + zeaxanthin
svycox.lup.gr4.2<-svycoxph(Surv(HSAITMOR, age_die, m_inf)~factor(s_lup.gr4)+factor(sex)+factor(raceth)+factor(phase)
                           +factor(educ)+chol+bmi+factor(smk)+factor(diag_dm), bpdsn45.dm)
summary(svycox.lup.gr4.2)

# lycopene
svycox.lyp.gr4.2<-svycoxph(Surv(HSAITMOR, age_die, m_inf)~factor(s_lyp.gr4)+factor(sex)+factor(raceth)+factor(phase)
                           +factor(educ)+chol+bmi+factor(smk)+factor(diag_dm), bpdsn45.dm)
summary(svycox.lyp.gr4.2)

# total antioxidant capacity
svycox.TAC.gr4.2<-svycoxph(Surv(HSAITMOR, age_die, m_inf)~factor(s_TAC.gr4)+factor(sex)+factor(raceth)+factor(phase)
                           +factor(educ)+chol+bmi+factor(smk)+factor(diag_dm), bpdsn45.dm)
summary(svycox.TAC.gr4.2)

#####
## +heart attack
# vitC
svycox.vitC.gr4.2<-svycoxph(Surv(HSAITMOR, age_die, m_inf)~factor(s_vitC.gr4)+factor(sex)+factor(raceth)+factor(phase)
                            +factor(educ)+chol+bmi+factor(smk)+factor(diag_heart_attack_attack), bpdsn45.heart)
summary(svycox.vitC.gr4.2)

# vitA
svycox.vitA.gr4.2<-svycoxph(Surv(HSAITMOR, age_die, m_inf)~factor(s_vitA.gr4)+factor(sex)+factor(raceth)+factor(phase)
                            +factor(educ)+chol+bmi+factor(smk)+factor(diag_heart_attack), bpdsn45.heart)
summary(svycox.vitA.gr4.2)

# vitE
svycox.vitE.gr4.2<-svycoxph(Surv(HSAITMOR, age_die, m_inf)~factor(s_vitE.gr4)+factor(sex)+factor(raceth)+factor(phase)
                            +factor(educ)+chol+bmi+factor(smk)+factor(diag_heart_attack), bpdsn45.heart)
summary(svycox.vitE.gr4.2)

# carotene
svycox.car.gr4.2<-svycoxph(Surv(HSAITMOR, age_die, m_inf)~factor(s_car.gr4)+factor(sex)+factor(raceth)+factor(phase)
                           +factor(educ)+chol+bmi+factor(smk)+factor(diag_heart_attack), bpdsn45.heart)
summary(svycox.car.gr4.2)

# beta cryptoxanthin
svycox.Bxp.gr4.2<-svycoxph(Surv(HSAITMOR, age_die, m_inf)~factor(s_Bxp.gr4)+factor(sex)+factor(raceth)+factor(phase)
                           +factor(educ)+chol+bmi+factor(smk)+factor(diag_heart_attack), bpdsn45.heart)
summary(svycox.Bxp.gr4.2)

# lutein + zeaxanthin
svycox.lup.gr4.2<-svycoxph(Surv(HSAITMOR, age_die, m_inf)~factor(s_lup.gr4)+factor(sex)+factor(raceth)+factor(phase)
                           +factor(educ)+chol+bmi+factor(smk)+factor(diag_heart_attack), bpdsn45.heart)
summary(svycox.lup.gr4.2)

# lycopene
svycox.lyp.gr4.2<-svycoxph(Surv(HSAITMOR, age_die, m_inf)~factor(s_lyp.gr4)+factor(sex)+factor(raceth)+factor(phase)
                           +factor(educ)+chol+bmi+factor(smk)+factor(diag_heart_attack), bpdsn45.heart)
summary(svycox.lyp.gr4.2)

# total antioxidant capacity
svycox.TAC.gr4.2<-svycoxph(Surv(HSAITMOR, age_die, m_inf)~factor(s_TAC.gr4)+factor(sex)+factor(raceth)+factor(phase)
                           +factor(educ)+chol+bmi+factor(smk)+factor(diag_heart_attack), bpdsn45.heart)
summary(svycox.TAC.gr4.2)

#####
## +COPD
# vitC
svycox.vitC.gr4.2<-svycoxph(Surv(HSAITMOR, age_die, m_inf)~factor(s_vitC.gr4)+factor(sex)+factor(raceth)+factor(phase)
                            +factor(educ)+chol+bmi+factor(smk)+factor(diag_COPD), bpdsn45.COPD)
summary(svycox.vitC.gr4.2)

# vitA
svycox.vitA.gr4.2<-svycoxph(Surv(HSAITMOR, age_die, m_inf)~factor(s_vitA.gr4)+factor(sex)+factor(raceth)+factor(phase)
                            +factor(educ)+chol+bmi+factor(smk)+factor(diag_COPD), bpdsn45.COPD)
summary(svycox.vitA.gr4.2)

# vitE
svycox.vitE.gr4.2<-svycoxph(Surv(HSAITMOR, age_die, m_inf)~factor(s_vitE.gr4)+factor(sex)+factor(raceth)+factor(phase)
                            +factor(educ)+chol+bmi+factor(smk)+factor(diag_COPD), bpdsn45.COPD)
summary(svycox.vitE.gr4.2)

# carotene
svycox.car.gr4.2<-svycoxph(Surv(HSAITMOR, age_die, m_inf)~factor(s_car.gr4)+factor(sex)+factor(raceth)+factor(phase)
                           +factor(educ)+chol+bmi+factor(smk)+factor(diag_COPD), bpdsn45.COPD)
summary(svycox.car.gr4.2)

# beta cryptoxanthin
svycox.Bxp.gr4.2<-svycoxph(Surv(HSAITMOR, age_die, m_inf)~factor(s_Bxp.gr4)+factor(sex)+factor(raceth)+factor(phase)
                           +factor(educ)+chol+bmi+factor(smk)+factor(diag_COPD), bpdsn45.COPD)
summary(svycox.Bxp.gr4.2)

# lutein + zeaxanthin
svycox.lup.gr4.2<-svycoxph(Surv(HSAITMOR, age_die, m_inf)~factor(s_lup.gr4)+factor(sex)+factor(raceth)+factor(phase)
                           +factor(educ)+chol+bmi+factor(smk)+factor(diag_COPD), bpdsn45.COPD)
summary(svycox.lup.gr4.2)

# lycopene
svycox.lyp.gr4.2<-svycoxph(Surv(HSAITMOR, age_die, m_inf)~factor(s_lyp.gr4)+factor(sex)+factor(raceth)+factor(phase)
                           +factor(educ)+chol+bmi+factor(smk)+factor(diag_COPD), bpdsn45.COPD)
summary(svycox.lyp.gr4.2)

# total antioxidant capacity
svycox.TAC.gr4.2<-svycoxph(Surv(HSAITMOR, age_die, m_inf)~factor(s_TAC.gr4)+factor(sex)+factor(raceth)+factor(phase)
                           +factor(educ)+chol+bmi+factor(smk)+factor(diag_COPD), bpdsn45.COPD)
summary(svycox.TAC.gr4.2)

#####
## +cancer
# vitC
svycox.vitC.gr4.2<-svycoxph(Surv(HSAITMOR, age_die, m_inf)~factor(s_vitC.gr4)+factor(sex)+factor(raceth)+factor(phase)
                            +factor(educ)+chol+bmi+factor(smk)+factor(diag_ca), bpdsn45.ca)
summary(svycox.vitC.gr4.2)

# vitA
svycox.vitA.gr4.2<-svycoxph(Surv(HSAITMOR, age_die, m_inf)~factor(s_vitA.gr4)+factor(sex)+factor(raceth)+factor(phase)
                            +factor(educ)+chol+bmi+factor(smk)+factor(diag_ca), bpdsn45.ca)
summary(svycox.vitA.gr4.2)

# vitE
svycox.vitE.gr4.2<-svycoxph(Surv(HSAITMOR, age_die, m_inf)~factor(s_vitE.gr4)+factor(sex)+factor(raceth)+factor(phase)
                            +factor(educ)+chol+bmi+factor(smk)+factor(diag_ca), bpdsn45.ca)
summary(svycox.vitE.gr4.2)

# carotene
svycox.car.gr4.2<-svycoxph(Surv(HSAITMOR, age_die, m_inf)~factor(s_car.gr4)+factor(sex)+factor(raceth)+factor(phase)
                           +factor(educ)+chol+bmi+factor(smk)+factor(diag_ca), bpdsn45.ca)
summary(svycox.car.gr4.2)

# beta cryptoxanthin
svycox.Bxp.gr4.2<-svycoxph(Surv(HSAITMOR, age_die, m_inf)~factor(s_Bxp.gr4)+factor(sex)+factor(raceth)+factor(phase)
                           +factor(educ)+chol+bmi+factor(smk)+factor(diag_ca), bpdsn45.ca)
summary(svycox.Bxp.gr4.2)

# lutein + zeaxanthin
svycox.lup.gr4.2<-svycoxph(Surv(HSAITMOR, age_die, m_inf)~factor(s_lup.gr4)+factor(sex)+factor(raceth)+factor(phase)
                           +factor(educ)+chol+bmi+factor(smk)+factor(diag_ca), bpdsn45.ca)
summary(svycox.lup.gr4.2)

# lycopene
svycox.lyp.gr4.2<-svycoxph(Surv(HSAITMOR, age_die, m_inf)~factor(s_lyp.gr4)+factor(sex)+factor(raceth)+factor(phase)
                           +factor(educ)+chol+bmi+factor(smk)+factor(diag_ca), bpdsn45.ca)
summary(svycox.lyp.gr4.2)

# total antioxidant capacity
svycox.TAC.gr4.2<-svycoxph(Surv(HSAITMOR, age_die, m_inf)~factor(s_TAC.gr4)+factor(sex)+factor(raceth)+factor(phase)
                           +factor(educ)+chol+bmi+factor(smk)+factor(diag_ca), bpdsn45.ca)
summary(svycox.TAC.gr4.2)

############################
### Supplemental Table 5 ###
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
### Supplemental Table 6 ###
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


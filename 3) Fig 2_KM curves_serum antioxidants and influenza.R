
##################################################################
# Kaplan-Meier curves with IPW: Figure 2
##################################################################
library(survival)
library(epiDisplay)
library(Hmisc)
library(ggplot2)
library(tidyverse)
library(dplyr)
#### load the survey package
library(survey)
library(haven)
options(survey.lonely.psu="adjust")

setwd("C:\\Users\\great\\OneDrive\\01 on going\\2020- antioxidant\\antioxidant and influenza mortality\\Github")
load("serum.rda")
names(serum)
dim(serum)

attach(serum)

summary(serum[which(m_inf==1),]$PERMTH_EXM)

library("survminer")
library("grid")


## check survtime
summ(PERMTH_EXM)
summ(age_die)

#####
###for vitC

serum.vitC <- transform(serum,
                        s_vitC.gr1 = ifelse(s_vitC < 0.4, 0, 1))
attach(serum.vitC)
table(s_vitC.gr1)
vitC.wt1 <- glm(s_vitC.gr1 ~ age+factor(sex)+factor(raceth)
                +factor(educ)+chol+bmi+factor(phase)+factor(smk), 
                family = binomial(), data = serum.vitC)
summary(vitC.wt1)

p.vitC.obs <- ifelse(s_vitC.gr1 == 0, 1 - predict(vitC.wt1, type = "response"),
                     predict(vitC.wt1, type = "response"))

## pred<-predict(cd.wt, type = "response")
## length(pred)

serum.vitC$w <- 1/p.vitC.obs
summary(serum.vitC$w)
sd(serum.vitC$w)
summ(serum.vitC$w)
sum(serum.vitC$w)
tapply(serum.vitC$w, s_vitC.gr1, sum)

## stabilized w
pd.vitC<-predict(vitC.wt1, type = "response")

# estimation of numerator of ip weights
numer.fit <- glm(s_vitC.gr1 ~1, family = binomial(), data=serum.vitC)
summary(numer.fit)
pn.vitC <- predict(numer.fit, type = "response")

serum.vitC$sw1 <- ifelse(s_vitC.gr1 == 0, ((1-pn.vitC)/(1-pd.vitC)),
                        (pn.vitC/pd.vitC))
summary(serum.vitC$sw1)
tapply(serum.vitC$sw1, s_vitC.gr1, sum)

### weighted survfit
### using age as the time scale
summary(age_die/12)
tapply(age_die/12, m_inf, max)

fit1 <- survfit(Surv(HSAITMOR/12, age_die/12, m_inf) ~ s_vitC.gr1, weights=sw1, data=serum.vitC)
summary(fit1)
km.nh3.age.vitC<-ggsurvplot(fit1, data = serum.vitC, xlab="Age", xlim=c(45,105), ylim=c(0,0.3),
                   conf.int=F, linetype = "strata", 
                   legend.labs = c("<0.4", "¡Ã0.4"),
                   break.time.by=10, ggtheme = theme_bw(), censor=F, fun="cumhaz",
                   title="A. Serum vitamin C (mg/dL)", risk.table = T)
km.nh3.age.vitC

#####
###for car

serum.car <- transform(serum,
                        s_car.gr1 = ifelse(s_car < 23, 0, 1))
attach(serum.car)
table(s_car.gr1)
car.wt1 <- glm(s_car.gr1 ~ age+factor(sex)+factor(raceth)
                +factor(educ)+chol+bmi+factor(phase)+factor(smk), 
                family = binomial(), data = serum.car)
summary(car.wt1)

p.car.obs <- ifelse(s_car.gr1 == 0, 1 - predict(car.wt1, type = "response"),
                     predict(car.wt1, type = "response"))

## pred<-predict(cd.wt, type = "response")
## length(pred)

serum.car$w <- 1/p.car.obs
summary(serum.car$w)
sd(serum.car$w)
summ(serum.car$w)
sum(serum.car$w)
tapply(serum.car$w, s_car.gr1, sum)

## stabilized w
pd.car<-predict(car.wt1, type = "response")

# estimation of numerator of ip weights
numer.fit <- glm(s_car.gr1 ~1, family = binomial(), data=serum.car)
summary(numer.fit)
pn.car <- predict(numer.fit, type = "response")

serum.car$sw1 <- ifelse(s_car.gr1 == 0, ((1-pn.car)/(1-pd.car)),
                         (pn.car/pd.car))
summary(serum.car$sw1)
tapply(serum.car$sw1, s_car.gr1, sum)

### weighted survfit
### using age as the time scale
summary(age_die/12)
tapply(age_die/12, m_inf, max)

fit2 <- survfit(Surv(HSAITMOR/12, age_die/12, m_inf) ~ s_car.gr1, weights=sw1, data=serum.car)
summary(fit2)
km.nh3.age.car<-ggsurvplot(fit2, data = serum.car, xlab="Age", xlim=c(45,105), ylim=c(0,0.3),
                            conf.int=F, linetype = "strata",
                            legend.labs = c("<23", "¡Ã23"), 
                            break.time.by=10, ggtheme = theme_bw(), censor=F, fun="cumhaz",
                            title="B. Sum of serum ¥á- and ¥â-carotene (¥ìg/dL)", risk.table = TRUE)
km.nh3.age.car

#####
###for lyp

serum.lyp <- transform(serum,
                       s_lyp.gr1 = ifelse(s_lyp < 11, 0, 1))
attach(serum.lyp)
table(s_lyp.gr1)
lyp.wt1 <- glm(s_lyp.gr1 ~ age+factor(sex)+factor(raceth)
               +factor(educ)+chol+bmi+factor(phase)+factor(smk), 
               family = binomial(), data = serum.lyp)
summary(lyp.wt1)

p.lyp.obs <- ifelse(s_lyp.gr1 == 0, 1 - predict(lyp.wt1, type = "response"),
                    predict(lyp.wt1, type = "response"))

## pred<-predict(cd.wt, type = "response")
## length(pred)

serum.lyp$w <- 1/p.lyp.obs
summary(serum.lyp$w)
sd(serum.lyp$w)
summ(serum.lyp$w)
sum(serum.lyp$w)
tapply(serum.lyp$w, s_lyp.gr1, sum)

## stabilized w
pd.lyp<-predict(lyp.wt1, type = "response")

# estimation of numerator of ip weights
numer.fit <- glm(s_lyp.gr1 ~1, family = binomial(), data=serum.lyp)
summary(numer.fit)
pn.lyp <- predict(numer.fit, type = "response")

serum.lyp$sw1 <- ifelse(s_lyp.gr1 == 0, ((1-pn.lyp)/(1-pd.lyp)),
                        (pn.lyp/pd.lyp))
summary(serum.lyp$sw1)
tapply(serum.lyp$sw1, s_lyp.gr1, sum)

### weighted survfit
### using age as the time scale
summary(age_die/12)
tapply(age_die/12, m_inf, max)

fit3 <- survfit(Surv(HSAITMOR/12, age_die/12, m_inf) ~ s_lyp.gr1, weights=sw1, data=serum.lyp)
summary(fit3)
km.nh3.age.lyp<-ggsurvplot(fit3, data = serum.lyp, xlab="Age", xlim=c(45,105), ylim=c(0,0.3),
                           conf.int=F, linetype = "strata",
                           legend.labs = c("<11", "¡Ã11"), 
                           break.time.by=10, ggtheme = theme_bw(), censor=F, fun="cumhaz",
                           title="C. Serum lycopene (¥ìg/dL)", risk.table = TRUE)
km.nh3.age.lyp


#####
###for TAC

serum.TAC <- transform(serum,
                       s_TAC.gr1 = ifelse(s_TAC < 0.8, 0, 1))
attach(serum.TAC)
table(s_TAC.gr1)
TAC.wt1 <- glm(s_TAC.gr1 ~ age+factor(sex)+factor(raceth)
               +factor(educ)+chol+bmi+factor(phase)+factor(smk), 
               family = binomial(), data = serum.TAC)
summary(TAC.wt1)

p.TAC.obs <- ifelse(s_TAC.gr1 == 0, 1 - predict(TAC.wt1, type = "response"),
                    predict(TAC.wt1, type = "response"))

## pred<-predict(cd.wt, type = "response")
## length(pred)

serum.TAC$w <- 1/p.TAC.obs
summary(serum.TAC$w)
sd(serum.TAC$w)
summ(serum.TAC$w)
sum(serum.TAC$w)
tapply(serum.TAC$w, s_TAC.gr1, sum)

## stabilized w
pd.TAC<-predict(TAC.wt1, type = "response")

# estimation of numerator of ip weights
numer.fit <- glm(s_TAC.gr1 ~1, family = binomial(), data=serum.TAC)
summary(numer.fit)
pn.TAC <- predict(numer.fit, type = "response")

serum.TAC$sw1 <- ifelse(s_TAC.gr1 == 0, ((1-pn.TAC)/(1-pd.TAC)),
                        (pn.TAC/pd.TAC))
summary(serum.TAC$sw1)
tapply(serum.TAC$sw1, s_TAC.gr1, sum)

### weighted survfit
### using age as the time scale
summary(age_die/12)
tapply(age_die/12, m_inf, max)

fit4 <- survfit(Surv(HSAITMOR/12, age_die/12, m_inf) ~ s_TAC.gr1, weights=sw1, data=serum.TAC)
summary(fit4)
km.nh3.age.TAC<-ggsurvplot(fit4, data = serum.TAC, xlab="Age", xlim=c(45,105), ylim=c(0,0.3),
                           conf.int=F, linetype = "strata",
                           legend.labs = c("<0.8", "¡Ã0.8"),  
                           break.time.by=10, ggtheme = theme_bw(), censor=F, fun="cumhaz",
                           title="D. Serum TAC (mg VCE/dL)", risk.table = TRUE)
km.nh3.age.TAC




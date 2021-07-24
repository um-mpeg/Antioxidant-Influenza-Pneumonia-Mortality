library(dplyr)
library(corrplot)

setwd("C:\\Users\\great\\OneDrive\\01 on going\\2020- antioxidant\\antioxidant and influenza mortality\\Github")
load("serum.rda")
load("diet.rda")

a <- merge(serum,diet, by="SEQN")
a1 <- select(a, s_vitC, s_vitA, s_vitE, s_car, s_Bxp, s_lup, s_lyp, s_TAC, d_vitC, d_reti, d_vitE, d_Bcar)
dim(a1)

cor1 <- rcorr(as.matrix(a1))

corrplot(cor1$r,  order="original", addCoef.col = "black", tl.col="white", sig.level=.05, insig="blank")

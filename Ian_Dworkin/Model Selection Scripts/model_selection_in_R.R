# Multiple regression examples


require(car)
require(arm)

 # We are trying to account for variation in # SCT given the length of different segments of the leg, femur, tibia and tarsus. The problem is they are all correlated.
 
setwd("/Users/ian/R/R scripts/Dll data/") 
dll.data = read.csv("dll.txt", header=TRUE)   #data frame input
dll.data$temp <- factor(dll.data$temp)
dll.data$replicate <- factor(dll.data$replicate)
attach(dll.data)



# Full model
lm.full <- lm(SCT ~ femur +tibia + tarsus, data=dll.data)

vcov(lm.full)
vif(lm.full)
par(mfrow=c(3,2))
confidence.ellipse(lm.full, which.coef=c(1,2))
confidence.ellipse(lm.full, which.coef=c(1,3))
confidence.ellipse(lm.full, which.coef=c(1,4))
confidence.ellipse(lm.full, which.coef=c(2,3))
confidence.ellipse(lm.full, which.coef=c(2,4))
confidence.ellipse(lm.full, which.coef=c(3,4))

summary(lm.full)
coefplot(lm.full, vertical=F)
anova(lm.full)
AIC(lm.full)

lm.no.tar <- lm(SCT ~ femur + tibia , data=dll.data)
lm.no.tib <- lm(SCT ~ femur + tarsus, data=dll.data)
lm.no.fem <- lm(SCT ~ tibia + tarsus, data=dll.data)
lm.tib    <- lm(SCT ~ tibia , data=dll.data)
lm.tar    <- lm(SCT ~ tarsus , data=dll.data)
lm.fem    <- lm(SCT ~ femur , data=dll.data)
lm.intercept    <- lm(SCT ~ 1 , data=dll.data)

aic.lm <- AIC(lm.full,lm.no.tar, lm.no.tib, lm.no.fem, lm.tib, lm.tar, lm.fem, lm.intercept)
aic.lm

# A little function to compute BIC#n is number of observations
# k is number of parameters
BIC <- function(lm.object,n,k) {x <- logLik(lm.object)[1] # extract logLik
	  bic <- -2*x + k*log(n)
	  print(bic)  }


# Centre the data to look at 
dll.data$femur.cent <- dll.data$femur -mean(na.omit(dll.data$femur))
dll.data$tibia.cent <- dll.data$tibia -mean(na.omit(dll.data$tibia))
dll.data$tarsus.cent <- dll.data$tarsus -mean(na.omit(dll.data$tarsus))

lm.full.cent <- lm(SCT ~ femur.cent + tibia.cent + tarsus.cent, data=dll.data)
summary(lm.full.cent)
coefplot(lm.full.cent, vertical=F)
vcov(lm.full.cent)
vif(lm.full.cent) # Does not change these VIF's
par(mfrow=c(3,2))
confidence.ellipse(lm.full.cent, which.coef=c(1,2))
confidence.ellipse(lm.full.cent, which.coef=c(1,3))
confidence.ellipse(lm.full.cent, which.coef=c(1,4))
confidence.ellipse(lm.full.cent, which.coef=c(2,3))
confidence.ellipse(lm.full.cent, which.coef=c(2,4))
confidence.ellipse(lm.full.cent, which.coef=c(3,4))


# Could have also used the scale() with scale=F option

lm.full.std <- lm(SCT ~ scale(femur) + scale(tibia) + scale(tarsus) , data=dll.data)
summary(lm.full.std)
vcov(lm.full.std)
vif(lm.full.std) # Does not change these VIF's
par(mfrow=c(1,1))
coefplot(lm.full.std, vertical=F)
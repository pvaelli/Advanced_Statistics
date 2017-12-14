#model selection - November 30th 2009

# This just uses our SCT data set to demonstrate A) The nice little AICtab & BICtab function
# B) Demonstrate what a different picture AIC and BIC can give you for your models!!!
# It is also worth noting that we treating "line" as a fixed effect, when it should really be considered a random effect, and we should nbe using a mixed model approach. We will talk about this next week.

 # Ben Bolker has written a nice function to compute AIC, BIC weights etc...
#AICtab
#BICtab #both in bbmle
require(bbmle)

 
setwd("/Users/ian/R/R scripts/Dll data/") 
dll.data = read.csv("dll.txt", header=TRUE)   #data frame input
dll.data$temp <- factor(dll.data$temp)
dll.data$replicate <- factor(dll.data$replicate)
dll.data <- na.omit(dll.data) # remember for model comparison you want to make sure you are using the same data set each time, so take care with missing observations.

# Below is just an example of the model set that I was considering. 
lm.null <- lm(SCT~ 1, data=dll.data)
lm.model.1 <- lm(SCT ~ genotype + temp + line, data=dll.data)
lm.model.1.t <- lm(SCT ~ genotype + temp + line + tarsus, data=dll.data)
lm.model.2.t <- lm(SCT ~ genotype*temp + line + tarsus, data=dll.data)
lm.model.2 <- lm(SCT ~ genotype*temp + line, data=dll.data)
lm.model.3.t <- lm(SCT ~ genotype*temp + line + line*temp + line*genotype+ tarsus, data=dll.data) # This was the model I was testing with respect to the design of the experiment
lm.model.3 <- lm(SCT ~ genotype*temp + line + line*temp + line*genotype, data=dll.data)
lm.model.4 <- lm(SCT ~ genotype*temp*line , data=dll.data) #3way interaction
lm.model.4.t <- lm(SCT ~ genotype*temp*line + tarsus , data=dll.data) #3way interaction
lm.model.full<- lm(SCT ~ genotype*temp*line*tarsus , data=dll.data) #fully specified model that has co-efficients far too complex to ever interpret, but is here for comparison

#AIC
aic.tab.lm <- AICtab(lm.null, lm.model.1, lm.model.1.t, lm.model.2.t, lm.model.2, lm.model.3.t, lm.model.3, lm.model.4, lm.model.4.t, lm.model.full, base=T, weights=T, nobs=nrow(dll.data))
aic.tab.lm
# The AIC suggests that model lm.model.4.t is the best approximating model. This is the second most complex model with 103 parameters!!
summary(lm.model.4.t)

### Note for AIC corrected  just use AICctab (or BICctab)

bic.tab.lm <- BICtab(lm.null, lm.model.1, lm.model.1.t, lm.model.2.t, lm.model.2, lm.model.3.t, lm.model.3, lm.model.4, lm.model.4.t, lm.model.full, base=T, weights=T, nobs=nrow(dll.data))
bic.tab.lm
# This suggests that a far simpler model is the best approximating model. This is not suprising given the very large sample size of the experiment. The next best model is actually pretty far away.
summary(lm.model.2.t)

# It is probably also worth comparing residuals between models
# sum of squares
sum(resid(lm.model.4.t)^2)
sum(resid(lm.model.2.t)^2)
#sd
sd(resid(lm.model.4.t))
sd(resid(lm.model.2.t))
par(mfrow=c(2,2))
plot(lm.model.2.t)

par(mfrow=c(2,2))
plot(lm.model.4.t)

#Neither of these was really the model I was interested in! Surprise!! We will see how treating line as a random effect can change things (since we will not be estimating a parameter for each level of line in the same way)
# Last Modified on December 3rd 2009 - Ian Dworkin




require(bbmle)
require(MASS)
require(car)
setwd("/Users/ian/R/R scripts/Dll data/") 
dll.data = read.csv("dll.txt", header=TRUE)   #data frame input
dll.data$temp <- factor(dll.data$temp)
dll.data$replicate <- factor(dll.data$replicate)
dll.data <- na.omit(dll.data)

#If you remember from our earlier discussions of probability distributions, it appeared as if the normal distribution was a better fit for the # of Sex Comb teeth over either poisson or negative binomial, despite the fact that it represents count data constrained 

fit.normal <- fitdistr(na.omit(dll.data$SCT),densfun="normal")

fit.poisson <- fitdistr(na.omit(dll.data$SCT),densfun="poisson")


fit.normal; fit.poisson
AIC(fit.normal); AIC(fit.poisson)
BIC(fit.normal);BIC(fit.poisson)
# but does this mean we should be using the normal distribution to model the residual variation?


# Now we are in a position to look at interactions
anova.interaction <- lm(SCT ~ temp*genotype, data=dll.data)




# Let's try a poisson family GLM
anova.poisson.GLiM <- glm(SCT ~ temp*genotype,family=poisson, data=dll.data)
summary(anova.poisson.GLiM)




AICtab(anova.interaction, anova.poisson.GLiM)
 # let's take a look with BIC

AICtab(anova.interaction, anova.poisson.GLiM, nobs=1918)
summary(anova.poisson.GLiM) # compare residual deviance with df... not in the same ball park. 

# Both AIC and BIC suggest that a model with normal residual variation is best.


# of course if there is over-dispersion we can try either the quasi-poisson or a negative binomial GLiM
anova.poisson.GLiM.OD <- glm(SCT ~ temp*genotype,family=quasipoisson, data=dll.data)
summary(anova.poisson.GLiM.OD) # less than one. So we will not get much help here...


# For a negative binomial
anova.nb.GLiM <- glm.nb(SCT ~ temp*genotype, data=dll.data) # This function is in the MASS library
summary(anova.nb.GLiM)


# Let's try a poisson family GLM
anova.poisson.GLiM.FULL <- glm(SCT ~ temp*genotype*line + tarsus,family=poisson, data=dll.data)

# For a negative binomial
anova.nb.GLiM.FULL <- glm.nb(SCT ~ temp*genotype*line + tarsus, data=dll.data) # This function is in the MASS library

#Gaussian model
anova.interaction.FULL <- lm(SCT ~ temp*genotype*line + tarsus, data=dll.data)


AICtab(anova.poisson.GLiM.FULL, anova.interaction.FULL, anova.nb.GLiM.FULL, weights=T)
BICtab(anova.poisson.GLiM.FULL, anova.interaction.FULL, anova.nb.GLiM.FULL, nobs=1931)

# just no evidence for this data set that a poisson based model makes sense....
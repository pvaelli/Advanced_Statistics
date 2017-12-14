# Mixed model analysis in R

# Using the Data set from Dworkin, I. 2005. Evolution and Development.

# Trying to account for variation in  the number of sex comb teeth on male flies as a function of  Genotype (Dll or wt), Rearing temperature, and line, as well as a replicate (block) effect

#Fixed effects for this analysis will be genotype and temperature.

# Random effects will be line, any interactions with line and the replicate effect.

# The design of this experiment is (sadly) unbalanced, as the numbers for each combination of treatment factors is unequal.

# There ar two different libraries for analyzing mixed models in R

# The original one nlme  has the function lme() for lmm's

# The newer one lme4 (which still has some bugs) uses lmer for lmm's



# libraries required (lme4 and sciplot will need to be installed first if you do not have them on your computer)
#require(nlme)

require(lme4)
#require(arm)
#require(sciplot)
#require(MCMCglmm)
# require(nlme)

# Inputting the data, and making sure variables are encoded correctly.
setwd("/Users/ian/R/R scripts/Dll data/") 
dll.data = read.csv("dll.txt", header=TRUE)   #data frame input
dll.data$temp <- as.factor(dll.data$temp)
dll.data$replicate <- factor(dll.data$replicate)

dll.data<- na.omit(dll.data)
dll.data.25.wt <- subset(dll.data, temp=="25", genotype="wt") # Just the data at 25 degrees C

# Let's start by examining the data
bwplot(SCT ~ temp:genotype,data=dll.data)
bwplot(SCT ~ line|temp:genotype,data=dll.data)
xyplot(SCT ~ tarsus|line, data=dll.data.25)





# We can try a simple model where we just account for line

lm.line.1 <- lm(SCT ~ line, data=dll.data)
summary(lme.line.1) # Let's think about what model we are fitting here.

lm.line.2 <- lm(SCT ~ 0 + line, data=dll.data)
summary(lm.line.2)


######

#########Using lmer() in lme4 for a simple mixed model
lmer.1 <- lmer(SCT ~ 1 + (1|line), REML=T, na.action=na.omit, data=dll.data)
print(lmer.1) # variance components
ranef(lmer.1)
coef(lmer.1)
coefplot(lmer.1)

##### Code for lme() in nlme - which we will not use for the moment
lmm.1 <- lme(SCT ~ 1, data=dll.data,random= ~1|line,  method="REML", na.action=na.omit)
summary(lmm.1)
AIC(lmm.1)
BIC(lmm.1)
qqnorm(lmm.1)
ranef(lmm.1)

logLik(lmm.1)
plot(lmm.1)
hist(resid(lmm.1))
sd(resid(lmm.1))

 # You can see we do not get individual estimates for line, but instead a variance component
 






# Let us first start where we left off with this data set. using a standard linear model, accounting for variation in # of SCT with the following model. Check out the set of models from the original model set if you wish to compare AIC's or BIC's for those models. 

# However since this was actually a designed experiment, we should really be focusing in on the effects that it was designed to test. Namely the line*genotype + line*temp.
# However we want to control for size related effects, so we use tarsus as a covariate. Whether we should tarsus as part of a random slopes model is up for discussion.


SCT.lm.1 <- lm(SCT ~ temp*genotype*line - line:genotype:temp + tarsus, data=dll.data)
summary(SCT.lm.1)
anova(SCT.lm.1)


hist(resid(SCT.lm.1))
sd(resid(SCT.lm.1)) # This should be similar to the Residual standard error
# You can also divide the  RSS of the model by the Residual df (which provides a Residual variance)
# Then take the sqrt of this residual variance and it should equal the Residual standard error of the model.


plot(SCT.lm.1)
lag.plot(resid(SCT.lm.1), lags=2)

acf(resid(SCT.lm.1))

# perhaps a better way of demonstrating this (need to use vetted data)
SCT.lm.1.omit <- lm(SCT ~ temp*genotype, data=dll.data)
bwplot(line ~ resid(SCT.lm.1),xlim=c(-8, 8), data=dll.data)

# Clearly some evidence of unmodeled variation






# That is we are estimating a whole set of parameters (the line means) we are probably not interested in. We just want to estimate the distribution of the effect of line, and factors that interact with line.


#  Thus we will treat variable such as line, line:genotype, and line:temp as random effects best described by a distribution.



# Sometimes we may be interested in fitting a seperate model for each level of a grouping factor (line in this case)

# we ca use lmList function for this

models.list <- lmList(SCT ~ genotype | line, data=dll.data)
ranef(models.list) # Random effects of line
pairs(models.list) # shows correlations between estimates

plot(intervals(models.list))
# Clearly lots of subject to subject variability




###### the lme4 library for mixed models

 lmm.4 <- lmer(SCT ~ genotype*temp + (1|line) + (1|line:genotype) + (1|line:temp) + tarsus, REML=T, data=dll.data)
 lmm.4
 summary(lmm.4)
ranef(lmm.4)


# If we wanted to a random slopes model

lmm.5.RS <- lmer(SCT ~ genotype*temp + (1|line) + (1|line:genotype) + (1|line:temp) + tarsus + (tarsus|line), REML=T, data=dll.data)
 lmm.5.RS
 
 #Based on BIC, this is a worse fit

#Same model but with no random intercept built into tarsus (which makes sense since we already have it for 1|line) 

lmm.6.RS <- lmer(SCT ~ genotype*temp + (1|line) + (1|line:genotype) + (1|line:temp) + tarsus + (0+tarsus|line), REML=T, data=dll.data)
 lmm.6.RS
 
 
#####Some other little tricks with random slopes
# THis models both random slopes and random intercepts ( we will use centered tarsus)
# we are also only using a small amount of the data

dll.data.25.wt$tarsus.centered <- dll.data.25.wt$tarsus - mean(dll.data.25.wt$tarsus)


 lmm.RS<- lmer(SCT ~ 1 + tarsus.centered + (1+ tarsus.centered|line), data=dll.data.25.wt)

 
 # Notice that there is in fact a lot of variance for slopes of tarsus
 # Also some correlation between the random effects estimates of random slopes and random intercepts
 
 # Now we fit a model where we allow slopes to vary, but not the intercept ( more to show you how technically, than for any other reason)
 lmm.RS.2 <- lmer(SCT ~ 1 + tarsus.centered + (0 + tarsus.centered|line), data=dll.data.25.wt)
 lmm.RS.2
 
 # The final model in this series that I think is useful to think about (or at least to know about) is where we fit both random intercept and random slope, BUT we force them to be uncorrelated with one another.
 
 lmm.RS.3 <- lmer(SCT ~ 1 + tarsus.centered + + (1|line)+ (0 + tarsus.centered|line), data=dll.data.25.wt)
 lmm.RS.3
 
 # THis is very similar to lmm.RS, but you will see that the variance components are uncorrelated. This can potentially affect the BLUP's as well. You really need to ask yourself though what your rational would be for forcing these to be uncorrelated.
 
 a <- ranef(lmm.RS)
 b <- ranef(lmm.RS.3)
 


# Nesting vs. blocking in lmer

#Let us suppose that our two replicates were measured at very different times, but everyone in replicate 1 was measured together, and everyone in replicate 2 was measured together. 
# In this case, it would be best to treat rep as a blocking effect, especially if all lines were reared in each replicate (complete block)
 
 lmm.7 <- lmer(SCT ~ genotype*temp + (1|line) + (1|line:genotype) + (1|line:temp) + tarsus + (1|replicate), REML=T, data=dll.data)
 lmm.7
# It turns out that replicate effects are so small that they are effectively zero if modeled in this fashion.


# However in this case (if I remember the experiment) each temperature treatment was reared in its own incubator, so arguably the above design is incorrect, and instead we should think of replicate nested within temperature

lmm.8 <- lmer(SCT ~ genotype*temp + (1|line) + (1|line:genotype) + (1|line:temp) + tarsus + (1|replicate:temp), REML=T, data=dll.data)
 lmm.8
# No real effect so we can probably drop replicate entirely.


# If we had nested replicate bottles within lines then we would model it slightly differently. We would want the replicate:line effect, or possibily the replicate:line:temp effect

lmm.9 <- lmer(SCT ~ genotype*temp + (1|line) + (1|line:genotype) + (1|line:temp) + tarsus +(1|replicate:line:temp), REML=T, data=dll.data)
 lmm.9
 # While it is non-zero, it is still pretty small and does not influence the model too much. THE BIC & AIC excluding replicate also suggest a better approximating model without rep, so we are probably pretty safe excluding it from the model.
 
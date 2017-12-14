#### 

library(lme4) 

library(afex)
library(pbkrtest)
library(MuMIn)
library(arm)

library(effects)

#### Block 1: Read in and set up data, 
#### Block 2:  Some simple models using both lmer (lme4 ) 
#### Block 3: Different parameterizations for random effect interactions 
#### Block 4:  "Random regression" 
#### Block 5: Big model

#Block 1
#setwd("/Users/ian/R/R scripts/Dll data/") 
#dll.data = read.csv("dll.txt", header=TRUE)   #data frame input
dll.data = read.csv("http://datadryad.org/bitstream/handle/10255/dryad.8377/dll.csv", header=TRUE)

dll.data$temp <- as.factor(dll.data$temp)
dll.data$replicate <- factor(dll.data$replicate)
dll.data<- na.omit(dll.data)
dll.data$tarsus.centered <- scale(dll.data$tarsus, center=T, scale=F) #(mean centers, but does not  standardize)
dll.data$genotype <- relevel(dll.data$genotype, "wt")

# First a relatively simple model

#####Block 1

model.1.ML.us <- lmer(SCT ~ 1 + genotype + (genotype| line), data=dll.data) # NOTE: REML=TRUE by default

model.2.ML.us <- lmer(SCT ~ 1 + genotype + (1+ genotype| line), data=dll.data) # equivalent to the previous model!

print(model.1.ML.us)
print(model.2.ML.us)

summary(model.1.ML.us)

model1_ML_profile <- profile(model.1.ML.us) # profile the likelihood surface
confint(model1_ML_profile) # For more complicated models this often does not work, so you need to use parametric or non-parametric boostrap in particular for the variances.

# Keep in mind that "counting" parameters for random effects is difficult and may be done differently
-2*logLik(model.1.ML.us)

AIC(model.1.ML.us) # deviance + 2*6 parameters

extractAIC(model.1.ML.us) # Conditional AIC from arm library. Also see tools in MuMIn

extractDIC(model.1.ML.us)
DIC(model.1.ML.us) # From the MuMIn library

r.squaredGLMM(model.1.ML.us) # marginal (fixed) and condition (fixed + random R2 measures for model)
### NOTE THIS r.squaredGLMM() IS A VERY NEW METHOD AND FUNCTION!! NOT CLEAR HOW ROBUST IT IS.


# The effects package allows for some useful visualization of fixed (but not random effects)
plot(effect("genotype", model.1.ML.us), grid=TRUE)
plot(allEffects(model.1.ML.us), grid=TRUE)

# random effects
ranef(model.1.ML.us) # Best Linear Unbiased Predictors (BLUPs)

# Usually helps to plot these out (and you can see the effects of the large -ve correlation for the random effects)
print(dotplot(ranef(model.1.ML.us, condVar=TRUE),
              scales = list(x = list(relation = 'free')))[["line"]])


## Some notes on Model "testing and comparison"
# there are some functions in car (Anova) and in afex (and probably some other libraries) to do anova tables. 
#i.e.
car::Anova(model.1.ML.us)


#Generally they only assess fit of fixed effect terms (conditioned on the random terms), but do not assess differences in models that differ for random effects.  

# However,are some conditional *ICs like extractAIC() (from the arm library) and I believe that MuMIn also has some. However keep in mind that counting the number of parameters for random effects can be really difficult if the data is not balanced. Therefore it is best to use "test statistics" which will not depend on knowing the number of random effect parameters being estimated (as a df). So commonly you will see (and this is what pbkrtest does) a likelihood ratio test, but evaluated using parametric bootstrap, permutation or non-parametric bootstrap (the latter two take some seriously thinking about the resampling procedure).

# I have never used the pbkrtest library, and generally write my own parametric bootstrap simulator either like we did in class, or using simulate(yourMerModel). simulate() is a generic function that calls class specific methods, in this case (for lmer) simulate.merMod


# First we can at least look at the likelihood ratio

# Fit model WITHOUT a random effect for genotype
model.0.ML.us <- lmer(SCT ~ 1 + genotype + (1| line), data=dll.data) # NOTE: REML=TRUE by default
extractAIC(model.0.ML.us)
extractAIC(model.1.ML.us)

REMLcrit(model.0.ML.us) # this is just the REML deviance. It used to work with just deviance()!
REMLcrit(model.1.ML.us)

# Likelihood Ratio
LR.model <-  -as.numeric(REMLcrit(model.1.ML.us)- REMLcrit(model.0.ML.us))
LR.model
# The question now becomes how many df these two models differ by. We are estimating 2 less terms for the random effects, the line variance for genotype and the covariance between genotype and the intercept.So one approach would be to compare this to a chisq distribution with 2 df 
pchisq(q = LR.model, df=2, lower=F) 

# Which is very similar to as the default anova method (although this refits under ML not REML)
anova(model.1.ML.us, model.0.ML.us)

# (reminder about all of the issues with LRT near boundary conditions, and whether this approach is at all appropriate.

# We can also say that we are now estimating 27 less BLUPs (the one for genotype), so maybe we should use 27df
pchisq(q = LR.model, df=27, lower=F) 


# In this case it may not matter, but you can see by the massive change in the magnitude of the p-value that it could very easily.  Best to rely on the parametric bootstrap (or an actually resampling approach)

LikRatioSim <- function(mod=model.0.ML.us) {
	y.sim <- simulate(mod, nsim = 1)  # one set of simulated data under simple model
	model.lower <- lmer(y.sim$sim_1 ~ 1 + genotype + (1|line), dll.data) # fit simulated data under simple model
	model.full  <- lmer(y.sim$sim_1 ~ 1 + genotype + (1 + genotype|line), dll.data) # fit simulated data w complex model 
	LRSim <-  -as.numeric(REMLcrit(model.full) - REMLcrit(model.lower)) # Likelihood Ratio
	return(LRSim)
	rm(y.sim, model.lower, model.full, LRSim)
}

n.sim = 100 # 100 simulations (possibly not enough simulations, but for time...)
LikRatParBoot <- replicate(n = n.sim, LikRatioSim())


# And can compare this to the observed LR
(length(LikRatParBoot[LikRatParBoot > LR.model])+1)/n.sim  # Our p-value
# Which in this case is consistent with our previous assessment

# Do keep in mind that the parametric bootstrap does make very strong assumptions, so do evaluate whether your data conforms to them.


# I have less experience with it, but you can also use the PBmodcomp function in pbkrtest (also see library RLRsim)

pb.b <- PBmodcomp(model.1.ML.us, model.0.ML.us, nsim=100) # from pbkr - will take a long time even for 100
summary(pb.b)

###  note: We can have a "no intercept version for random effects which may help in the interpretation. This is actually much easier using MCMCglmm. This is identical though to the model above.

model.1.ML.noRandomIntercept <- lmer(SCT ~ 1 + genotype + (0 + genotype| line) , data=dll.data) 
summary(model.1.ML.noRandomIntercept) 
extractAIC(model.1.ML.noRandomIntercept)
r.squaredGLMM(model.1.ML.noRandomIntercept)




#### Different ways of parameterizing the model for random effect "interactions"
###### Block 3
# This is a fairly different parameterization where fundamentally we are asking about additive variances (i.e. the variance for line + variance for the interaction). As we will go over on the board, the previous model has a very different variance structure. Generally speaking this will NOT be the way you would want to model it unless you have very particular a priori reasons to model the variance this way.


model1_alt_parameterization <- lmer(SCT ~ 1 + genotype + (1|line) + (1|line:genotype), 
    data=dll.data)
# Here we are assuming that :
# (1|line) ~ N(0, var(line)) and 
# (1|line:genotype) ~ N(0, var(line:genotype))

summary(model1_alt_parameterization)
extractAIC(model1_alt_parameterization)
extractDIC(model1_alt_parameterization)
r.squaredGLMM(model1_alt_parameterization)




######### Block 4: Random Regression
# Random regression example


# Let's look at the data
# First we consider just variation from line to line for SCT ~ tarsus, but pooling across genotype & temp
xyplot(SCT ~ tarsus.centered | line, data = na.omit(dll.data), 
  xlab = " tarsus length (centered)", ylab = "# SCT", lwd=3, pch=21, fill="grey",
  type = c("r", "p"),  # Here we have asked for points (p), and a regression line ("r")
  auto.key = list(points = TRUE, lines = TRUE, columns = 2)) # generating the legend

# with more detail
xyplot(SCT ~ tarsus.centered | line, data = na.omit(dll.data), 
  xlab = " tarsus length (centered)", ylab = "# SCT", lwd=3, pch=21, fill="grey",
  type = c("r", "p"),  # Here we have asked for points (p), and a regression line ("r")
  groups=genotype:temp,  # adds colours for particular levels of the grouping variables
  auto.key = list(points = TRUE, lines = TRUE, columns = 2)) # generating the legend


# This gives an alternative view of the data
xyplot(SCT ~ tarsus.centered | genotype + temp, data = na.omit(dll.data), 
  xlab = " tarsus length (centered)", ylab = "# SCT", lwd=3, pch=21, fill="grey",
  type = c("r", "p"),  
  groups=line)


# In this case, we do not really have enough data to estimate each of the regression lines, but it gives us some ideas


#### So we will have a continuous predictor (tarsus length), and ask how the slope and intercept vary across lines.


# model for lmer() in the lme4 library:
model.tarsus.ML.us <- lmer(SCT ~ 1 + tarsus.centered  + (1 + tarsus.centered | line), 
    data=dll.data) 
summary(model.tarsus.ML.us)

plot(effect("tarsus.centered", model.tarsus.ML.us), grid=TRUE)
plot(allEffects(model.tarsus.ML.us), grid=TRUE)

# random effects
ranef(model.tarsus.ML.us) # blups

# Usually helps to plot these out (and you can see the effects of the large -ve correlation for the random effects)
print(dotplot(ranef(model.tarsus.ML.us, condVar=TRUE),
              scales = list(x = list(relation = 'free')))[["line"]])

# or we can force co-variances between the  slopes and  intercepts for the random effect (line) to be zero

#equivalent model for lmer() in the lme4 library:
model.tarsus.ML.idh <- lmer(SCT ~ 1 + tarsus.centered  + (1| line) + (0 + tarsus.centered | line),  
    data=dll.data) 
print(model.tarsus.ML.idh)

extractAIC(model.tarsus.ML.us)
extractAIC(model.tarsus.ML.idh)

# We could evaluate what happens if we do not allow a random effect for the contribution for tarsus by line
model.tarsus.ML.NoRandomSlope <- lmer(SCT ~ 1 + tarsus.centered  + (1| line), 
    data=dll.data) 
extractAIC(model.tarsus.ML.NoRandomSlope)
extractAIC(model.tarsus.ML.us) 

# While I have not done the Parametric Bootstrap, there is at least pretty reasonable evidence to retain it.
   
############## Block 5 - the big model.

# Maybe we want to fit a complex model. This is pretty close to what we would be interested in ultimately (although we have not accounted for heterogeneity of residual variances here)

model.tarsus.ML.BIG <- lmer(SCT ~ 1 + tarsus.centered + genotype*temp
  + (1 + tarsus.centered + genotype*temp | line), data=dll.data)
  
  
extractAIC(model.tarsus.ML.BIG)   
print(model.tarsus.ML.BIG)
# keep in mind we have estimated 11 variances & co-variances from our one random effect with 27 levels.
# Not likely that this is really reasonable to expect an overall good fit for the random effects.

r.squaredGLMM(model.tarsus.ML.BIG)

car::Anova(model.tarsus.ML.BIG)
plot(allEffects(model.tarsus.ML.BIG), grid=F)

ranef(model.tarsus.ML.BIG)

# Usually helps to plot these out
print(dotplot(ranef(model.tarsus.ML.BIG, condVar=TRUE),
              scales = list(x = list(relation = 'free')))[["line"]])

## So what might be a more realistic model to fit (in terms of random effect df?)

# Remove one random effect term
model.tarsus.ML.BIG.alt <- lmer(SCT ~ 1 + tarsus.centered + genotype*temp
  + (1 + tarsus.centered + genotype + temp | line), data=dll.data)
print(model.tarsus.ML.BIG.alt)
  
extractAIC(model.tarsus.ML.BIG)  
extractAIC(model.tarsus.ML.BIG.alt)  

# Suggests that there likely is some really evidence for retaining those (co)variances. Worth testing with a parametric bootstrap

anova(model.tarsus.ML.BIG, model.tarsus.ML.BIG.alt) # keep in mind this is a straight up LRT fit with ML

# Let's look at a parametric bootstrap (from pbkrtest)
pb.b <- PBmodcomp(model.tarsus.ML.BIG.alt, model.tarsus.ML.BIG.alt, nsim=5) # from pbkr - will take a VERY long time even for 100
summary(pb.b)

#                              
# example using glmmML library and glmmboot?                                 
#Example using MASS::glmmPQL?
# Example using glmmadmb?                                   
# Ian Dworkin November 2014


#### Longitudinal data analysis (repeated measures)
# Example originally from Douglas Bates 2009, with many modifications!

library(lattice)
library(Matrix)
library(lme4)
library(arm)
library(effects)
library(MuMIn)

# #this is data looking at reaction times of  a set of individuals (Subjects) who are subjected to sleep deprivation. Their reaction times are examined after each data of sleep deprivation. 
# Thus we are repeatedly sampling the same individuals over a time period.
# The sleepstudy data set is in lme4 
print(xyplot(Reaction ~ Days | Subject, sleepstudy, aspect = "xy",
                    layout = c(9,2), type = c("g", "p", "r"),
                    index.cond = function(x,y) coef(lm(y ~ x))[1],
                    xlab = "Days of sleep deprivation",
                    ylab = "Average reaction time (ms)"))


# While we are interested in the relationship between Reaction ~ Days, clearly individuals will vary in their response, so how can we deal with this?

# We could "ignore" the information at the individual level, and pool the data as follows.
plot(Reaction ~ Days, data=sleepstudy, main="Reaction time as a function of sleep deprivation")

lm.pooled <- lm(Reaction ~ Days, data=sleepstudy)

summary(lm.pooled)

abline(lm.pooled, lwd=3)
# CI's via Monte Carlo Simulation

simmie <-function(){ # function to generate one simulation, and then plot the fitted line
    sim.lm <- sim(lm.pooled, n.sims=1)
    abline(a=sim.lm@coef[,1], b=sim.lm@coef[,2], col="grey") }

replicate(95, simmie())  # replicate the simulation function 95 times

# While it seems like we are fitting the data, clearly there is a lot of variation that have not accounting for (and we know it since we have ignored it).


# We can then do an analysis where we fit a seperate regression for each subject in the study.

# We could do this with a for loop or an apply function, but the easiest way is to us the lmList() in nlme & lme4

print(plot(confint(lmList(Reaction ~ Days | Subject, sleepstudy),
                   pooled = F), order = 1))
# THis just makes it clear how much variation there really is.

# What we need is an intermediate between pooled and no-pooling, where we can account for the stochastic effects of between subject variability, but still assess the overall Reaction ~ Days relationship.                


# Using lmer it is quite easy to fit a repeated measure longitudinal data set
#  We want to fit random effects for each individual for both there "intercept" (reaction time prior to sleep deprivation), 
#  As well as for the slope (rate at which their reaction time gets worse)
fm1 <- lmer(Reaction ~ 1 + Days + (1+ Days|Subject), sleepstudy)

# Let's discuss in class how we might interpret the results from this model
print(fm1)

# For the fixed effects, we can use the effects library to generate some plots (accounting for random terms)
# We can use the allEffects (we only have one fixed effect anyway) from the effects library
effects_fm1 <- allEffects(fm1)
plot(effects_fm1, grid=F)

# This is using a method specific to the effects library for plotting (plot.efflist() I think)
# If you do not like the look of the plot (although you can change aspects of it), you can extract estimates and CIs for specific fixed effects as follows (and make your own plot)
as.data.frame(effects_fm1$Days) 


# By accounting for inter-subject variation (and variation within  subject across Days) we have less uncertaintly with the slope of Days, and plus we get the variation due to subject, that is seperated from residual.


# From the model summary, there is clearly a considerable amount of variation among subjects, both in their initial response time, and how response time changes with sleep deprivation (with at least one subject getting a bit better). 

# Indeed think about the variances (or standard deviation) of the random effects relative to the magnitudes of the fixed effect (including the intercept which is interesting when we consider subject-to-subject variation)

#We can look at the random effects explicitly (the best linear unbiased predictors, or "blups")
ranef(fm1)

# Usually helps to plot these out
print(dotplot(ranef(fm1, condVar=TRUE),
              scales = list(x = list(relation = 'free')))[["Subject"]])


# How do we make inferences? With mixed models it can be tricky, but bootstrapping (parametric or non-parametric) can help.
# Many of the difficulties in making inferences for mixed models has to do with figuring our the number of "estimated parameters in the model" (degrees of freedom), which is a result of a mixed model lying somewhere between a no-pooling and complete pooling approach. We will return to this issue a fair bit below.


# For the fixed effects it is often a bit easier (and you can figure out to some degree with the estimates and standard errors). The Anova function in car can help (but be forewarned about the assumptions of the wald or F-test it uses, and the issues with degrees of freedom discussed below)

car::Anova(fm1, test.statistic="Chisq") # Wald test. review our discussion of how appropriate this might be for MLE (as opposed to LRT)
car::Anova(fm1,test.statistic="F") # Also see issues below about degrees of freedom

# For random effects it can be even trickier (because df are even harder to really know).

# Let's determine whether it is worth keeping Days in (1 + Days| Subject).
# First we will do it on the basis of a Likelihood Ratio Test (LRT). 
# Since we are comparing two models which only differ in random effects we will stick with REML (REstricted MaxLik), which lmer uses by default


# We can start by performing a LRT between the model above (fm1), and a model without this effect

fm0 <- lmer(Reaction ~ 1 + Days + (1|Subject), sleepstudy)
summary(fm0)
LR.model <-  -as.numeric(REMLcrit(fm1) - REMLcrit(fm0)) #REMLcrit is the deviance for the REML fit
LR.model

# FYI remember that the deviance is -2logLik so the following is equivalent
# LR.model.alt <- c(2*(logLik(fm1)-logLik(fm0)))

# So we have computed the likelihood ratio between the two models, but to compare this value to the chi-square distribution we need to know the number of parameters that differ between the models (df). So where should we start?

# If you look at the logLik for the first model fm1
logLik(fm1)
# It tells you that you have 6 parameters. Two fixed effects (intercept and slope of Days), and four random effects ( subject level variation for intercept, slope and the covariance between those parameters, as well as residual variation)

# you may wish to look at the output for the model to confirm this yourself
print(fm1)

# you do the same for the simpler model
logLik(fm0)
# you have only 4 parameters since you are not estimating a subject level variance for the Slope of Days (and thus no covariance term either)

# Thus from one point of view these two models differ by 2 df. So we can use this for our LRT

pchisq(q = LR.model, df=2, lower=F) 
# 2 df because we have to estimate variance in Days and covariance between Days and intercept


# There are two issues with this approach.
#As we discussed a while back in class, the  LRT is not trustworthy near a boundary (like zero), which is essentially what we are asking here (is the variance for Days associated with subjects = 0 )

# The second issue is how many parameters are different. If we consider just the two additional parameters we estimated (think about why there are two, we have estimated more than just the variance in Days associated with subjects), then df = 2.
# But there is the "hidden" df of used for each random effect (in this case the blups associated with Days|Subject). These are the 18 subjects. So we have estimated somewhere between 2 and 20 additional parameters (in some sense), but we do not know where. Again this is because a mixed model can be considered somewhere in between the no-pooling and complete pooling models.

# We could (I guess ) try the 20 additional parameters approach

pchisq(q = LR.model, df= 20, lower=F)  # Which has a huge effect on the test, and is almost certainly way too conservative.

# So this is where simulations can come in really handy.

# We can use the simulate function (which simulates response variables from the fitted model so it is a PARAMETRIC bootstrap). We simulate under our "lower" model (which is fm0)

S1 <- simulate(fm0, nsim = 2)
dim(S1)  # So each simulated data set represents a column
# ONE WARNING: it is not clear to me if simulate is allowing sigma^2 to vary (let alone other variances), so this may not be sampling quite correctly. You may be able to  use the sim() in arm, but it may be more straight forward to modify the code we wrote for the simulation section.  YOU CAN USE THIS approach for class though!

# We can write a function like this
LikRatioSim <- function(mod=fm0) {
	y.sim <- simulate(mod, nsim = 1)  # one set of simulated data
	model.lower <- lmer(y.sim$sim_1 ~ 1 + Days + (1|Subject), sleepstudy)
	model.full  <- lmer(y.sim$sim_1 ~ 1 + Days + (1+ Days |Subject), sleepstudy)
	LRSim <-  -as.numeric(REMLcrit(model.full) - REMLcrit(model.lower))
	return(LRSim)
	rm(y.sim, model.lower, model.full, LRSim)
}

n.sim = 100
LikRatParBoot <- replicate(n = n.sim, LikRatioSim())

# And can compare this to the observed LR
(length(LikRatParBoot[LikRatParBoot > LR.model])+1)/n.sim

# Same thing (to get p value), just cleaner (and more clever, thanks Megan!)
mean(c((LikRatParBoot > LR.model), 1))

# while I was unaware of it when I wrote this function there is also the refit() function in lmer
# see ?"simulate-mer" for an example that is very similar to the one above.


# Think about what is going on, and why this can help/work. 


# We can use a similar method to get approximate confidence intervals on the variances as well. 
# I would probably say non-parametric bootstrapping may be useful in this case, but care must be taken on HOW we would do the sampling (and do we really have enough subjects for effective resampling).  

# Permutation tests are possible, but you would have to decide how to do the shuffling (since it is only one of the effects that we want a "null" distribution for).


fm1.alt <- lmer(Reaction ~ 1 + Days + (1|Subject) + (0+ Days|Subject), sleepstudy)


### R^2 for mixed models
# One common question is how to look at the coefficient of determination (R^2) for the mixed model. Like we have with this issue for Bayesian analysis (Gelman 2006, also look in Gelman and Hill pages 473-476), this is not really well worked out and is an active area of research.  One VERY APPROXIMATE approach, for which I have NO IDEA how well (or poorly) it behaves is to examine the response variable VS the fitted values (so use at your own peril), although it almost certainly over-estimates the amount of variance accounted for by the model!

# As I pointed out in the previous script and screencast. There is now some functionality in the MuMIn library

plot(sleepstudy$Reaction, fitted(fm1))

# And we can also examined the coefficient of determination from the linear model of the original response variable onto the fitted values.
summary(lm(sleepstudy$Reaction ~ fitted(fm1)))


r.squaredGLMM(fm1) # mariginal (fixed) and condition (fixed + random R2 measures for model)
### NOTE THIS r.squaredGLMM() IS A VERY NEW METHOD AND FUNCTION!! NOT CLEAR HOW ROBUST IT IS.

# We can also try fitting and making inferences for this model using MCMC simulations from the posterior instead!

require(MCMCglmm)

prior = list(
             R=list(V=1,nu=-2),
             G=list(G1=list(V=diag(2), nu=0.02)))

# Let's go over this on the board.
             
fm1.MCMC <- MCMCglmm(Reaction ~ Days, 
   random=~us(1+Days):Subject, 
   data=sleepstudy, 
   burnin=5000, nitt=150000, thin=50, verbose=F, 
   prior=prior
   )
   
summary(fm1.MCMC$VCV) # Just for variance components
summary(fm1.MCMC$Sol) # just for fixed effects 
summary(fm1.MCMC)


plot(fm1.MCMC$VCV)
acf(fm1.MCMC$VCV)

summary(fm1.MCMC$VCV)[[1]][1]  # extracting specific values from the VCV
summary(fm1.MCMC$Sol)[[1]][1]

# Note their is (currently) no vcov method associated with MCMCglmm objects, so you need to extract them like above or use something like
colMeans(fm1.MCMC$VCV) # to get the means from the posterior distribution for each parameter.

# Let's compare

# The fixed (location) effects are pretty similar, as are the co-variances between the fixed parameters
cov2cor(vcov(fm1))

# How about for the random effects. Clearly not as spot on, but how different really??



# we can also fit the model, fixing the covariances between Days:Subject and  1:Subject to 0

fm1.MCMC.alt <- MCMCglmm(Reaction ~ Days, 
   random=~idh(1+Days):Subject, data=sleepstudy, 
   burnin=2000, nitt=90000, thin=20, verbose=F, prior=prior
   )
summary(fm1.MCMC.alt)
# to compare these models we can use a Bayesian variant of *IC, called deviance information criteron (DIC), which is essentially computed by calculating the deviance of the model at each MCMC iteration, as well as the deviance at the mean values for the posterior for each estimated parameter. Same basic rules of thumb apply as with AIC & BIC. It is worth noting that DIC is influenced by the choice of prior (and here we have used a prior with a zero co-variance....)

fm1.MCMC$DIC
fm1.MCMC.alt$DIC



# hypothesis testing for mixed models - a few R packages that primarily use parametric bootstrap on lmer objects.

# glmmML (parametric bootstrap - I think it does its own model fitting.)
# RLRsim (also a parametric bootstrap)
# pbkrtest (parametric bootstrap)
# car has anova.lmer and anova.lme, i.e. Anova(lme.object). By default it is using a type Wald's test, so keep in mind our discussion of Wald (VS a pure LRT) back when we went over inferences using MLE.


# see the example in ?"simulate-mer" in the lme4 package to roll your own, using a combination of simulate() and refit(). (from glmm.wikidot)

# also lmmfit (for some simple R^2 measures for simple mixed models)
# tools in the effects library are useful for plotting.

 # Take a look here for more details and examples
# http://glmm.wikidot.com/faq
# http://glmm.wikidot.com/software
# and for some examples.
# http://glmm.wikidot.com/examples
# Introduction to Logistic Regression using our own likelihood support functions as well as the canned function glm()
# Last updated December November 22nd 2011.

####################
# inverse logit function so that you can interpret your co-efficients in term of the probability p of an observation being = 1
# p = Pr(y[i] = 1)

invlogit <- function(x){1/(1+exp(-x))}
# how about a null model
####################


# We are going to use some data looking at the relationship between wing size, and how this influences survival in the presence of predators for Drosophila.

setwd("/Users/ian/Projects_current/mantid selection/predation_wing/")
pred <- read.csv("centroid.csv", header=T)
str(pred)
summary(pred)
attach(pred)

#Pulling out the right side wing data only for the survivors and those not so lucky...
dead.survived <- subset(pred, Side != "L" & Treatment !="C")

# creating dummy variable for survival for logistic regression
dead.survived$surv <-ifelse(dead.survived$Treatment == "D",0,1) # If the flies died then label as 0, otherwise they surivived and we  give them a dummy value of 1.

# Standardizing the size variable (common for studies of natural selection)
dead.survived$size.std <- (dead.survived$Log.Centroid.Size - mean(dead.survived$Log.Centroid.Size))/sd(dead.survived$Log.Centroid.Size)

# Remember when you standardize a variable you need to interpret it in terms of standard deviations
sd(dead.survived$Log.Centroid.Size) # ~ 0.1

# thus a change of 1 unit for size.std (the standardized size variable) corresponds to a change of 0.1 log centroid size units. Centroid size is a size measure for the wing for this data set.

######## The likelihood approach
y <- dead.survived$surv
N <- 1 # Each "trial" has one observation (Bernoulli Trial). That is we are asking what is the probability of an individual surviving!!!
x <- dead.survived$size.std


# Negative log Likelihood calculator (Support function)
logistic.reg.fn.1 <- function(a,b) {
	p.pred = exp(a + b*x)/ (1 + exp( a + b*x )) # right out of lecture
	-sum(dbinom( y, size=N, prob=p.pred, log=T )) # We can plop this right into our NLL calculator.
	}
# y is the data
# N is not the sample size for the overall trials, but the Bernoulli trial (i.e. what is the probability of the individual having survived or been eaten)

# As you can see this is really just an extension of what we did for general linear models, except our "deterministic" part of the model is the inverse logit function and we use the binomial distribution instead of the normal.

###########
# Or you can write it like this... Does not change anything!
logistic.reg.fn.2 <- function(a,b) {
	p.pred = 1/ (1 + exp(-( a + b*x )))
	-sum(dbinom( y, size=N, prob=p.pred, log=T ))
	}
###############
# note you can also use the built in function plogis to save yourself from writing everything out
logistic.reg.fn.3 <- function(a,b) {
	p.pred = plogis(a + b*x)
	-sum(dbinom( y, size=N, prob=p.pred, log=T ))
	}
####################

# you can write your optimization function the same as always

require(bbmle)
log.reg.optim.1 <- mle2(logistic.reg.fn.1, start=list(a=0.5,b=0.1))
log.reg.optim.2 <- mle2(logistic.reg.fn.2, start=list(a=0.5,b=0.1))
log.reg.optim.3 <- mle2(logistic.reg.fn.3, start=list(a=0.5,b=0.1))
log.reg.optim.1
log.reg.optim.2
log.reg.optim.3


###
profile.log.reg <- profile(log.reg.optim.1)
plot(profile.log.reg)
summary(log.reg.optim.1)
confint(profile.log.reg)
AIC(log.reg.optim.1)

# let's plot this curve (yes curve, we are using a logistic function)
par(mfrow=c(1,3))
curve( (exp(coef(log.reg.optim.1)[1] + coef(log.reg.optim.1)[2]*x))/ (1 + exp( coef(log.reg.optim.1)[1] + coef(log.reg.optim.1)[2]*x)), -3, 2, ylab = "probability of survival", ylim=c(0,1))

# It is much easier to write with the plogis() or the invlogit() we wrote.
curve(plogis(coef(log.reg.optim.1)[1] + coef(log.reg.optim.1)[2]*x),-3, 2, ylab = "probability of survival", ylim=c(0,1))

# Gee, this certainly does not look like a curve, why?
curve(plogis(coef(log.reg.optim.1)[1] + coef(log.reg.optim.1)[2]*x),xlim=c(-50, 50), ylab = "probability of survival",ylim=c(0,1))


# Look at the co-efficients.. very small



# Using glm() 
selection.glm <- glm(surv ~ size.std, family=binomial(link="logit"), data=dead.survived)
# If you look in ?glm or ?family you will see what choices you have for the glm(). 
summary(selection.glm)
# Explain residual deviance null deviance.

coef(selection.glm)
AIC(selection.glm)
anova(selection.glm, test="Chisq")
logLik(selection.glm)

par(mfrow=c(1,2))
curve( coef(selection.glm)[1] + coef(selection.glm)[2]*x, -3, 2, ylab= " What should be on this axis")

curve( coef(selection.glm)[1] + coef(selection.glm)[2]*x, -30, 30, ylab= " What should be on this axis")

#  We are looking at the log odds, not the probability 
# ylab= "log odds"
par(mfrow=c(1,2))
curve( coef(selection.glm)[1] + coef(selection.glm)[2]*x, -3, 2, ylab= "log odds")

curve( coef(selection.glm)[1] + coef(selection.glm)[2]*x, -60, 60, ylab= "log odds")

# On the log odds scale it will be a linear relationship


# let's use the inverse logit function (or plogis)
par(mfrow=c(1,2))
curve( plogis(coef(selection.glm)[1] + coef(selection.glm)[2]*x), -30, 30, ylab= " Survival probability")

curve( invlogit(coef(selection.glm)[1] + coef(selection.glm)[2]*x), -60, 60, ylab= " Survival probability")


selection.glm.null <- glm(surv ~ 1, family=binomial(link="logit"), data=dead.survived)
AIC(selection.glm.null)
logLik(selection.glm.null)

anova(selection.glm.null, selection.glm, test="Chisq")



# For comparison, let's take a look at the model if we performed a simple linear regression instead of a logistic regression.

selection.lm.1 <- glm(surv ~ size.std, family=gaussian, data=dead.survived) # Wait why not use lm()?
 # no reason, just demonstrating that you can do it with glm

# Same model with lm
selection.lm.2 <- lm(surv ~ size.std, data=dead.survived)
summary(selection.lm.1)
summary(selection.lm.2)
AIC(selection.lm.1)  # linear regression model
AIC(selection.glm)   # logistic regression model.
# Despite the fact that the overall model is not a good fit, the logistic regression in overall a better fit.


# let's compare the fit by looking at the parameter estimates
plot(jitter(dead.survived$surv, factor=0.20)~ dead.survived$size.std, xlim=c(-3,2), ylim=c(0,1)) # fit of raw data


curve( plogis(coef(selection.glm)[1] + coef(selection.glm)[2]*x), -3, 2, ylab= "Prob of survival", lwd=5,add=T)

curve(coef(selection.lm.1)[1] + coef(selection.lm.2)[2]*x,-3, 2, ylab= "Prob of survival", add=T, col="red", lwd=2)

# almost identical fits.. Not surprising for these very small parameters.



 #####  Let's look at the more realistic model including a quadratic
 selection.glm.quad <- glm(surv ~ size.std + I(size.std^2), family=binomial(link="logit"), data=dead.survived)
 summary(selection.glm.quad)
 anova(selection.glm,selection.glm.quad, test= "Chisq")
 
 # And the linear regression model for comparison?
 selection.lm.quad <- lm(surv ~ size.std + I(size.std^2),  data=dead.survived)
 summary(selection.lm.quad)
 AIC(selection.glm.quad, selection.lm.quad, selection.glm, selection.lm.1)
 
 # What does this tell us?
 plot(dead.survived$surv~ dead.survived$size.std, xlim=c(-3,2), ylim=c(0,1))
 
curve( plogis(coef(selection.glm.quad)[1] + coef(selection.glm.quad)[2]*x + coef(selection.glm.quad)[3]*I(x^2)), -3, 2, ylab= "Prob of survival", lwd=5,add=T)


curve(coef(selection.lm.quad)[1] + coef(selection.lm.quad)[2]*x + coef(selection.lm.quad)[3]*I(x^2),-3, 2, ylab= "Prob of survival", add=T, col="red", lwd=2)

# Not really much of a difference....


# For confidence intervals, again you have the obvious choices.  the asymptotic normal intervals, simulation, bootstrap, Likelihood profiles and credible intervals. If you are plotting your fitted function, it is worth also plotting the intervals back onto it like we did earlier.

# For Bayesian options
require(MCMCpack)
bayes.logistic.1 <- MCMClogit(surv ~ size.std + I(size.std^2),data=dead.survived, burnin=3000, mcmc=100000, thin=10)
summary(bayes.logistic.1)
HPDinterval(bayes.logistic.1)
acf(bayes.logistic.1) # auto-correlation to be dealt with.
plot(bayes.logistic.1)
#OR

#arm
require(arm)
bayes.logistic.2 <- bayesglm(surv ~ size.std + I(size.std^2),data=dead.survived, family=binomial(link="logit"))
mbl.2.sim <- sim(bayes.logistic.2, n.sims=30000)
summary(bayes.logistic.2)
confint(bayes.logistic.2)

coef.sim <- coef(mbl.2.sim)

#also see polr() and bayespolr() for ordered multinomial logit models.

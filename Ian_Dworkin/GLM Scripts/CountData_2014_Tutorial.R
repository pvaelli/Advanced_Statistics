# Count data
library(MASS) # for glm.nb
library(bbmle)
library(emdbook) # for the dzin function

table(quine$Days)

quine.pois1 <- glm(Days ~ Sex, data = quine, family="poisson")
summary(quine.pois1)
-2*logLik(quine.pois1)

# Should we worry about over-dispersion?


# For the likelihood, we would write it like

#design matrix
X <- model.matrix(~ Sex, data=quine)

quine.pois.NLL <- function(a,b, y=quine$Days){
	 det.model <-exp(a + b*X[,2])
	 -sum(dpois(y, lambda=det.model, log=T)) 
}

quine.pois.MLE <- mle2(quine.pois.NLL, start=list(a=2, b=0.16))

summary(quine.pois.MLE)
prof.quine.pois <- profile(quine.pois.MLE)
plot(prof.quine.pois)
confint(prof.quine.pois)


# If we wanted to use a quasi-likelihood, we could do this with glm
quine.QuassiPois1 <- glm(Days ~ Sex, data = quine, family="quasipoisson")
summary(quine.QuassiPois1)

# But we will not really discuss this here....



# Instead let us utilize a Negative binomial errors approach, to allow for over-dispersion

# Here we assume that y is distributed according to a poisson, where lambda itself varies according to a gamma distribution. Using the pre-built function in MASS, glm.nb

quine.nb1 <- glm.nb(Days ~ Sex, data = quine)
summary(quine.nb1)
plot(quine.nb1)
car::Anova(quine.nb1)



# Fitting the negative binomial model as a full MLE

# we can use dnbinom() , using the "mu"  and "size" (over-dispersion) parameterization (see Bolker page 124 )

quine.nbinom.NLL <- function(a,b, dispersion, y=quine$Days){
	 det.model <- exp(a + b*X[,2])
	 -sum(dnbinom(y, mu=det.model, size=dispersion, log=T)) 
}

#Confirm we get get the same deviance as glm.nb
-2*quine.nbinom.NLL(a=2.7229, b=0.1649, dispersion=1.0741)


# Now we fit the model using mle2, and compare to glm.nb
quine.nb.MLE <- mle2(quine.nbinom.NLL , start=list(a=2, b=0.16, dispersion=1))
summary(quine.nb.MLE)


prof.quine.zinb <- profile(quine.nb.MLE)
plot(prof.quine.zinb)
confint(prof.quine.zinb)

############
# While it will just repeat what we have just done, we can also do this by using the poisson and gamma distributions (page 124 in Bolker).
# Again, we can get a negative binomial as a result of poisson sampling, where lambda varies according to a gamma distribution
# lambda is ~gamma with shape  = dispersion parameter (k or size), and mean = mu
# your random variable (in this example quine$Days) is then poisson distributed with mean lambda
################



# Zero inflated distributions (binomial, poisson, neg. binom)
# One other way of considering this data is as a mix of a two distributions, one all zero's, and one with some positive mean, but including some zeroes. IN other words the negative binomial distribution is inflated with zeroes. 



# We can deal with this using a zero inflated distribution (page 285)
# Ben Bolker has already done the work for us (in the emdbook package)

emdbook::dzinbinom

# So let us modify our support function to allow for this.
quine.zinbinom.NLL <- function(a,b, dispersion, zprob, y=quine$Days){
	 det.model <- exp(a + b*X[,2])
	 -sum(dzinbinom(y, mu=det.model, size=dispersion, zprob=zprob, log=T)) 
}

# First we double check that we get the same result as the negative binomial when zprob=0
-2*quine.zinbinom.NLL(a=2.7229, b=0.1649, dispersion=1.0741, zprob=0)


# Now we can fit the model. 

# Since I do not know what a good starting point for zprob would be I compute the proportion of zeroes
with(quine, length(Days[Days==0])/length(Days))
# And use this as a starting value

# I will also use the MLE from the negative binomial model as  starting values for all of the other parameters

quine.zinb.MLE <- mle2(quine.zinbinom.NLL , start=list(a=2.72, b=0.1649, dispersion=1.07, zprob=0.06 ))
summary(quine.zinb.MLE)
prof.quine.zinb <- profile(quine.zinb.MLE)
plot(prof.quine.zinb)
confint(prof.quine.zinb)


# How do these models compare to each other?
AICctab(quine.zinb.MLE,quine.nb.MLE, quine.pois.MLE, 
        nobs=length(quine$Days), weights=T, base=T )


# Is the differences in AIC just a function of the additional parameters (i.e no change in the deviance)
deviance(quine.zinb.MLE)
deviance(quine.nb.MLE)

# Of course you may want to examine predicted vs observed to get a sense of whether there is a good global fit!

# There is also an R library which has many useful functions called pscl that includes various zero inflated distributions. I have not used it, but it looks very promising. The PDF associated with the library is available on ANGEL.
# Also see ZIGP, gamlss, and flexmix. MCMCglmm and glmmADMB for mixed models.

crapMLE <- function(a,b, sigma, y=quine$Days){
	 det.model <-a + b*X[,2]
	 -sum(dnorm(y, mean=det.model, sd=sigma, log=T)) 
}

normal_MLE_fit <- mle2(crapMLE, start=list(a=0, b=0, sigma=1.34))


AICctab(quine.zinb.MLE,quine.nb.MLE, quine.pois.MLE, normal_MLE_fit,
        nobs=length(quine$Days), weights=T, base=T )


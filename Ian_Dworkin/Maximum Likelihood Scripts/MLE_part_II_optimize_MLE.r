# Introduction to maximum likelihood estimation - PART II

# Written by Ian Dworkin,  modified October 14th/2010

#(using the functions of Ben Bolker  http://www.zoo.ufl.edu/bolker/)
# BOOK: Ecological models and data in R. 2008 PUP.

#  setwd("C:/Documents and Settings/Ian Dworkin/My Documents/Work/R-Project/Dll data")  # WINDOWS
setwd("/Users/ian/R/R scripts/Dll data/")   # MAC
dll.data = read.csv("dll.csv", header=TRUE) 
dll.data = na.omit(dll.data)


# not sure if I need these, but...
# dll.data$temp <- factor(dll.data$temp)
# dll.data$replicate <- factor(dll.data$replicate)

# let's look at the distribution of the number of sex comb teeth, as well as 
# its mean and sd
mean(dll.data$SCT)
sd(dll.data$SCT)
hist(dll.data$SCT, xlim = c(5, 20))


# The logic of numeric likelihood - Brute Force
### getting a likelihood for a single parameter estimate
# Assume SCT# is distributed poisson, and we want to see the likelihood for a parameter estimate
# of lambda = 11


# Let us write a function to calculate the likelihood of a set of parameters given the data.

# Since the number of SCT bristles is discrete and is constrained to be positive, let us start with a poisson distribution with lambda = l.
# k = SCT        # the data vector. 
pois1 = function(l, k=dll.data$SCT){
  -sum(dpois(k,l,log=TRUE))     # negative log likelihood calculator, notice log=T. This allows us to add instead of multiply
  }




# trial1
# what is the likelihood of lambda = l (11 in this case)
pois1_lambda_11 = pois1(l=11, k=dll.data$SCT)
pois1_lambda_11 # outputs the neg log-likelihood

#trial2, lambda= 12
pois1_lambda_12 = pois1(l=12,k=dll.data$SCT)
pois1_lambda_12
# this would be bloody slow, try it looped.

#### NOTE: WE HAVE NOT YET PERFORMED ANY OPTIMIZATION!!!!
 
 
 # So how do we find the best estimate?
 
 # let's try a brute force approach.
 
 # with a for loop. - BRUTE FORCE!!!!  
 N= 100
 l=9 # initial value for lambda
 loglik = numeric(N)      # empty vector of length N to hold the loglik's for given parameter value
 value = numeric(N)       # empty vector of length N to hold the parameter values
 for (i in 1:N) {
   loglik[i] = pois1(l,k=dll.data$SCT); value[i] =l;  
   l = l+0.1
   }
 plot(loglik~value, ylab = '-Log Likelihood') # this gets us in the neighbourhood, but not close enough
 
 # let's try smaller steps, closer to the approximate MLE for the parameter
 N= 100
 l=10.9
 loglik = numeric(N)      # empty vector of length N to hold the loglik's for given parameter value
 value = numeric(N)       # empty vector of length N to hold the parameter values
 for (i in 1:N) {
   loglik[i] = pois1(l,k=dll.data$SCT); value[i] =l;  
   l = l+0.01
   }
 plot(loglik~value, ylab = '-Log Likelihood')
 
dat.ml.fit <- cbind(loglik, value) 
which.min(dat.ml.fit[,1])
dat.ml.fit[which.min(dat.ml.fit[,1]),] # This extracts the row with the lowest -loglik

 # ok this tells us that the MLE for SCT# is ~ 11.13 or so
 # but this approach is very approximate, and extremely laborious!!!!
 
 
############# USING R's built in optimization functions  ########

# this first part does nothing terribly interesting, just estimates mean & sd
# but it is used for illustrative purposes

# first using fitdistr    - MLE for univariate distributions
 library(MASS)
 q2 = fitdistr(SCT, densfun = "poisson")
 q2
 -logLik(q2);
 AIC(q2)
 # this provides a similar estimate to what we observed for the method of moments,
 # and our brute force approach, both in terms of parameter MLE, and the log-likelihoods themselves
  

# if you want you can check if another distribution, such as normal/gaussian
q1 = fitdistr(SCT, densfun = "normal")
q1
-logLik(q1)
AIC(q1)
mle.q1 = coef(q1)
# Also provides very similar estimates to those above.   AIC is smaller. Explain .


 # 2 parameters for the normal, one for poisson... try the negative binomial ()

q.neg.binom <- fitdistr(SCT, densfun = "negative binomial", start = list(size=50000,mu=11), method="BFGS")
q.neg.binom
AIC(q.neg.binom)

# THis is not giving as good a fit?? perhaps not optimizing well (play with size, start at low numbers and go up.)

# In any case it is clear that a normal distribution is still a better fit, why?

# it may be easier to look at the poisson & negative binomial by drawing random samples
# numbers represent MLE (except size for neg binom)
sim.neg.binom <- rnbinom(length(SCT),size=1000, mu=11.13)
sim.poisson <- rpois(length(SCT),lambda=11.13)
sim.normal <- rnorm(length(SCT),mean=11.13, sd=1.62)
par(mfrow=c(2,2))
hist(sim.normal, xlim=c(0,20), ylim=c(0,550))
hist(sim.neg.binom, xlim=c(0,20), ylim=c(0,550))
hist(sim.poisson, xlim=c(0,20), ylim=c(0,550))
hist(SCT, xlim=c(0,20), ylim=c(0,550))
# none are perfect, but normal certainly comes closest.. The poisson and neg. binom both seem to have to variances that are too great We could certainly try gamma as well...


# let's overlay the observed and the theoretical given the MLE for the normal.
par(mfrow=c(1,1))
 plot(hist(dll.data$SCT),freq=F,main="",xlab="SCT #",
     ylab="Probability density",col="darkgray", xlim=c(5,20))
n <- nrow(dll.data)
curve(dnorm(x,mean(dll.data$SCT),sd(dll.data$SCT)),lty=2,add=TRUE)

# or alternatively.. computing the density of the data, and overlaying the theoretical
plot(density(SCT,adjust=2)) #check out what happens when adjust=1 or less...
curve(dnorm(x,mean(dll.data$SCT),sd(dll.data$SCT)),lty=2,add=TRUE,col="red")


 ########################
# using the MLE(stats4 library) wrapper for the function optimize. 
# This is a general and flexible way to perform MLE.
library(stats4)

# we can do this for the poisson distribution, using the pois1 function we wrote above
# with a few small changes for ease.

pois1 = function(l){              # just has the lamdba parameter in it.
  -sum(dpois(SCT,l,log=TRUE))     # negative log likelihood calculator
  }
  
mle.pois = mle(pois1, start = list(l=11)) # starting list for parameters, in this case poisson only has lambda.
summary(mle.pois)
-logLik(mle.pois)
# lots of useful things you can do with the mle object, see ?mle
plot(profile(mle.pois)) # likelihood profile
confint(mle.pois) # confidence intervals for 
AIC(mle.pois)


# for normal distribution
normNLL1 =function(m, sd1) {
    -sum(dnorm(SCT, mean=m, sd = sd1, log=T))
       }
mle.1 = mle(normNLL1, start = list(m = 11, sd1 = 2)) # starting list has both a mean and a standard deviation
summary(mle.1)
-logLik(mle.1) # the same as we observed for fitdistr.  
plot(profile(mle.1)) # likelihood profile
confint(mle.1)
AIC(mle.1)

# You will notice that you are getting the "In dnorm(x, mean, sd, log) : NaNs produced" warning
mle.2 = mle(normNLL1, start = list(m = 11, sd1 = 2), method="L-BFGS-B", lower=0.0002)
summary(mle.2)
AIC(mle.2)
confint(mle.2)
# No difference, just no warning. All we have done is set a lower constraint for the parameter search.





# In addition there is the bbmle library by Ben bolker which has an mle function
#  (mle2 so as not to confuse it with mle in stats4) which is also a wrapper for 
# optimize, and also has some nice functions
##################
# What if we want to test if genotype (wt VS Dll) has an effect on SCT number, and 
# what their MLE's would be? (i.e. an ANOVA for genotype)

#Before that let's explore it numerically.
# let us compare the estimates by genotype

# subsets of the data

dll.SCT = dll.data[dll.data$genotype=="Dll", "SCT"]   # notice " == "
# selecting observations for the Dll genotype only
l1 = length(dll.SCT)  # length of vector of values

 # and for wild type (wt)
wt.SCT = dll.data[dll.data$genotype=="wt", "SCT"]
l2 = length(wt.SCT)

# likelihood function
dll.pois = function(l, k1){
  -sum(dpois(k1,l,log=TRUE))     # negative log likelihood calculator
  }
k1 = dll.SCT

wt.pois = function(l, k2){
  -sum(dpois(k2,l,log=TRUE))     # negative log likelihood calculator
  }
k2 = wt.SCT


# numerical optimization
# For Dll
 N= 250
 l=10.0
 loglik1 = numeric(N)      # empty vector of length N to hold the loglik's for given parameter value
 value = numeric(N)       # empty vector of length N to hold the parameter values
 for (i in 1:N) {
   loglik1[i] = dll.pois(l,k1); value[i] =l;  
   l = l+0.01
   }
 #plot(loglik1~value, ylab = '-Log Likelihood') 
 
 #For wt
  N= 250
 l=10.0
 loglik2 = numeric(N)      # empty vector of length N to hold the loglik's for given parameter value
 value = numeric(N)       # empty vector of length N to hold the parameter values
 for (i in 1:N) {
   loglik2[i] = wt.pois(l,k2); value[i] =l;  
   l = l+0.01
   }
 #plot(loglik2~value, ylab = '-Log Likelihood', col='red')
 
 par(mfrow=c(1,2))
 plot(loglik1~value, ylab = '-Log Likelihood', main = 'MLE(SCT) for Dll')
 plot(loglik2~value, ylab = '-Log Likelihood', col='red', main = 'MLE(SCT) for wt')
 
 # we can also look at the MLE's for each subset 
 require(MASS)
 q.Dll = fitdistr(dll.SCT, densfun = "poisson")
 q.Dll
s1 = -logLik(q.Dll)
 
 
 q.wt = fitdistr(wt.SCT, densfun = "poisson")
 q.wt
 s2 = -logLik(q.wt)
 
 # let's look at the changes in logLik relative to that for the MLE
 loglik1B = loglik1 - s1[1]
 loglik2B = loglik2 -  s2[1]
 par(mfrow=c(1,2))
 plot(loglik1B~value, main = "MLE(SCT) for Dll genotype", ylab = "delta logLik")
 plot(loglik2B~value, main = "MLE(SCT) for wt  genotype", ylab = "delta logLik", col = 'red')

# how about confidence intervals? 
# According to Edwards (1992) an approximate 95% confidence intervals is 2 loglik units
# on either side of the mle
 
 s1_min = s1[1] - 2
 s1_max = s1[1] + 2
 s2_min = s2[1] - 2
 s2_max = s2[1] + 2
 
 par(mfrow=c(1,2))
 plot(loglik1B~value, main = "MLE(SCT) for Dll genotype", ylab = "delta logLik")
 abline(h=2)
 plot(loglik2B~value, main = "MLE(SCT) for wt  genotype", ylab = "delta logLik", col = 'red')
 abline(h=2)
 
 # this approach is fine for quick and dirty, but there is a more formal method.
 cutoffs = c(0, qchisq(c(0.95, 0.99), 1)/2) # based on chisquare, the 95 and 99% CIs. How did I figure out df's for this distribution?
 par(mfrow=c(1,2))
 plot(loglik1B~value, main = "MLE(SCT) for Dll genotype", ylab = "delta logLik")
 abline(h = cutoffs, lty = 1:3)
 text(rep(0.5, 3), cutoffs + 0.2, c("min", "95%", "99%"))
 plot(loglik2B~value, main = "MLE(SCT) for wt  genotype", ylab = "delta logLik", col = 'red')
 abline(h = cutoffs, lty = 1:3)
 text(rep(0.5, 3), cutoffs + 0.2, c("min", "95%", "99%"))
 
 # Even at a course eyeball level, clearly the CI's do not overlap.
 
 # you could also do this, without using the delta logLik just by adding the the minimum logLik
 # at the MLE...
 
 
 # LRT, using the approximate chi2 distribution of  LR - comments about how
 # poor an approximation this is at the boundary. 
 
 
 
 # Some other packages for numerical optimization.
 
 # In addition to optim, mle (stats4), mle2(bbmle), fitdistr (MASS) there are a few more options
 
 #optimize (stats) - for one dimensional optimization.
 #costrOptim - constrained optimization (more constraints than L_BFGS-B)
 #nlminb - optimization using port routines
 # nlm - nonlinear optimization
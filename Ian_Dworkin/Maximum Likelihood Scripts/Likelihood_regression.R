# Ian Dworkin
# Last modified Nov 10th 2012

# A little script to demonstrate how writing out the likelihood function gives the correct estimate (and -logLik). This also includes performing it by MLE.

n = 100
x1 = rnorm(n, 7, 2) # continuous predictor variable
y1 = rnorm(length(x1), 2*x1 + 1, 1) # continuous dependent variable with some added error
plot(y1~x1)

model1.lm <- lm(y1~x1)
summary(model1.lm)
confint(model1.lm)  # Parametric confidence intervals for the model parameters
vcov(model1.lm) # This will give us the variance and covariances between the estimated parameters.  The sqrt of the variances for the parameter should be ~the standard errors.
#i.e.
sqrt(vcov(model1.lm)[1,1])
# should be ~equal to the SE of the first estimated parameter.


# comparing MLE of models

logLik(model1.lm) # logLikelihood of model using lm()
AIC(model1.lm)
BIC(model1.lm)

 # now add a function to perform the analysis by MLE
 
 
# using mle2() in the bbmle library. mle2() is a wrapper for optim()
require(bbmle)
 
 # Can either write the likelihood calculator, then optimize...
 
linregfun1 = function(a, b, sigma) { 
      Y.pred = a + b * x1 
      -sum(dnorm(y1, mean = Y.pred, sd = sigma, log = TRUE)) } 
 
 
mle2.model <- mle2(linregfun1, 
    start= list(a= 14, b=0, sigma=1))
    
   # We get an error, although it produces the correct estimate.
warnings()    
   #This has to do with "testing" values that are impossible
  
summary(mle2.model)  # gives parameters, -2logLik etc..

# We can get the -logLik, deviance and AIC easily enough.
-logLik(mle2.model)  # in case you do not like to divide by 2 ;)
deviance(mle2.model) 
AIC(mle2.model) 

AIC(model1.lm) # for a point of comparison



### Confidence intervals
# By default this will compute the profile confidence intervals
confint(mle2.model)
 
# For complex problems, profiling takes time. It is Better to do this!
profile.mle2.model <- profile(mle2.model)
confint(profile.mle2.model)

# Plotting the profile can be very useful to see if there is anything wonky going on.
par(mfrow=c(1,3))
plot(profile.mle2.model, abs=T, conf = c(99, 95, 90, 80, 50)/100)  

  
vcov(mle2.model) # As with an lm object produces variances and covariances of the parameter estimates (but from the fisher information matrix)

sqrt(vcov(mle2.model)[1,1]) # SE of parameters.  
# Why are they lower than the OLS? Think about the MLE of the parameter variance.
 
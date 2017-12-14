# November 1st 2011
require(bbmle)


# Thinking about optimization....

setwd("/Users/ian/R/R scripts/Dll data/")   
dll.data = read.csv("dll.csv", header=TRUE) 
dll.data = na.omit(dll.data)


# In the last class, we used maximum likelihood to estimate parameters for a linear model.

# Support function
  linregfun1 = function(a=0, b=0, sigma=0.5, x1=dll.data$tarsus, y1= dll.data$SCT) { 
 Y.fixed = a + b*x1 
 -sum(dnorm(y1, mean = Y.fixed, sd = sigma, log = TRUE)) 
 } 
 


mle2.model.1  <- mle2(linregfun1, start= list(a=0,b=0,sigma=1), method="BFGS")
summary(mle2.model.1)

# We get warnings, and no estimates so the search did not do a good job of optimization. This is most obvious if you try profiling the object
plot(profile(mle2.model.1))

# As we talked about in groups in class, much of this is due to starting values. If we provide some reasonable values, things will work.
# We can use the overall mean of SCT for a starting values for the intercept, and we can use the overall sd(SCT) for sigma
mle2.model.take.2  <- mle2(linregfun1, start= list(a=mean(dll.data$SCT),b=0,sigma=sd(dll.data$SCT)), method="BFGS")
summary(mle2.model.take.2 )
plot(profile(mle2.model.take.2))

# But you can imagine there will be situations where we have some harder time getting reasonable starting values (although always start with method of moments as a first guess).
# Some like to use OLS to get starting values, which is potentially also a useful idea (if somewhat redundant)

# One thing we can do is utilize different optimization methods to get a set of "initial" estimates, and then use those estimates as starting values for a second round of optimization. 

# Let's try simulated annealing (SANN), also known as the Metropolis algorithm.

# Using different algorithms for optimization.
mle2.model.SANN  <- mle2(linregfun1, start= list(a=0,b=0,sigma=0.5), method="SANN")
summary(mle2.model.SANN) # Still somewhat off, but worth trying
plot(profile(mle2.model.SANN)) # instead of plotting it, it spits out an error, that let's you know that they found a better fit while profiling...

# let's try those values now...
mle2.model.take.3  <- mle2(linregfun1, start= list(a=7.69,b=18.42,sigma=1.56), method="L-BFGS-B", lower=c(a=-10, b=-10, sigma=0.0001))
summary(mle2.model.take.3 )
plot(profile(mle2.model.take.3))

### In general if you are having problems, I recommend
# A) Find the method of moments estimators, and using those as guesses for starting values for both fixed and random parameters in the model
# B) If BFGS or L-BFGS-B is not working, try SANN, then use those "estimates" as starting values for an optimization using BFGS.


 
# Since we are estimating 3 parameters, the Nelder-Mead Simplex may work as well.

mle2.model.Nelder  <- mle2(linregfun1, start= list(a=0,b=0,sigma=0.5), method="Nelder-Mead")
summary(mle2.model.Nelder) # nope
plot(profile(mle2.model.Nelder))
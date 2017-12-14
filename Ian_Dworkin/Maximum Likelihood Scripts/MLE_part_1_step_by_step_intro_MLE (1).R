# Ian Dworkin last updated October 23rd/2014
# this tutorial provides a step by step introduction to what likelihood is, how we compute it, and then how we find the maximum likelihood

# The concept of likelihood.

# Let's start with a set of observations Y (y1,y2,..., yN)
# In probability we generally think about the problem of having some random variable Y which is unknown, but some some probability function (like a normal distribution), where the parameters of the model (such as the mean and variance) are treated as "fixed" or "known"

# Likelihood turns this idea around. Instead we think of the observed data Y as "known" but the parameters of the model as unknown. The basic principle of Maximum Likelihood estimation (MLE) is that we find estimates for the parameters (such as the mean and variance) that maximizes the likelihood given the observations.

# We can think about this in terms of probability. Given a sample of observations Y (y1, y2,...yN)  where each element of Y is independent of each other,  we need to find a solution that maximizes the likelihood of the joint probability function f(Y;parameters.). The estimates of the parameter values that satisfy this (maximization) are called the "Maximum Likelihood Estimators".


# There are two major points. First  L(parameters|Y) is proportional to p(Y|parameters), which as we will see below makes this tractable, since we know how to find the probabilities. Second we are not interested in the absolute likelihood or probability, but the relative probability. That is we want to find the set of parameters that provides the highest probability of the observations, given the observations themselves. 



# So how do we go about doing this?


# let's take a look at the first 10 observations of a data set for femur lengths from the fruit-fly Drosophila melanogaster
femur.sample <- c(0.5903, 0.5504, 0.5884, 0.5956, 0.5767, 0.6183, 0.5817, 0.5725, 0.5680, 0.5554)
mean(femur.sample)
min(femur.sample)
max(femur.sample)


# Maximum likelihood starts with the principle that we should compare the relative probabilities of the data, given a set of parameters.


# But let's start by examining the likelihood for particular sets of parameters , given our data.

# let us assume that  we think that we can describe these 10 observations by a normal distribution. Given a set of parameter values for the normal distribution ( a "mean" and a "variance") what will be the probability of getting these 10 observations.

# let us first write out the normal distribution (we will see that there are functions to take of this later)

normal <- function(mean, variance,x){
	(1/(sqrt(2*pi*variance)))*exp(-(x-mean)^2/(2*variance))
	}
	
# Let us arbitrarily pick a number to be our first guess at the mean of our sample, I chose the min(femur.sample) which equals 0.55. Let's also assume that the variance=1 (this value is not a good estimate, but for the moment that does not matter). We will see that our initial guess was not a good one, but we can  improve upon this using methods to find the MLE.

 # Given the mean = 0.55, and the variance=sd=1, what is the "probability" of an observation from our sample being 0.59 (which happens to be the first observation in our sample). (Note: since this is a probability density we are examining the height of the density at this point, not the area under the curve, so it is approximate, but proportional to the probability in a region around 0.59. Again, we are only interested in the relative probabilities, so this is sufficient).
 
normal(mean=0.55, variance=1, x=0.59) # this is not an efficient way to code this, but it works.
 # this gives us a probability density of ~0.399
 
 # It turns out that there is a pre-built function in R to do this for the normal (and most other) distribution
 
dnorm(x=0.59, mean=0.55, sd=1) # See it gives the same number, so let's use this from here on out.
 
 
 # Let us now ask what the probability of observing 0.55 (the second observation from our sample.... i.e. femur.sample[2]    )
dnorm(x=0.55, mean=0.55, sd=1)
 
 # We of course are not particularly interested in the probability of observing each of these individual observations, but we are (to start with)
 # interested in the probability of observing both of them given our "estimated" mean and variance.
 
 # Since we treat each observation as independent, the probability of observing both of them is found by multipling each of the individual probabilities (joint probabilities).
 dnorm(x=0.59, mean=0.55, sd=1) * dnorm(x=0.55, mean=0.55, sd=1)
 
 # how about if we wanted to know the probability of all 10 of our observations in our sample. We can use an R shortcut...
 
 sample.prob <- dnorm(x=femur.sample, mean=0.55, sd=1)
 sample.prob  # the individual probabilties for each of the 10 observations in the sample
 
 # now we can multiply them together
 prod(sample.prob) # this is the joint probability of observing our sample. Which is the likelihood!!!!
 
 # We know that the values we guessed for for the mean and standard deviation were probably not so good. Let's calculate them the old fashion way (method of moments).
 mean(femur.sample) # ~0.58
 sd(femur.sample)   # ~0.02
 
 # let's input these
 sample.prob.2 <- dnorm(x=femur.sample, mean=0.58, sd=0.02)
 prod(sample.prob.2)  # the probabilty of observing the data seems better with these estimates...
 
 
 
# congratulations you have now calculated the likelihood of your data given estimates of your parameters.
  
  
  
  
  ###### HOW MIGHT WE  IMPROVE OUR FIT... STARTING TO THINK ABOUT  FINDING A MAXIMUM ######
 # We now move on to thinking about how we might calculate the MAXIMUM likelihood estimates. Where we want to optimize the likelihood given the criteria of a maximum relative probability, which means the smallest -log likelihood.
 
# Let's take a look across a range of possible value for the mean of the distribution. 


mean.x <- seq(0.55,0.61, by=0.001)

likelihood.x <- function(x) {
  sample.prob2 <- dnorm(x=femur.sample, mean=x, sd=1)
  prod(sample.prob2) 
  }
  
  
likelihood.x(x=0.55)  # making sure it gives us the same results as before
prod(sample.prob) # Our result from earlier as a comparison.

lik.1 <- sapply(mean.x, likelihood.x) # sapply just feeds in each element from the mean.x vector one at a time into the prob.x function. 
plot(lik.1 ~ mean.x, type="b", ylab="Likelihood", xlab= "values of x") # your very first likelihood plot!!!!!!
# You can see from this that the approximate MLE is around mean=0.58. The Methods of moments estimate is  also ~0.58, but this is hardly exact.

# This approach is called "grid" or "numerical" or "brute force" searching, and as we will see it is not actually very efficient.




###### REMEMEBER WE ARE ONLY INTERESTED IN THE RELATIVE LIKELIHOOD #####


# if we go back to the original function that we wrote

normal <- function(mean, variance,x){
	(1/(sqrt(2*pi*variance)))*exp(-(x-mean)^2/(2*variance))
	}

# You may notice that the first part of the expression (1/(sqrt(2*pi*variance))) does not depend upon x (our observed data).

normal2 <- function(mean, variance,x){
	exp(-(x-mean)^2/(2*variance))
	}

#For purposes of the approach we are using we can re-write this function without that part.
likelihood.x.2 <- function(mean=x, variance=1,y=femur.sample) {
  prod(exp(-(y-mean)^2/(2*variance))) 
  }

lik.2 <- sapply(mean.x, likelihood.x.2) # sapply just feeds in each element from the mean.x vector one at a time into the prob.x function. 

# Let us compare the results with and without the "constant"
par(mfrow=c(2,2))
plot(lik.1 ~ mean.x, type="b", main='Likelihood curve')
plot(lik.2 ~ mean.x, type="b", main='Likelihood curve, without the "constant"')

##### Log likelihood and Negative log Likelihoods #####

#Log Likelihood
# By convention, we log transform the individual probabilties (per observation). This allows us to add instead of multiplying, which is computationally more efficient. 
 
 log.likelihood.x <- function(x) {
  sample.prob2 <- dnorm(x=femur.sample, mean=x, sd=1, log=T) # log transformed probabilties, not the observations
  sum(sample.prob2) 
  }

log.lik.1 <- sapply(mean.x, log.likelihood.x)
plot(log.lik.1 ~ mean.x, type="b", main='Log Likelihood curve')


# Another common convention is  to use the negative log likelihood, so that the smallest value (the minimum) represents the best estimate.

neg.log.likelihood.x <- function(x) {
  sample.prob2 <- dnorm(x=femur.sample, mean=x, sd=1, log=T) # log transformed probabilties, not the observations
  -sum(sample.prob2) # Here is the negative sign
  }

neg.log.lik.1 <- sapply(mean.x, neg.log.likelihood.x)
plot(neg.log.lik.1 ~ mean.x, type="b", main='-Log Likelihood curve')

# Some people call this a "badness of fit" curve, since higher values mean worse fits.

# These all provide the same information, but by convention people use the log or negative log likelihood.

# Clearly there has to be a better way to find the maximum likelihood than this brute force approach... Which we will begin to look at in the next tutorial.

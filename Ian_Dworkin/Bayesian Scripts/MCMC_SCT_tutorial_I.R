# Last updated November 3rd 2011
# MCMC for simulating from posterior for a poisson distribution for the number of sex comb teeth in Drosophila (which is discrete and constrained to be positive)


# This is not a very efficient algorithm, but I have sacrificed efficiency for clarity ( I HOPE). In particular this is a very simple version using the Metropolis sampler algorithm.

# Currently I am doing this without an explicit function for the prior, just the likelihood.
# However that does not mean there is not a prior that is implied. Indeed, hopefully it will be clear from the lecture that implicitly I am using a flat prior, where the probability is uniform across the entire range of possible parameter values. When this is the case, the prior will divide out when calculating A. THis approach is useful for teaching purposes, however it can lead to instability in the behaviour. So it is best to explicitly write in the prior (which we will do in the next R tutorial)




# let's take a look at the first 10 observations of a data set for the number of sex-comb teeth from the fruit-fly Drosophila melanogaster
SCT.sample.2 <- c(9,13,11,12,14,11,12,10,12,13)
length(SCT.sample.2)
mean(SCT.sample.2)
min(SCT.sample.2)
max(SCT.sample.2)


######
# You may want to take a look about at the step by step introduction to MLE as well (in case you forgot)

# We are going to create a little support function "lik" to compute the likelihood
# note we are not log transforming the probabilities. However with larger sample sizes, the probabilities become vanishingly small and as far as the computer is concerned = 0. So for larger sample sizes you need to log transform the probabilities.
 
 lik <- function(x) {
 sample.prob.2 <- dpois(x=SCT.sample.2, lambda=x)
 prod(sample.prob.2) # multiplying probabilities, not adding since we did not take the log on the previous line
 }
 
 # try a few values to convince yourself that the highest probability is close to mean(SCT.sample), and the probability decreases away from this. 
 
# The likelihood we have written at the moment is the "target function"
###################

# Before we go on to do the actual computation, let's remind ourselves of the algorithm for the Metropolis sampler.
 ####
 ###FILL ME IN
 
########################

# Next, we will program a Metropolis Sampler scheme to sample
# from a distribution proportional to the target distribution 

# A is the acceptance probability for "accepting the proposed value for the parameter given the "current value".
n=10000 # # of monte carlo iterations
x = rep(0,n) # a vector with n values of 0. This replaces any previous values of x.

#Initialize with our starting value
x[1] = 16

# now if you type x, the first number in the vector  is 16, followed by 9999 0s.



# Here is the meat of the sampler, written as a for loop. Remember that the current value x[t] only depends on  what happened in the previous iteration x[t-1].


for(i in 2:n) {
  currentx = x[i-1]
  proposedx = currentx + rnorm(1,mean=0,sd=2) # the jump from the old value to the proposed new value
  ifelse(proposedx >0, proposedx,0) # constrains the proposedx (lambda) to be 0 or greater

  A = lik(proposedx)/lik(currentx) 

  if(runif(1)<A){
   x[i] = proposedx       # accept move with probability min(1,A). 
     } else {
   x[i] = currentx        # otherwise "reject" move, and stay where we are
}
}



  # Let's examine what we have done above
  
  #First we are making our initial value of x (3) our "currentx". 
  
  #We then make a proposed new value by taking our current value and adding some random value to it ( in this case from a normal distribution with mean=0 and sd=1).
  
  # We computer our acceptance probability as the ratio of the proposed values to the current value (this is A).
  
  #  We then accept the proposed value with probability min(1,A). 
  # This means that if A > 1 (the proposedx gives a higher probability of fitting the data) we accept the proposedx.
  # If A < 1, we accept the proposedx with probability A. Computationally this means we compare A to a random number between 0 and 1 (from a uniform distribution on 0 to 1).
  
  # If A is greater than  this random value then we accept the move to our proposed value, otherwise we reject it and stay where we are. 
  
  #We will always accept currentx when it is a better fit (A >=1), and we still accept it if the ratio A is greater than runif(1).
  
  # the idea here is that if we are in a low probability region of our target distribution, we will tend to accept the new 
  #  proposed values more often (and thus move). However if we are in a region of higher probability, new proposed values will    
  # tend to be only move us out of this region, so they will be rejected more often. The consequence of this approach is that 
  # more "mass" from the MCMC will be in regions of higher probabilities, approximating the posterior distribution.

# note that x is a realisation of a Markov Chain (it only depends on x[i-1])

# we can make a few plots of x:
plot(x, type="l")
hist(x,prob=T) # scales histogram so that density equals=1
acf(x) #plotting correlation between successive values of X - lag is # steps

mean(x)
sd(x)

#let us look at the first few hundred iterations to see if there if we require some burn in time.
plot(x[1:200], type="l")


# we can remove the burn in (i.e. to remove the effect of the starting values)

#  Remove the first 1000 iterations
x.burn <- x[1001:10000]
length(x.burn)
plot(x.burn, type="l")

# This can reduce the autocorrelation
acf(x.burn)
hist(x.burn)
mean(x.burn)
var(x.burn)

#Posterior mode (to compare to likelihood)
approx_mode.fun <- function(x)
{
   den <- density(x) # approximate density of posterior
   den.1 <- den[[2]] # extract heights of approximate density, double brackets [[]] extracts from list 
   max.den <- as.numeric(which.max(den.1)) # find the index for the maximum height of the density
   den.2 <- den[[1]] # extract data values corresponding to density
   den.2[max.den] # Approximate mode 
}

approx_mode.fun(x.burn)
mean(SCT.sample.2) # compare with mean. Pretty similar even without thinning.

 # exercise1 .. playing with starting values, or jump size. DO things converge to the same stationary distribution?
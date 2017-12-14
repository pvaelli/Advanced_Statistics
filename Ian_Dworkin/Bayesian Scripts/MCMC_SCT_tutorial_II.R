# Step by step introduction to MCMC
# Part II

# Here we are going to use the complete data set. This will require us to make a few changes to how we implement the algorithm.

###########
#complete data
setwd("/Users/ian/R/R scripts/Dll data/")
sct <- read.csv("dll.csv", h=T)
sct <-na.omit(sct$SCT)
mean(sct)
length(sct) # far more observations


# As sample sizes get larger, the probabilities become vanishingly small and as far as the computer is concerned = 0. So for larger sample sizes you need to log transform the probabilities.

#posterior calculator
log.post<- function(x) {
 sample.prob.2 <- dpois(x=sct, lambda=x, log=T)
 sum(sample.prob.2) + dunif(x,min=0,max=1000) #second part is the prior
 # adding log probabilities since we took the log on the previous line. Notice WE ARE NOT using the neg. log lik. just the log lik.
 }
 
 # the log.post we are using is the "target function"
 
 # A is the acceptance probability
N=40000
x = rep(0,N) # a vector with N values of 0. This replaces any previous values of x.

x[1] = 5

# now if you type x, the first number in the vector  is 5, followed by 9999 0s.



for(i in 2:N){
 currentx = x[i-1]
 proposedx = currentx + rnorm(1, mean=0, sd=0.1) # the jump from the old value to the proposed new value
 ifelse(proposedx > 0, proposedx, 0) # constrains the proposed x (lambda) to be 0 or greater
 A = exp(log.post(proposedx) - log.post(currentx)) # for using log transformed probabilities
    # you may wish to prove to your self that the above is the equivalent to initial ratio
 if(runif(1) < A){
  x[i] = proposedx       # accept move with probability min(1,A). 
   } else {
   x[i] = currentx        # otherwise "reject" move, and stay where we are
 }
}



# we can make a few plots of x:
plot(x[1:300], type="l") # how long to approach the stationary distribution
plot(x, type="l")
hist(x, xlim=c(10.5,11.5), breaks=100) 
acf(x) #plotting correlation between successive values of X - lag is # steps
acf(x, plot=F)

mean(x)
sd(x)

# we can remove the burn in (i.e. to remove the effect of the starting values). This makes a difference for virtually everything.
x.burn <- x[1001:N] # We are treating the first 10000 steps as burn in.
length(x.burn)
plot(x.burn, type="l") # This seems like it is mixing well... but look closely at say 100 values (this will depend on jump probabilities)
plot(x.burn[11000:11100], type="l") # for sd=0.1 this is not bad
acf(x.burn)
hist(x.burn)
mean(x.burn)
sd(x.burn) # the standard error of the iterations

# Exercise # How do things change when you use all of the iterations, vs burn-in after 50, 100, 500, 1000, 10000 iterations? 

# Exercise. How do things change when using different jump sizes for the proposedx?

### thinning the herd
acf(x.burn)
# If you think that the acf's are too high, but the chain is mixing well you can use a thinning approach

# here is a home-made thinning
thin <- seq(1,39000, by=10) # this produces a sequence of numbers from 1 to 30000 by values of 10
x.burn.thin <- x.burn[thin] # pulls out every 10th value
length(x.burn.thin)
acf(x.burn.thin)
mean(x.burn.thin)
plot(x.burn.thin, type="l")
plot(x.burn.thin[1:100], type="l")
hist(x.burn.thin)
var(x.burn.thin)
plot(density(x.burn.thin, adjust=1.5), xlim=c(10,14))
# in this case it did not make much of a difference

curve(dunif(x,min=10, max=14),10,14, add=T, col="red") # prior
###
# we can use some of the pre-built functions in the coda package to help look
require(coda)
 
x.mcmc <- mcmc(x.burn.thin)
plot(x.mcmc)
length(x.mcmc)
effectiveSize(x.mcmc) # uses information from chain acf to figure out how much information is actually there.

x.big.mcmc <- mcmc(x.burn)
effectiveSize(x.big.mcmc)
raftery.diag(x.mcmc)

summary(x.mcmc) # compare to our original results
HPDinterval(x.mcmc)  # confidence intervals on our estimates

#let's compare to the confidence intervals from the likelihood profile.
require(stats4)
log.lik.SCT<- function(l=7) {
 sample.prob.2 <- dpois(sct, lambda=l, log=T)
 -sum(sample.prob.2)} 
mle.sct <- mle(log.lik.SCT, start=list(l=8), method="BFGS")
confint(mle.sct)


####
#Posterior mode (to compare to likelihood)
# THis code is really ugly. Turn it into a nice function.
den <- density(x.burn.thin, adjust=1.5) # approximate density of posterior
den.1 <- unlist(den[2]) # extract heights of approximate density (unlist to remove list class)
max.den <- as.numeric(which.max(den.1)) # find the maximum point for the height of the density
den.2 <- unlist(den[1]) # extract data values corresponding to density
den.2[max.den] # Posterior mode 
mean(x.burn.thin) #posterior mean

mean(sct) # Method of moments mean. Pretty similar even without thinning.

fitdistr(sct,densfun="poisson") #MLE


### add more from my original MCMC tutorial.

 # exercise1 .. playing with starting values. DO things converge to the same stationary distribution?
 
 # Clearly some issues with auto-correlation. We could change the jump size (try 2), or use one of the better algorithms.
 
 
 #####
 ### Exercises
 ####
 
 


# Ok, so far we have done the analysis as if the the prior was completely flat, such that the probability from the prior distribution will always be a constant, which can be divided out (so if effect, it is only the likelihood that is left)


# THe first question is what is a good prior distribution for lambda?  Lambda has to be zero or positive, but we want a continuous distribution. Thoughts?

# two reasonable choices may be gamma (which is often used)  or log-normal. Which one is the conjugate prior for the Poisson distribution? Why?


curve(dgamma(x,shape=39.47, scale=0.283),0,30)  # Where do you think I got these parameter values from?

log.post<- function(x) {
 sample.prob.2 <- dpois(x=sct, lambda=x, log=T)
 sum(sample.prob.2)+ dgamma(x,shape=39.47,scale=0.283) #second part is the prior
 # adding log probabilities since we took the log on the previous line. Notice WE ARE NOT using the neg. log lik. just the log lik.
 }


#run everything as we did above
plot(density(x.burn.thin), xlim=c(5,15))

curve(dgamma(x,shape=39.47, scale=0.283),5,15, add=T, col="red") 
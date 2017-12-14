# Updated Oct 21st 2010
# demonstration of binomial likelihood, beta prior and beta posterior (conjugate pair)


# Say we come back to an example where we are modeling the probability of a fly being eaten by its predator. We set up our assays with 25 flies and 3 predators, and are interested in the survival probability P(p| D).

# A good starting place for our likelihood function is the binomial distribution, which will have N=25 (number of flies in the arena), k=# flies that survived. What we want to estimate is p= probability of an individual surviving.

# Our data may look like this
surv <- c(8, 11, 8, 10, 9, 9, 10, 10, 8, 8) # survivors per arena
N    <- c(25, 25, 25, 25, 25, 25, 25, 25, 25, 25) # total number of flies per arena
mean(surv/N) # mean of 0.364

# Our basic likelihood calculator/support function
predation.Likelihood <- function(p=0.5){
 prod(dbinom(x=surv, size=N, prob=p))}

p <- seq(0,1,by=0.001)

pred.lik <- sapply(p, predation.Likelihood)
par(mfrow=c(2,2))
plot(pred.lik~p, main="Likelihood")

lik2 <- cbind(pred.lik, p)

# Find the the parameters which give the maximum probability of the data.
lik2[which.max(lik2[,1]),]
######
# we can also use mle2() in bbmle to find the MLE for p
require(bbmle)
predation.NLL <- function(p=0.5){
 -sum(dbinom(x=surv, size=N, prob=p, log=T))} # We use the Negative log likelihood instead, as the optimizer does not do well without it.
mle.p <- mle2(predation.NLL, start=list(p=0.5))
mle.p
# plot(profile(mle.p))
##########

###But we want to take a Bayesian approach to the analysis.


#To calculate
# Let's start with the following prior distribution
curve(dbeta(x,shape1=60, shape2=60),x=c(0,1), main="beta dist a=b=60, highly subjective prior") ## This is a VERY informative prior distribution

# This is our un-normalized "posterior" calculator, which is just the prior*likelihood,
predation.unnorm.post <- function(p){
	prod( prod(dbinom(x=surv, size=N, prob=p)) * dbeta(p, shape1=60, shape2=60))} 
	
# We can also just use the likelihood function from above directly
predation.unnorm.post.alt <- function(p){
	prod( predation.Likelihood(p) * dbeta(p, shape1=60, shape2=60))} 
	
	
# Calculating the un-normalized posterior probability for each value of p
pred.posterior.unnorm <- sapply(p, predation.unnorm.post.alt)
plot(pred.posterior.unnorm ~ p, main="Un-normalized posterior distribution") 

# here we can calculate the p(D) easily, as the sum of all of the posterior probabilties at each data point
pred.posterior <- pred.posterior.unnorm/sum(pred.posterior.unnorm) 

#plotting the posterior probability distribution
plot(pred.posterior ~ p, main="posterior distribution")
pred.post.2 <- cbind(pred.posterior, p) 
pred.post.2[which.max(pred.post.2[,1]),]# mode ~ 0.408
abline(v=0.408)
# Even with a highly informative/subjective prior, the posterior is only affected somewhat....


#try again with a "FLAT" prior
# Flat priors are sometimes called 'weak' or 'uninformative'.
# The point being is that do not really influence the posterior as they have equal probability throughout the range of the parameter space (theta). However prior distributions depend on the scale (such as a log transformation). So what is an uninformative prior on one scale, may be informative on another.


# Flat prior with beta
par(mfrow=c(2,2)) 
plot(pred.lik~p, ylab="Likelihood or P(D|theta)", main="Likelihood")

curve(dbeta(x,shape1=1, shape2=1),x=c(0,1), xlab="p", main="beta distribution(a=b=1), flat PRIOR") # Flat / uninformative prior distribution.  Question to the class: How else may you specify the same prior in this case?

#calculator for likelihood*prior (un-normalized posterior)
predation.unnorm.post.flat <- function(p){
	prod( prod(dbinom(x=surv, size=N, prob=p)) * dbeta(p, shape1=1, shape2=1))}
	
	
pred.posterior.unnorm.flat <- sapply(p, predation.unnorm.post.flat)

pred.posterior.flat <- pred.posterior.unnorm.flat/sum(pred.posterior.unnorm.flat) # normalize the posterior

plot(pred.posterior.flat ~ p, main="Posterior Distribution") # 
pred.post.3 <- cbind(pred.posterior.flat, p) 
pred.post.3[which.max(pred.post.3[,1]),]# mode ~ 0.364
abline(v=0.364)


# As an exercise, it is absolutely worth your time playing around with the prior, to see what happens. You can also increase your sample size, and see how strong the prior influences the posterior distribution when you have just a little bit of data, versus a lot of data.




# Let's think about this more carefully for the above case

# In particular let's use the conjugate prior for a binomial likelihood (distribution of the observed data).
# As we discussed in lecture the beta distribution is the conjugate prior for the binomial. Our posterior distribution will also be a beta distribution.
# As a reminder, the great advantage of using the conjugate prior is that we do not need to perform any integration to find the posterior. 
# We use the observed data to UPDATE the parameters of the conjugate prior to find the posterior.

# beta has two shape parameters a,b

# Let us start with the two different priors
par(mfrow=c(2,2))
curve(dbeta(x,shape1=60, shape2=60),x=c(0,1), main="beta prior distributions",col="blue", ylab="")
curve(dbeta(x,shape1=1, shape2=1),x=c(0,1), col="red", add=T)

plot(pred.lik~p, main="Likelihood")# a quick plot of the likelihood

# For comparison, the "approximate posterior from earlier in the script"
plot(pred.posterior ~ p, main="posterior distribution, numerical", col="blue", ylab="scaled posterior density")
points(pred.posterior.flat ~ p, col="red", add=T ) # 

# From what I showed in class the hyper-parameters for the prior are a=b=60
# We update these to get our new values of a & b (often denoted a' and b')

# for the highly informative prior
a.prime <- 60 + sum(surv) -1
b.prime <- sum(N) - sum(surv) + 60 -1

normalizing.constant.1 <- sum(dbeta(p, shape1=a.prime, shape2=b.prime)) # to normalize the posterior

curve(dbeta(x,shape1=a.prime, shape2=b.prime)/normalizing.constant.1, x=c(0,1), xlab="p", ylab="scaled posterior density", main="Posterior distributions - using conjugacy", col="blue", lwd=2)



# For flat prior
a.prime.2 <- 1 + sum(surv) -1
b.prime.2 <- sum(N) - sum(surv) + 1 -1
normalizing.constant.2 <- sum(dbeta(p, shape1=a.prime.2, shape2=b.prime.2)) # to normalize the posterior
curve(dbeta(x,shape1=a.prime.2, shape2=b.prime.2)/normalizing.constant.2 ,x=c(0,1), col="red", add=T, lwd=2)




# Should we try a normal example, with unknown mean, but known variance?
# Gibbs sampler for simple normal mean hierarchical model
# Based on Brian Neelon's script - Nov 19, 2007  http://www.duke.edu/~neelo003/r/
# Modified by Ian Dworkin August 31st 2010
###################################

# Conjugate priors (mu~ N(m0,1/t0), sigma2<-Inv-chi(a, b))

set.seed(250)

######## 
# Data #
########
n <- 100 # number of observations
nsim <-5000
mu <-3 # So the real mean for the data is mu (3 in this case)
sigma <- 2 # THe real standard deviation is sigma (3 in this case)
sigma2 <- sigma^2 # the real variance
y<-rnorm(n,mu,sigma)

####################
# Prior Hyperparms #
####################
m0<-0
t0<-.001 # t0 is the precision which is the reciprocal of the sd in this case.
a<-b<-.001

# Let's think about what the priors look like
# The prior on the mean for our model will be ~N (mean=m0, sd=1/t0)
curve(dnorm(x, mean=m0, sd= 1/t0),-100, 100)



#########
# Inits #
#########
tau<-1
mu<-0
mtau<-a+n/2

#################
# Store Results #
#################
nsim<-15000
Tau<-Mu<-Sigma<-rep(0,nsim)

#########
# Gibbs #
#########
for (i in 1:nsim) {
  v <- 1/(tau*n + t0) 
  m <- v*(tau*sum(y) + t0*m0)
  mu <- Mu[i] <- rnorm(1,m,v)

  tau <- Tau[i] <- rgamma(1,mtau,b+sum((y-mu)^2)/2)
  Sigma[i] <- 1/tau
}
#######################
# Posterior summaries #
#######################
mean(Mu[3001:nsim])	# Discard first 1500 as burn-in
mean(Sigma[3001:nsim])
plot(501:nsim,Mu[501:nsim],type="l",col="lightgreen")
abline(h=mean(Mu[501:nsim]))

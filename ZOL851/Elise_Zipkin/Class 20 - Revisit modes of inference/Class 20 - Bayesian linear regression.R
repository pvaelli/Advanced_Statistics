#Program a linear regression using a Bayesian analysis "from scratch" in R

#Simulate the data
# y = alpha + beta * x + eps
alpha <- 0
beta <- 5

sigma <- 10
sampleSize <- 31

# Create a sequence of independent x-values 
x <- (-(sampleSize-1)/2):((sampleSize-1)/2)

# Create dependent values according to the linear predictor
eps = rnorm(n=sampleSize, mean=0, sd=sigma)
y <-  beta * x + alpha + eps

#Take a look at the data
plot(x,y, main="Test Data")

###################################################################################################

# Step 1 : Derive the likelihood function from the model

#For estimating parameters in a Bayesian analysis, we need to derive the likelihood function. 
#The likelihood is the probability (density) with which we would expect the observed data to 
#occur conditional on the parameters of the model. Given that our model
#y = alpha + beta*x + N(0,sigma) takes the parameters (alpha, beta, sigma) as an input, we 
#have to return the probability of obtaining the test data above under this 
#model. To do this, we simply calculate the difference between predictions y = alpha + beta*x 
#and the observed y, and then we have to look up the probability densities (using dnorm) for 
#such deviations to occur.

#Create a the likelihood function
#Note this is the specific likelihood for a linear regression 
#but modifying pred can allow this likelihood to take on a different linear predictor for the deterministic portion
#and modifying the density function (dnorm) can allow for a different distribution for the stochastic portion
likelihood <- function(param){
  m = param[1]
  b = param[2]
  sd = param[3]
  
  pred = m*x + b
  singlelikelihoods = dnorm(y, mean = pred, sd = sd, log = T)
  sumll = sum(singlelikelihoods)
  return(sumll)   
}


############################
#Side note - just for fun and to get an idea of how the likelihood function works
# Example: plot the likelihood profile of the slope beta while holding alpha and sigma fixed at their true values
slopevalues <- function(x){
  return(likelihood(c( x, alpha, sigma)))
}

slopelikelihoods <- lapply(seq(3, 7, by=.05), slopevalues )
plot (seq(3, 7, by=.05), slopelikelihoods , type="l", xlab = "values of slope parameter beta", ylab = "Log likelihood")

#You might have noticed that we return the logarithm of the probabilities in the likelihood function, 
#which is also the reason why we sum the probabilities of all our datapoints (the logarithm of a product 
#equals the sum of the logarithms). 

#Why do we do this? You don't have to, but it's strongly advisable because likelihoods, where a lot of 
#small probabilities are multiplied, can get ridiculously small pretty fast (something like 10 ^ -34). 
#At some stage, computer programs are getting into numerical rounding or underflow problems then. 
#So, bottom-line: when you program something with likelihoods, using logarithms is a safe bet.

###################################################################################################

#Step 2 - Specify priors for all three parameters: alpha, beta, sigma


# Prior distribution
prior <- function(param){
  m = param[1]
  b = param[2]
  sd = param[3]
  mprior = dnorm(m, sd = 5, log = T)
  bprior = dnorm(b, sd = 5, log = T)
  sdprior = dunif(sd, min=0, max=30, log = T)
  return(mprior+bprior+sdprior)
}

###################################################################################################

#Step 3 - Specify how the posterior is calculated

#The product of prior and likelihood is the actual quantity the MCMC will be working on. This function is called
#the posterior (or to be exact, it's called the posterior after it's normalized, which the MCMC will do for us). 
#Again, here we work with the sum because we work with logarithms.

posterior <- function(param){
  return (likelihood(param) + prior(param))
}

###################################################################################################

#Step 4 - Develop the MCMC algorithm

######## Metropolis algorithm ################

#Now, here comes the actual Metropolis-Hastings algorithm. One of the most frequent applications of this algorithm 
#(as in this example) is sampling from the posterior density in Bayesian statistics. In principle, however, the 
#algorithm may be used to sample from any integrable function. The aim of this algorithm is to jump around 
#in parameter space, but in a way that the probability to be at a point is proportional to the function we sample 
#from (this is usually called the target function). In our case this is the posterior defined above.

#This is achieved by:

# 1. Starting at a random parameter value
# 2. Choosing a new parameter value close to the old value based on some probability density that is called the 
     #proposal function
# 3. Jumping to this new point with a probability p(new)/p(old), where p is the target function, and p>1 means 
      #jumping as well

proposalfunction <- function(param){
  return(rnorm(3,mean = param, sd= c(0.2,0.8,0.5)))
}

run_metropolis_MCMC <- function(startvalue, iterations){
  chain = array(dim = c(iterations+1,3))
  chain[1,] = startvalue
  for (i in 1:iterations){
    proposal = proposalfunction(chain[i,])
    
    probab = exp(posterior(proposal) - posterior(chain[i,]))
    #Note:  p1/p2 = exp[ log(p1)-log(p2) ]
    
    if (runif(1) < probab){
      chain[i+1,] = proposal
    }else{
      chain[i+1,] = chain[i,]
    }
  }
  return(chain)
}


#Intitate the chain with random start values
startvalue = c(rnorm(1,0),rnorm(1,5),rnorm(1, 10))

#Set the number of iterations to run
niter = 50000

chain = run_metropolis_MCMC(startvalue, 50000)

#Specify the burin - the number of chains to discard from the beginning
burnIn = 10000
acceptance = 1-mean(duplicated(chain[-(1:burnIn),]))

#Again, working with the logarithms of the posterior might be a bit confusing at first, in particular when you look 
#at the line where the acceptance probability is calculated (probab = exp(posterior(proposal) - posterior(chain[i,]))). 
#To understand why we do this, note that p1/p2 = exp[log(p1)-log(p2)].

#The first steps of the algorithm may be biased by the initial value, and are therefore usually discarded for the 
#further analysis (burn-in time). An interesting output to look at is the acceptance rate: how often was a proposal 
#rejected by the metropolis-hastings acceptance criterion? The acceptance rate can be influenced by the proposal 
#function: generally, the closer the proposals are, the larger the acceptance rate. Very high acceptance rates, 
#however, are usually not beneficial: this means that the algorithms is "staying" at the same point, which results 
#in a suboptimal probing of the parameter space (mixing). It can be shown that acceptance rates between 20% and 30% 
#are optimal for typical applications.


###################################################################################################

#Step 5 - Pull out the results.  Look at summaries of the parameter estimates.  Plot the chains and histrograms
#of the parameter estimates.


#Make some plots to look at convergence and our estimated parameter values
par(mfrow = c(2,3))

hist(chain[-(1:burnIn),2],nclass=30, main="Posterior of alpha", xlab="True value = red line")
abline(v = mean(chain[-(1:burnIn),2]), lwd=2, col="blue")
abline(v = alpha, col="red", lwd=2 )

hist(chain[-(1:burnIn),1],nclass=30, , main="Posterior of beta", xlab="True value = red line" )
abline(v = mean(chain[-(1:burnIn),1]), lwd=2, col="blue")
abline(v = beta, col="red", lwd=2 )

hist(chain[-(1:burnIn),3],nclass=30, main="Posterior of sd", xlab="True value = red line")
abline(v = mean(chain[-(1:burnIn),3]), lwd=2, col="blue")
abline(v = sigma, col="red", lwd=2 )

plot(chain[-(1:burnIn),2], type = "l", xlab="True value = red line" , main = "Chain values of alpha" )
abline(h = alpha, col="red", lwd=2 )
abline(h = mean(chain[-(1:burnIn),2]), col="blue", lwd=2 )

plot(chain[-(1:burnIn),1], type = "l", xlab="True value = red line" , main = "Chain values of beta")
abline(h = beta, col="red", lwd=2 )
abline(h = mean(chain[-(1:burnIn),1]), col="blue", lwd=2 )

plot(chain[-(1:burnIn),3], type = "l", xlab="True value = red line" , main = "Chain values of sd" )
abline(h = sigma, col="red", lwd=2 )
abline(h = mean(chain[-(1:burnIn),3]), col="blue", lwd=2 )

#Pull out the parameter values and look at summary
alpha.est = (chain[-(1:burnIn),2])
beta.est = (chain[-(1:burnIn),1])
sigma.est =(chain[-(1:burnIn),3])

summary(alpha.est)
summary(beta.est)
summary(sigma.est)

#Compare our results to standard MLE estimates
summary(lm(y~x))



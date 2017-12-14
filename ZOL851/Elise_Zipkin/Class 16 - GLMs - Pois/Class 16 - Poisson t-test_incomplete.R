#We are interested in whether the abundance (population counts) of spotted salamanders 
#is different in arable  versus grassland habitat types.  You walk the perimeter of 30
#wetlands in each habitat type and count the number of spotted salamanders that you observe
#at every location.

#Generate an indicator for land use
n.sites <- 30
x <- gl(n = 2, k = n.sites, labels = c("arable", "grassland"))
n <- 2*n.sites

#Assume mean salamander counts in arable and grassland areas be 2 and 5 salamanders, 
#respectively. 

#Define the parameters according to an effects parameterization.
a=log(2)
b=log(5)-log(2) #effects parameterization! beta is the difference between alpha
# and our mean count for grassland (5)

#Expected density lambda is given by:
  
lambda <- exp(a+b*(as.numeric(x)-1))
      # x has levels 1 and 2, not 0 and 1

#Generate the count data, C, with Poisson noise
C <- rpois(n, lambda=lambda)        # Add Poisson noise using random Poisson function
aggregate(C, by = list(x), FUN = mean)		# The observed means


#Use a boxplot to blot the data

boxplot(C ~ x, col="grey",  xlab="Land Type", ylab="Salamander count", las=1)

#######################################################################

#Fit Poisson t-test in R

#Use the glm function and the "family" indicator
#Do ?glm and look up the different families

?glm

poisson.t.test <- glm(C ~ x, family=poisson)
summary(poisson.t.test)

alpha=poisson.t.test$coefficient[1]
beta=poisson.t.test$coefficient[2]

#Inspect the design matrix of that model
#What kind of parameterization does the glm function use?

model.matrix(poisson.t.test)


#Pull out the coefficients of the model - the parameters that we estimated
#alpha = a
#beta = b
#What is the interpretation of the coefficients?

exp(alpha) #log counts of habitat A; take exponent to obtain value
exp(alpha+beta) # beta is the log difference between counts
# not the same as exp(alpha)+exp(beta), but is the same as exp(alpha)*exp(beta)

#What is the mean estimated count for salamanders in arable habitat
#Need to transform the coefficients back to the normal scale


#How does that compare to the value we used to generate the data


#What is the mean estimated count for salamanders in grassland habitat
#Need to transform the coefficients back to the normal scale



#How does that compare to the value we used to generate the data




#### Note!
#exp(a + b) does not equal exp(a) + exp(b)
#Check this out for yourself and think about the importance of this for future analyses

## What are the confident intervals for the parameters

confint(poisson.t.test)

#Based on this and the p-values, what can you conclude about the differences in counts
#of salamanders in wetlands in arable and grassland habitats?
#Remember values are on the log scale

#Pull out the residuals

residual = poisson.t.test$residuals
predicted = poisson.t.test$fitted.values

#Is there a pattern in the residuals?
plot(predicted, residual, main = "Residuals vs. predicted values", 
     las = 1, xlab = "Predicted values", ylab = "Residuals")
abline(h = 0)

#####################################
#You can also run this model using a means parameterization

#Use the glm function and remove the intercept
poisson.t.test1 <-  glm(C ~ x-1, family=poisson) # minus 1 removes the intercept and makes this a means parameterization instead of effects parameterization
summary(poisson.t.test1)					

#Inspect the design matrix of that model

model.matrix(poisson.t.test1)

#Look at the confidence intervals of the parameter estimates.
#Can take the exp to explicitly look at the 95% CI for mean counts in the two habitat types

exp(confint(poisson.t.test1))

#######################################################################

#Fit Poisson t-test with JAGS

#Load the correct library
library("R2jags")

# Write model into R working directory
sink("Poisson.t.test.txt")
cat("
model {

 # Priors
 alpha ~ dnorm(0,0.001)
 beta ~ dnorm(0,0.001)

 # Likelihood
 for (i in 1:n) {
     C[i] ~ dpois(lambda[i]) 
     log(lambda[i]) <- alpha + beta *x[i]
 }
}
",fill=TRUE)
sink()

# Package data
win.data <- list(C = C, x = as.numeric(x)-1, n = length(x))

# Inits
inits <- function(){ list(alpha=rlnorm(1), beta=rlnorm(1))}

# Stuff to get estimates from
params <- c("alpha", "beta")
#params <- c("lambda","alpha", "beta")

# MCMC settings
nc <- 3  	# Number of chains
ni <- 3000		# Number of draws from posterior
nb <- 1500		# Number of draws to discard as burn-in
nt <- 2		# Thinning rate


out <- jags(win.data, inits, params, "Poisson.t.test.txt", n.chains = nc, 
             n.thin = nt, n.iter = ni, n.burnin = nb, 
             working.directory = getwd())


#Look at the results

jagsout <- as.mcmc.list(out$BUGSoutput)
plot(jagsout)
summary(jagsout)





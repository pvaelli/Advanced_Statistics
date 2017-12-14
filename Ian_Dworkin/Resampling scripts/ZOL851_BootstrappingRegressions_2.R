##### Using the non-parametric bootstrap for inference with linear models.

# Written by Ian Dworkin. Last updated October 27th 2012

########
# Chunk 1: Getting the Data in
# Chunk 2: The basic linear model we are examining in this tutorial
# Chunk 3: non-parametric bootstrap: Pairs (Random) approach
# Chunk 4: non-parametric bootstrap: Residual (Fixed) approach
# Chunk 5: the Monte Carlo Simulation approach (also called Parametric bootstrap), for comparison 
# Chunk 6: Comparing distributions (expected and observed) for the coefficients
# Chunk 7: Using the boot library for these models (as an alternative to my function calls) - THIS ALLOWS Bias corrected confidence intervals.
# Chunk 8: Using the simpleboot library for these models (as an alternative to my own functions)

require(car)

# Chunk 1: Getting the Data in
# Getting the data in.
setwd("/Users/ian/R/R scripts/Dll data/") 
dll.data = read.csv("dll.txt", header=TRUE)   #data frame 
dll.data = na.omit(dll.data)
dll.data$temp <- factor(dll.data$temp)
dll.data$replicate <- factor(dll.data$replicate)
str(dll.data) 
summary(dll.data)


# We are also going to define "wt" as the reference level for the genotype factor
dll.data$genotype <-relevel(dll.data$genotype, ref="wt")
str(dll.data)



# Chunk 2: The basic linear model we are examining in this tutorial
# The basic model to examine.
observ.lm <- lm(SCT ~ genotype*temp, data=dll.data)
summary(observ.lm)
confint(observ.lm)

#quick look to remind ourselves of the co-variance among parameters
cov2cor(vcov(observ.lm))
vif(observ.lm)

par(mfrow=c(4,4))
confidenceEllipse(observ.lm)


# Chunk 3: non-parametric bootstrap: Pairs (Random x) approach
##### Pairs bootstrapping for a linear regression using sample #####

N = 5000  # of Resampling iterations

# Pairs/ random variable resampling
BootstrapRandomX <- function(dat=dll.data, mod.formula=formula(SCT ~ genotype*temp)){
  dat.boot <- dat[sample(x = nrow(dat), size = nrow(dat), replace=T),] # samples along index  
  boot.lm <- lm(mod.formula, data=dat.boot)
  coef(boot.lm)}

# perform N bootstrap iterations
vector.boot <- t(replicate(N, BootstrapRandomX()))

# standard error of the estimates via bootstrap
apply(vector.boot, MARGIN = 2, sd) 


# Percentile CIs (transposed so it is oriented like the output of confint)
t(apply(vector.boot, MARGIN = 2, quantile, probs=c(0.025, 0.975)))

# asymptotic normal CIs
confint(observ.lm)

par(mfrow=c(2,2))
MultipleHistograms <- function(X=vector.boot){

    for (i in 1:ncol(X)) {
	    hist(X[,i], freq=F,
	        main = colnames(X)[i],
	        xlab = colnames(X)[i])}}

MultipleHistograms()	    


# Look at the correlations among parameter estimates from the bootstrap iterations
pairs(paired.bootstrap)
	        	    
#### The Percentile Confidence Intervals can be biased in some cases (see Efron and Tibshriani). Therefore the stat gods developed the Bias-Corrected (BC) and acceletered (a) non-parametric bootstrap confidence intervals (BCa) to adjust for these biases. At the bottom of the script I show how to compute the BCa, but in general it is easier (but slower) to use the boot() in the boot library. We will do an example with this library in a bit.



# Chunk 4: non-parametric bootstrap: Residual (Fixed) approach
# I think of this as half way between the non-parametric (pairs/random) bootstrap and the parametric (monte carlo) bootstrap. The reason is that you end up using the fitted (deterministic) component of the model, but then you combine these with the empirical residuals. It is important to remember that the assumption of exchangability among observations (in this case among residuals) is in play, unlike the pairs/random-x bootstrap.

##### The residual/ Fixed effect Bootstrap
# Here we are performing an analysis analogous to what we did for Monte Carlo simulations to generate confidence intervals.
# The steps for the residual bootstrap method
#1 - Fit model (as normal)
#2 - Extract Residuals from the model
#3 - Bootstrap the residuals from the model (r*)
#4 - add the boostrapped residuals (r*) back onto the fitted component of the model. i.e.  y*[i] ~ a + b*x[i] +r*[i] 
# using Y* refit the model, extract co-efficients
# - Repeat N times.

# For the model above (obser.lm) extract residuals
resid.model.1 <- resid(observ.lm)
plot(density(resid.model.1, bw=0.5)) # Reasonably normal looking residuals, so the results should be pretty comparable to the Monte Carlo approach.

# We may also want to check this....
par(mfrow=c(1,2))
plot(resid.model.1 ~ dll.data$genotype)
plot(resid.model.1 ~ dll.data$temp)

# Write a function to perform the bootstrapping
# To make this easy we are extracting everything from the model call to lm
BootstrapFromResiduals <- function(mod.object = observ.lm, dat = dll.data) {
	
	resids = mod.object$resid # extracts residuals from model
	fittedValues = mod.object$fitted # extracts fitted values
	matr <- model.matrix(mod.object)
	# generating new values for each y[i], by adding  bootstrapped resids to fitted values.
	Y <- fittedValues + sample(resids, length(resids), replace=T) 
	
	# Using model.matrix for the predictors (not pretty, I know)
	model.boot <- lm( Y ~ 0 + matr, data=dat ) # refit model with new Y values
	
	coef(model.boot) # Extract the co-efficients
	}

residual.boot.N <- t(replicate(N, BootstrapFromResiduals()))

# Plot
par(mfrow=c(2,2))
MultipleHistograms(X=residual.boot.N)

pairs(residual.boot.N)


t(apply(residual.boot.N, MARGIN = 2, quantile, probs=c(0.025, 0.975)))
confint(observ.lm)


# Chunk 5: the Monte Carlo Simulation approach (also called Parametric bootstrap), that we used last week.
########
#For comparison here is the code for generating confidence intervals using a Monte Carlo simulation (also called the Parametric bootstrap)

	
SimulationUnderModel <- function(model = observ.lm){
	#extract design matrix
	matr <- model.matrix(model)
	rse = summary(model)$sigma 
	df=model$df 
	
	# incorporate uncertainty in RSE
	rse.sim <- rse*sqrt(df/rchisq(1, df = df)) 
	
	# Simulate data (response) conditional on the simulated RSE.
	y.sim <- rnorm(n = nrow(matr), 
	    mean = matr%*%coef(model), sd=rse.sim)
	# 0 + design matrix (since the intercept is already in the design matrix)
	lm.sim <- lm(y.sim ~ 0 + matr) # fit model with simulated response
	coef(lm.sim)}
		
sim.coef <- t(replicate(N, SimulationUnderModel()))

t(apply(sim.coef, MARGIN = 2, quantile, probs=c(0.025, 0.975)))


# Chunk 6: Comparing distributions (expected and observed) for the coefficients
########
##### Let's compare how we do with each of the methods we have used so far.
# 1) Asymptotic normal (the standard way), #2) Monte Carlo simulation (remember the approach we are using is also sometimes called the parametric bootstrap) #3) Pairs bootstrap #4) residual bootstrap.


# We can also compare the pairs (random) and residual (fixed) bootstrapping
par(mfrow=c(1,1))
plot(density(residual.boot.N[,2], bw=0.1), main="comparing bootstrap (and other) methods for parameter uncertainty", lwd=2, lty=2) #residual resampling
lines(density(vector.boot[,2], bw=0.1), col="red", lwd=2, lty=2) # pairs bootstrap method

# From the Monte Carlo simulation last week
lines(density(sim.coef[,2], bw=0.1), col="purple", lwd=2, lty=1) 

legend(x = 0.6, y= 2.7, legend=c("Residual Boot", "Pairs Boot", "Monte Carlo Normal"), 
    col=c("black", "red", "purple"), lty=c(2,2,1), lwd=2)
######




################ Chunk 7: Using the boot library

# The code I have written works just fine, but it does not correct for small biases which are known to occur with bootstrapped estimates. The library in R that can be used for bootstrapping is called boot. While the syntax for this library can be somewhat confusing, it is an extremely powerful and flexible library.


# The R book (Crawley) gives some description of its uses ( pages 320-322, 418-421, 523 for generalized linear models & 681-683 for non-linear models). Also see the PDF on ANGEL "Fox_appendix-bootstrapping" for a fair bit of detail.

# The trick to using this library is to write a function call to get the coefficient (or whatever is of interest to you). The function call should have at least two arguments, the first being the values we want to resample (YourObjectToBeResampled), and the second argument being an index that is used by boot() to sample. It is somewhat slow though!
# Efron, B. and Tibshirani, R. (1993) An Introduction to the Bootstrap. Chapman & Hall.
require(boot)

# The basic arguments for it would look like
boot(data = YourObjectToBeResampled, statistic = YourFunction, R = NumberResamplingEvents )


# Here is an example

# YourObjectToBeResampled
BootstrapFunctionRegression <- function(data=dll.data, index) {
	data <- data[index,] # We will sample along rows of the data frame
	model.boot <- lm( SCT ~ genotype*temp, data=data)
	coef(model.boot)
}

# R is the number of resampling events, which needs to be greater than number of rows in the dataset.
bootstrappedModel <- boot(dll.data, BootstrapFunctionRegression, R = 5000)
bootstrappedModel

# There are some plotting options for the bootstrap
plot(bootstrappedModel, index = 1)


plot(bootstrappedModel, index = 4)
boot.ci(bootstrappedModel , conf = 0.95, type = c("basic", "bca", "perc"), index=2) # index for which coefficient...

confint(observ.lm)

# We can also examine the joint distriution of the intercept (i.e wild-type), and the coefficient for the mutant effect
plot(bootstrappedModel$t[,1], bootstrappedModel$t[,2], 
    xlab = "wild-type estimate", ylab= "treatment contrast for Dll", pch=16)


# This may be more informative
dataEllipse(bootstrappedModel$t[,1], bootstrappedModel$t[,2], xlab = "wild-type estimate", ylab= "treatment contrast for Dll", levels=c(0.5, 0.95, 0.99), robust=T)


#CHUNK 8 CHUNK 8 CHUNK 8 CHUNK 8: Using the simpleboot library for these models. NOTE some people have had difficulty getting this to work, so if you have trouble, just stick with my functions for now.

require(simpleboot) # need to install this first

# random-x/pairs non-parametric bootstrap
regression.boot <-lm.boot(lm.full, R=5000)
summary(regression.boot)
#plot(regression.boot)
#percentile confidence intervals
perc(regression.boot) # percentile confidence intervals

#a messy way to extract the samples for plotting
regression.boot.CI <- samples(regression.boot, name = "coef")
regression.boot.CI.slope <-regression.boot.CI[2,]
hist(regression.boot.CI.slope) # distribution of the resampled values for the slope


# or we can resample the residuals
regression.boot.res <-lm.boot(lm.full, R=1000, rows=F) #rows=F is the option to use the residual approach.
summary(regression.boot.res)
perc(regression.boot.res) # percentile confidence intervals

#a messy way to extract the samples for plotting
regression.boot.res.CI <- samples(regression.boot.res, name = "coef")
regression.boot.res.CI.slope <-regression.boot.res.CI[2,]
hist(regression.boot.res.CI.slope) # distribution of the resampled values for the slope


####BCa
# There are two pieces to this (and we will see that the boot library can take care of all of these, although it is slow). First is to adjust for lack of symmetry. 
#For this example we will do this for the slope (the second column in vector.boot)

# To correct for bias (BC)
z <- function(x=vector.boot[,2], value = coef(observ.lm)[2], N = 1000){
	x <- (sum(x <= value))/(R + 1)
	z <- qnorm(x)
	z
} 


# The acceleration (a) factor (uses the third moment of the standard normal).
# Here we use the jackknife procedure (leave one out, and re-estimate for each observation in the original sample).

jackie <- function(dat = dll.data, mod.formula=formula(SCT ~ genotype)){
	n <- nrow(dat)
	t <- rep(NA, length=n)
	for (i in 1:n){
		model.1 <- lm(mod.formula, data=dat[-i,])
		t[i] <-coef(model.1)[2]
	}
	
	t.mean <- mean(t)
	bias   <- mean(t) - t
	
	a   <- sum(bias^3) / (6*(sum(bias^2))^1.5)
	return(a) 
}
a1 <- jackie()
z1 <- z()


# Then all of this information is used to adjust what percentiles we actually use to get the 95% coverage we expect. That is what this correction factor below does.
CorrectionFactor <- function(a = a1, z = z1, alpha=0.05 ){
	fact1 <- z + 
	  ( 
	     (z + qnorm(alpha/2) )/ 
	     (1 - a*(z + qnorm((alpha/2)))) 
	  ) 
	
	fact2 <- z + 
	  ( 
	     (z + qnorm((1 -(alpha/2) )))/ 
	     (1 - a*(z + qnorm((1 -(alpha/2) )))) 
	  )   
	
	 
	a1 <- pnorm(fact1)
	a2 <- pnorm(fact2)
	return(c(a1=a1, a2=a2))    
}

CF <- CorrectionFactor()
BCa.CI <-quantile(vector.boot[,2], probs=c(CF[1],CF[2]))
BCa.CI  # Corrected 95% CIs

# THis is all a pain in the butt.. So to use this we can just use the boot library instead!

########



# Exercises to do in groups
# A) Using the pairs/random resampling approach at the top (with our home made functions). Perform the analysis , but with genotype as a covariate. 

# B) Then do the same with BOTH genotype and tarsus as covariates.

# C) Confirm these results using the simpleboot random/pairs resampling, and then with residual resampling for each set of models (~ genotype  and ~ genotype + tarsus). Compare the results? What do you think is going on? You made need to treat genotype as numeric (i.e. as.numeric(genotype)) for the simpleboot library.

# D) Last week with the Monte Carlo Simulations we plotted Monte Carlo fits for the model. Do the same for at least one of the models (I would suggest the one with just ~ tarsus as a covariate)


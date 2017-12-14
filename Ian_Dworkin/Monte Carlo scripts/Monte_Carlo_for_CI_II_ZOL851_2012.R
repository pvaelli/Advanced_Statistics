#ZOL851 - Using Monte Carlo simulations to generate confidence intervals part II - more general approaches.
# updated October 10th 2012

# In the previous tutorials we generated confidence intervals using monte carlo methods, but the method we used did not necessarily allow us to include the correlations among predictors. Below we address that using a couple of different approaches.

# Setting a few options for clarity of printing on screen.
options(digits=3, sig.digits=3, show.signif.stars=F)


#Read data into R 
#dll.data = read.csv("/Users/ian/R/R scripts/Dll data/dll.csv", header=TRUE)
dll.data = read.csv("http://datadryad.org/bitstream/handle/10255/dryad.8377/dll.csv",
   header=TRUE)
dll.data$temp <- factor(dll.data$temp)
dll.data$replicate <- factor(dll.data$replicate)
dll.data$genotype <- relevel(dll.data$genotype, "wt")  # Setting the wild-type (wt) as reference



# We will remove rows with missing data for these examples.
dll.data <- na.omit(dll.data)
# word of warning, sometimes omitting data creates variants of the object class which can cause problems in other functions (like predict as it turns out). If you are having issues check mode() and class() of your variables.

# center our covariate (for interpreting the intercept later on)
dll.data$tarsus.cent <- dll.data$tarsus - mean(dll.data$tarsus)


# Plot the data (with some jitter for SCT to improve clarity)
plot(jitter(dll.data$SCT) ~ dll.data$tarsus.cent, xlab="tarsus", ylab="SCT")

dll.model <- lm(SCT~ 1 + tarsus.cent, data=dll.data)
summary(dll.model)$coef[,1:2]
abline(dll.model, lwd=3) # This will plot the parameter estimates as a line onto the plot

# Plot the asymptotic normal Confidence intervals

tarsus.frame <- data.frame(tarsus.cent = seq(min(dll.data$tarsus.cent), max(dll.data$tarsus.cent),
   length.out=100)) 
# creating a data frame with "values" of our predictors. 
#We generate a sequence of 100 numbers (length.out), from minimum to maximum value for our covariate.

conf.int <- predict (dll.model, interval="conf", newdata=tarsus.frame) # generating confidence intervals for the 100 values of "tarsus"

#Here we plot sex comb teeth as a function of tarsus length
par(mfrow=c(1,1))
plot(dll.data$SCT~dll.data$tarsus.cent, 
    xlab="tarsus length", ylab="number of sex comb teeth", 
    main="Calculated asymptotic Normal Confidence Intervals for Observed data")

#Here we add the confidence bands to the previous plot 
matlines(x = tarsus.frame$tarsus.cent, conf.int, lty=c(1,2,2), lwd=2, col="red")

# Approximate CI's for model via monte carlo simulation, Approach simulating the response variable.
# In this version we are not simulating the covariates at all, but using OBSERVED values. 
# I think this is generally the easiest approach, as you do not need to worry about explicitly modeling co-variation among the co-variates (and thus parameter estimates), as it is already in the real data.


SimulationUnderModel <- function(model = dll.model, covariate = dll.data$tarsus.cent){
	#extract what we need from the model from which we wish to simulate
	a = coef(model)[1]
	b = coef(model)[2]
	rse = summary(model)$sigma 
	df=model$df 
	
	rse.sim <- rse*sqrt(df/rchisq(1, df = df)) # incorporate uncertainty in RSE
	# Simulate data (response) conditional on the simulated RSE.
	y.sim <- rnorm(n = length(covariate), mean = a + b*covariate, sd=rse.sim)
	lm.sim <- lm(y.sim ~ 1 + covariate) # fit model with simulated response
	coef(lm.sim)}

N = 1000

simulated.coefficients <- replicate(N,SimulationUnderModel())
simulated.coefficients <- t(simulated.coefficients)

#confidence intervals (simulated and asymptotic)
quantile(simulated.coefficients[,2], c(0.025, 0.975))
confint(dll.model)[2,]

# plot(SCT ~ tarsus.cent, data=dll.data) # if you got rid of the previous plot
for (i in 1:95) {abline(a=simulated.coefficients[i,1],b=simulated.coefficients[i,2], lwd=1, col="grey")} 

#let's take a look at the distribution of the simulated values, versus the estimates
par(mfrow=c(1,2))	
hist(simulated.coefficients[,2]) # histogram of simulated slopes
abline(v=coef(dll.model)[2], col="red", lwd=2, lty=2) # estimated slope
hist(simulated.coefficients[,1])	# histogram of simulated intercepts
abline(v=coef(dll.model)[1], col="red", lwd=2, lty=2) # estimated intercept

#examining the covariance between the parameter estimates
par(mfrow=c(2,1))
plot(simulated.coefficients[,2]~simulated.coefficients[,1], 
    main="covariation between simulated parameter estimates",
    ylab = "simulates slopes",
    xlab = "simulated intercepts") 
# scatter plot demonstrating covariation between simulated parameter values
cov(simulated.coefficients) # covariance between simulated parameter values
vcov(dll.model) # estimated covariance between parameters
# compare to simulated values also check out var(slope) & var(intercept)

# Here is the same idea but from the joint confidence intervals from the parametric estimates.
car::confidenceEllipse(dll.model, main="co-variance for the parametric estimates")


##########Explicitly modeling the covariation among predictors while constructing approximate confidence intervals using Monte Carlo Simulation
 
# (similar to  sim() in arm library) 


library(MASS)	#load library for function mvrnorm()

# one important thing to remember is that the estimated parameter values for b0 and b1 may covary. We need to use this co-variation, so we will assume a joint bivariate normal distribution. mvnnorm().

n <- nrow(dll.data) #number of observations

#Here we generate a random intercept and slope for the observed model and plot the curve.
# Think about the difference here. We are not simulated the response and refitting, but simulating from the distribution of slopes and intercepts instead.


ModelSimParameters <- function(model=dll.model) {
  
	# Here we are inputting the estimated co-efficients from the model in mu, and 
	# Sigma is the variance covariance matrix for the parameter estimates.

	par.mat <- mvrnorm(n=1, # Number of simulated values
	    mu=c( coef(model)[1], coef(model)[2]), # vector of "means", in this case mean and intercept
	    Sigma=vcov(model))   # matrix of variances and co-variances
	return(par.mat)    		
}

# run it once.
ModelSimParameters()

# Let's do N simulations

simulated_slope_intercept <-replicate(N, ModelSimParameters())
simulated_slope_intercept <- t(simulated_slope_intercept)

#confidence intervals (simulated and asymptotic)
quantile(simulated_slope_intercept[,2], probs=c(0.025, 0.975))
quantile(simulated.coefficients[,2], c(0.025, 0.975))
confint(dll.model)[2,]


#Create a scatter plot of the response as a function of the predictor
par(mfrow=c(1,1))
plot(jitter(axodata$EOG.norm) ~ axodata$GSI, 
    xlab="tarsus length centered", ylab="number of sex comb teeth", 
    main="Simulated Confidence Intervals for Observed Model", pch=16)

for (i in 1:95) {
	abline(a=simulated_slope_intercept[i,1], b=simulated_slope_intercept[i,2], lwd=1, col="grey")} 

#Here we add the confidence bands to the previous plot 
matlines(tarsus.frame$tarsus.cent, conf.int, lty=c(1,2,2),lwd=2, col="red") # just adding back the asymptotic normal confidence intervals.


########Gelman and Hill have their own function sim() to perform all of this, but ours works just as well for basic problems
# However for more complex model (GLim and mixed models) the sim() is quite useful.####
# See chapter 7 of Gelman and Hill to see the approach that they use for the simulation (it is a bit different than ours)

library(arm)

#We use the sim function to generate simulations from the observed data
dll.sim <- sim(dll.model, n.sims=N)
sim.slopes <- dll.sim@coef[,2] # extracts the simulated slopes
quantile(sim.slopes,c(0.025,0.975)) # simulated con
confint(dll.model)


#Next we again plot sex comb teeth as a function of tarsus length
par(mfrow=c(1,1))
plot(jitter(dll.data$SCT) ~ dll.data$tarsus.cent, 
    xlab="tarsus length centered", ylab="number of sex comb teeth", 
    main="Simulated Confidence Intervals for Observed Model", pch=16)

#Here we plot the curves generated from the above simulations from sim()
for (i in 1:95){
	curve (dll.sim@coef[i,1] + dll.sim@coef[i,2] * x, add=T, col="gray")
	}

#asymptotic CIs
matlines(tarsus.frame$tarsus.cent, conf.int, lty=c(1,2,2),lwd=2, col="red") 

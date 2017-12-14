#ZOL851 Tutorial - Using Monte Carlo simulations for generating distributions under "null" models
# October 9th 2012, Ian Dworkin

# In the last tutorial we developed the tools to generate approximate confidence intervals on parameters in a simple regression model, using simulated data. Now we will switch back to some observed data.

# However before we venture back into generating approximate MC confidence intervals, we are going to consider how we can use Monte Carlo simulations to generate the distributions of parameter estimates under null models. While these can be used in a variety of useful ways, we will start with the obvious application of using them to generate p-values.

#Read data into R 
#dll.data = read.csv("/Users/ian/R/R scripts/Dll data/dll.csv", header=TRUE)
dll.data = read.csv("http://datadryad.org/bitstream/handle/10255/dryad.8377/dll.csv",
   header=TRUE)
dll.data$temp <- factor(dll.data$temp)
dll.data$replicate <- factor(dll.data$replicate)
dll.data$genotype <- relevel(dll.data$genotype, "wt")  # Setting the wild-type (wt) as reference

# For the example below only need to use a subset of the columns, so we will work with those, and remove only missing data with respect to those variables.
dll.data.subset <-na.omit(dll.data[, c("SCT", "tarsus")])
str(dll.data.subset)

# center our covariate (for interpreting the intercept later on)
dll.data.subset$tarsus.cent <- scale(dll.data.subset$tarsus, scale=F, center=T)


# As a reminded of what was discussed in class, let's review what we need to do.
#1 Fit the null model with the observed data.
#2 Using the estimates from the null model, simulate the data (response values) under the null
#3 Estimate the "full" model using the simulated data under the null model.
#4 Repeat this a large number (N) times. Store results for the parameters of interest
#5 Fit the full model with the real data. Compare these estimates to the distributions generated under simulation.

SimulationUnderNull <- function(model = dll.null, covariate = dll.data.subset$tarsus.cent){
	a = coef(model)[1] # intercept for null model
	rse=summary(model)$sigma # residual squared error for null model
	df=model$df # df for the null model
	
	# First we simulate the rse based on the observed rse and the df (from a chisq distribution)
	rse.sim <- rse*sqrt(df/rchisq(1, df = df)) 	
	
	# Now we generate simulated vector of the response variable under the null
	y.sim <- rnorm(n = length(covariate), mean = a, sd=rse.sim)
	
	# we take the simulated response vector and model it as a function of our co-variate(s) of interest
	lm.sim <- lm(y.sim ~ covariate)
	coef(lm.sim)} # grab coefficients

#Fit model under null
dll.null  <- lm(SCT~ 1, data=dll.data.subset) 
dll.model <- lm(SCT~ 1+ tarsus.cent, data=dll.data.subset) # running the full model
observed.coef <- coef(dll.model)[2] # extracting the slope of the model using the coef(). 

#let's take a quick look at some important parameters.
coef(dll.null)
coef(dll.model)
summary(dll.null)$sigma
summary(dll.model)$sigma

#let's try running just a single iteration of the simulation under the null.
SimulationUnderNull(model=dll.null, covariate=dll.data.subset$tarsus.cent)

# Ok, now let's do a lot of iterations
N=1000
simulated.coefficients <- replicate(N,SimulationUnderNull())
simulated.coefficients <- t(simulated.coefficients)
hist(simulated.coefficients[,2]) # compare to the observed slope

quantile(simulated.coefficients[,2], c(0.025, 0.975)) 
# 95% percentiles UNDER THE NULL (NOTE: NOT confidence intervals for the estimates)

simulated.slopes <- simulated.coefficients[,2]
p.value.2 <- length(simulated.slopes[simulated.slopes >= observed.coef])/N # how many simulated co-efficients are greater than the observed values.

p.value.2  # Remember that this p-value is not really zero....

##

# Now we have done this for about the simplest null model possible. But the power of this approach is to use for far more complex "null" models where you are trying to understand the distribution of effects under various models. We can also use this approach for other parameters of interest, like the coefficient of determination (R^2), which is one of the exercises we will do in class for a somewhat more interesting model

summary(model.object)$r.sq
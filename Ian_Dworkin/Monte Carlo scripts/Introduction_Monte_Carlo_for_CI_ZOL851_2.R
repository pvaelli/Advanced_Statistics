
# Ian Dworkin ZOL851: Tutorial on using monte carlo simulations for generating confidence intervals. 
# This tutorial puts the logic together more formally, but for a very simple example.
# Updated October 9th 2012

 # Here is an example of a simple Monte Carlo simulation, using fake data....
x <- -10:10 # generating a co-variate of "known" values.
b0 <- 0
b1 <- 1.5  # slope of the model

# Now we generate our response variable y, based on x.
# For the purposes of this tutorial, think of this vector y as the observed values for our response
y <-rnorm(length(x), b0 + b1*x, 10)
 
# We plot the relationship between y and x 
plot(y~x) # Scatterplot of y on x


lm.1 <- lm(y~x) # linear model of y ~ x. Before looking, think what the intercept should be? The slope?
coef(lm.1)
confint(lm.1)

vcov(lm.1) # This is the co-variance between parameter estimates. For the moment it does not matter much. Why?
summary(lm.1)

# And plot the fitted line (estimated from this simple regression)
abline(lm.1, lwd=3)


#####ASYMPTOTIC NORMAL CONFIDENCE BANDS  (FOR COMPARISON)  ##################
# Before we use monte carlo simulation to construct confidence intervals, let us start by using the asymptotic normal CIs (i.e. what you would normally construct based on the SE and the t-distribution).

#We use predict() to construct the asymptotic normal confidence intervals for the observed data. 

# creating a mock data frame to put in values for the predictor values which will be use to "predict values" of the response , including the confidence intervals.
x.frame <- data.frame(x=x) # the name of the predictor here needs to be the same as in the model frame 

# We are generating values of x using the actual x values, but we could also generate this just using a sequence, going from the minimum to maximum values of x as well. Take a look at x.frame to see what I mean. What IS IMPORTANT is that x is sorted!!!!

# generating confidence intervals. Here we are generating confidence intervals for the values of x in "x.frame" based on the model (lm.1). For each value of x it is providing a fitted value (first column) for y (i.e. the expected value, which should be the same as the abline), and the upper and lower confidence intervals. 
conf.int <- predict (lm.1, interval="conf", newdata=x.frame) 

# Note you can also use interval="pred" for CIs for future predicted observations of the response. These will have greater uncertainty (and this wider CIs).

#take a quick look at the predicted values (and their CIs)
head(conf.int)

# Now we plot the lines. We could use lines() to plot each of the three lines (estimates, lower CI, upperCI).
# Instead here we use the matlines (matrix of lines) function to plot all simultaneously.
# In this came for each column of y, it is plotting y vs x.
matlines(x = x.frame$x, y = conf.int, lty=c(1,2,2), lwd=2, col="red") 



#### Function to perform the Monte Carlo Simulation (Not accounting for uncertainty in Residual variation)
# Here is a little function that simulates a linear model with intercept a, slope b and standard deviation rse. It then performs a linear regression, and extracts co-efficients. This actually produces confidence intervals that are a bit too narrow (not appropriately conservative). Why? What is missing?
# a = intercept, b = slope, rse equals residual standard error for model.
# The default is to use the intercept and slope from the above model (lm.1) extracted from the coef(lm.1) and the rse from summary(lm.1)$sigma

SimReg1 <- function(mod.input = lm.1){
	a = coef(mod.input)[1] # intercept from model
	b = coef(mod.input)[2] # slope from model
	rse = summary(mod.input)$sigma  # residual standard "error" from model
	y.sim <- rnorm(n = length(x), mean = a + b*x, sd=rse) # generate simulated responses
	lm.sim <- lm(y.sim ~ x) # fit model of simulated responses onto predictor
	coef(lm.sim)} # extract coefficients (intercept and slope)

# number of simulations to perform
N <- 1000 

# Using replicate to repeat simulation N times.
simulated.coef <- replicate(N, SimReg1())

# For convenience, we transpose the output of simulations to have 2 columns (intercept and slope) and N rows
simulated.coef <- t(simulated.coef)

sd(simulated.coef[,1]) # This should approximate the standard error for the intercept (think about why)
sd(simulated.coef[,2]) # This should approximate the standard error for the slope

# compare to
summary(lm.1)$coef[,1:2]

# Approximate monte carlo CIs for the slope
quantile(simulated.coef[,2], c(0.025, 0.975))

#compare to asymptotic CIs for the slope
confint(lm.1)[2,]

# As with any parameter estimates we should consider if those estimated in this simulation are correlated.
cor(simulated.coef[,1], simulated.coef[,2])
plot(simulated.coef[,1], simulated.coef[,2])  

# no real covariation between intercept and slope in this case
# I did this on purpose. Think about how I did this (and why it is unlikey for real data)


###### Confidence bands: Asymptotic Normal vs. Monte Carlo for this simple model.
plot(y~x) # Scatterplot of y on x

# plotting the lines using the simulated intercepts and slopes.

for (i in 1:95){ # Why am I only using the first 95 simulations for the plot?
	curve (simulated.coef[i,1] + simulated.coef[i,2] * x, add=T, col="grey")
	}

# add the asymptotic normal CIs back on to the plot
matlines(x = x.frame$x, y = conf.int, lty=c(1,2,2),lwd=2, col="red") # Asymptotic Normal


	
##### Correcting our Simulation to allow for uncertainty in RSE
# As alluded to, the monte carlo CIs generated above may be too narrow, as we have not accounted for uncertainty in the RSE.

# Correcting our simulation to correctly propogate error (sample from RSE and then sample values conditional upon that).
# Remember that the RSE we get is just an estimate!!! It represents a RANDOM sample from a distribution as well.
# rse= residual standard error, df = residual degrees of freedom for model.
# We only need slightly modify our original function (df.sim, rse.sim)
# See Gelman and Hill for description of how we are incorporating uncertainty into the RSE

SimReg2 <- function(mod.input = lm.1){
	a = coef(mod.input)[1] 
	b = coef(mod.input)[2] 
	df.sim <- mod.input$df # residual degrees of freedom 
	rse = summary(mod.input)$sigma  # residual standard "error" from model
	rse.sim <- rse*sqrt(df.sim/rchisq(1, df = df.sim)) # incorporates UNCERTAINTY in RSE.
    y.sim <- rnorm(n = length(x), mean = a + b*x, sd=rse.sim) 
	lm.sim <- lm(y.sim ~ x) 
	coef(lm.sim)}

# And run simulation as before	
regression_simulation_adjusted <- replicate(N, SimReg2())
regression_simulation_adjusted <- t(regression_simulation_adjusted)

# Add lines to the previous plot
for (i in 1:95){
	curve (regression_simulation_adjusted[i,1] + regression_simulation_adjusted[i,2] * x, 
	    add=T, col="blue", lwd=0.8)}
################

# Compare the CIs for slope for the two different monte carlo functions with the asymptotic
quantile(regression_simulation_adjusted[,2], c(0.025, 0.975))
quantile(simulated.coef[,2], c(0.025, 0.975))
confint(lm.1)[2,]


### Note: Remember that we have NOT YET accounted for co-variation in parameter estimates. 
# We will cover this in one of the next screencasts

# A couple of alternatives
# Instead of using replicate can also use the for loop and call the function (which ever you prefer)

monte_carlo_estimates <- matrix(NA, ncol=2, nrow = N) # initialize a matrix to store simulated values

for(i in 1:N) {
	monte_carlo_estimates[i,] <- SimReg2()}

hist(monte_carlo_estimates[,2])


#######  sim() function in the arm library does most of this as well!!!!
# Gelman and Hill have also written some functions to do such simulations, but in very different ways. 
# See the notes in the powerpoint and the readings....

require(arm)
model.sim <- sim(lm.1, n.sims=N)
sim.slopes <- model.sim@coef[,2] # extracts the simulated slopes
quantile(sim.slopes,c(0.025,0.975)) # simulated confidence intervals
confint(lm.1)

# We can plot on the simulated lines to get confidence bands
plot(y~x)
for (i in 1:95){
	curve (model.sim@coef[i,1] + model.sim@coef[i,2] * x, add=T, col="red")
	}
	
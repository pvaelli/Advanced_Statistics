#ZOL851 -In Class tutorial.
# updated October 10th 2014

# In the previous tutorials we generated confidence intervals using monte carlo methods, but the method we used did not necessarily allow us to include the correlations among predictors. Below we address that using a couple of different approaches.

# Setting a few options for clarity of printing on screen.
options(digits=3, show.signif.stars=F)
library(effects)

#Read data into R (change for your own computer)
wing_data = read.csv("~/ZOL851_FS2014/Lectures/Module_7 Simulations and Power/Oct11_2012_classExercise/AfroTempSizeData.csv", header=TRUE)

str(wing_data)
wing_data <- na.omit(wing_data)
wing_data$Lat.std <- scale(abs(wing_data$Lat), center=T, scale=T) # also used abs()
wing_data$Elev.std <- scale(wing_data$Elevation, center=T, scale=T)
wing_data$Lon.std <- scale(wing_data$Lon, center=T, scale=T)

wing_model <- lm(CsizeScaled ~ (Elev.std + Lat.std + Temp + Sex)^2 + Lon , data=wing_data)
summary(wing_model)
confint(wing_model)

# We are particularly interested in the interaction between elevation and latitude
plot(effect("Elev.std*Lat.std", wing_model))

# How much correlation among parameter estimates?
cov2cor(vcov(wing_model))

# How about with elevation and latitude?
with(wing_data, cor(Lat.std, Elev.std))
with(wing_data, plot(jitter(Lat.std), jitter(Elev.std)))


# Here is a more general way of doing the simulation using the observed predictors.

# I have re-written this function so it is more general (although I may have sacrificed some clarity).
# I am explicity using the design matrix for the model (model.matrix())
# the line matr%*%coef(model) 
# essentially just does predictor1*coef1 + predictor2*coef2.....
# I then use the design matrix explicitly for the model call as well

SimulationUnderModel <- function(model = wing_model){
	#extract design matrix
	matr <- model.matrix(model)
	rse = summary(model)$sigma 
	df=model$df 
	
	# incorporate uncertainty in RSE
	rse_sim <- rse*sqrt(df/rchisq(1, df = df)) 
	
	# Simulate data (response) conditional on the simulated RSE.
	y_sim <- rnorm(n = nrow(matr), 
	    mean = matr%*%coef(model), sd=rse_sim)
	# 0 + design matrix (since the intercept is already in the design matrix)
	lm_sim <- lm(y_sim ~ 0 + matr) # fit model with simulated response
	coef(lm_sim)}


N = 1000 # how many iterations of simulations
simulated_coefficients <- replicate(N,SimulationUnderModel())
simulated_coefficients <- t(simulated_coefficients)

#confidence intervals (simulated and asymptotic)
apply(simulated_coefficients, MARGIN=2, quantile, probs=c(0.025, 0.975))
confint(wing_model)



##########Explicitly modeling the covariation among predictors while constructing approximate confidence intervals using Monte Carlo Simulation
 
library(MASS)	#load library for function mvrnorm()

# one important thing to remember is that the estimated parameter values for b0 and b1 may covary. We need to use this co-variation, so we will assume a joint bivariate normal distribution. mvnnorm().

n <- nrow(wing_data) #number of observations

#Here we generate a random intercept and slope for the observed model and plot the curve.
# Think about the difference here. We are not simulating the response variable (i.e. simulated data) and refitting, but instead we are simulating from the distribution of slopes and intercepts.


ModelSimParameters <- function(model=wing_model, N=1000) {
  
	# inputting the estimated coefs from the model in mu
	# Sigma is the variance covariance matrix for the parameter estimates.
	# we used estimated var-covar matrix for parameters for sigma

	par_mat <- mvrnorm(n=N, # Number of simulated values
	    mu=c(coef(model)), # vector of coefs
	    Sigma=vcov(model))   # matrix of variances and co-variances
	return(par_mat)}

# We can run all the simulations in one call to mvrnorm (no need to use replicate or a for loop)
simulated_slope_intercept <- ModelSimParameters(model=wing_model, N=N)

#confidence intervals (simulated and asymptotic)
apply(simulated_slope_intercept, MARGIN=2, quantile, probs=c(0.025, 0.975))
confint(wing_model)


### NOTE: We could also use this approach to generate co-variation among predictors and generate simulated values at the observation level (not just the parameter level). This can be VERY useful for teasing apart the effects of correlations among observations.



########Gelman and Hill have their own function sim() to perform all of this, but ours works just as well for basic problems
# However for more complex model (GLim and mixed models) the sim() is quite useful.####
# See chapter 7 of Gelman and Hill to see the approach that they use for the simulation (it is a bit different than ours)

library(arm)

#We use the sim function to generate simulations from the observed data
wing_sim <- sim(wing_model, n.sims=N)
sim_slopes <- dll_sim@coef[,2] # extracts the simulated slopes
apply(wing_sim@coef, MARGIN=2, quantile, probs=c(0.025, 0.975))
confint(wing_model)




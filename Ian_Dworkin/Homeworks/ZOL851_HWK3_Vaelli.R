### Patric Vaelli
## October 28, 2014
# ZOL851 Homework Assignment 3

#Calling libraries
require(arm)
require(car)
require(lattice)
#require(asbio)
require(simpleboot)

###### Reading the data into R ######
setwd("/Users/pvaelli/Desktop")
axodata <- read.csv("all_animals.csv", header=TRUE)
axodata <- na.omit(axodata)
axodata$Axolotl_Number <- factor(axodata$Axolotl_Number)
axodata$NS <- axodata$Nutritional_State
axodata$EOG.norm <- scale(axodata$EOG, center=FALSE, scale=TRUE)

#Releveling
#axodata$Stage <- relevel(axodata$Stage, "Baseline")
axodata$Neuropeptide <- relevel(axodata$Neuropeptide, "Ringers") 
axodata$Odorant <- relevel(axodata$Odorant, "IAA")
axodata$Odorant <- factor(axodata$Odorant, c("IAA", "Food", "Male", "Female"))

###### SUBSETTING THE DATA  ######
#Subset to include Ringers and GnRH as our treatment (not to be confused with Stage = Treatment, which is when the neuropeptide is applied; we will examine the effect of Neuropeptide vs ringers control in the "Treatment Stage" of the experiment)
#axodata <- subset(axodata, Neuropeptide=="GnRH" & Stage=="Treatment" & Odorant!="IAA")


###################################### QUESTION 1 ######################################

######## LINEAR MODELS ########
lm.peptide <- lm(EOG.norm ~ Neuropeptide, data=axodata)
lm.GSI <- lm(EOG.norm ~ GSI, data=axodata)
lm.both <- lm(EOG.norm ~ GSI + Neuropeptide, data=axodata)
lm.full <- lm(EOG.norm ~ GSI + Neuropeptide + GSI*Neuropeptide, data=axodata)



summary(lm.both)
confint(lm.full)
summary(lm.full)


######## MONTE CARLO SIMULATIONS ON FULL MODEL ########
SimulationUnderModel <- function(model = lm.full){
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
	coef(lm.sim)
	return(coef(lm.sim))	
}

#Running the simulation 50x
sim.coef.50 <- t(replicate(50, SimulationUnderModel()))
t(apply(sim.coef.50, MARGIN = 2, quantile, probs=c(0.025, 0.975)))

#Now running 1000x
sim.coef.1000 <- t(replicate(1000, SimulationUnderModel()))
t(apply(sim.coef.1000, MARGIN = 2, quantile, probs=c(0.025, 0.975)))

#Now running 5000x
sim.coef.5000 <- t(replicate(5000, SimulationUnderModel()))
CI.1 <- t(apply(sim.coef.5000, MARGIN = 2, quantile, probs=c(0.025, 0.975)))
print(CI.1)

sim.coef.5000 <- t(replicate(5000, SimulationUnderModel()))
CI.2 <- t(apply(sim.coef.5000, MARGIN = 2, quantile, probs=c(0.025, 0.975)))
print(CI.2)

sim.coef.5000 <- t(replicate(5000, SimulationUnderModel())) 
CI.3 <- t(apply(sim.coef.5000, MARGIN = 2, quantile, probs=c(0.025, 0.975)))
print(CI.3)


######## PERMUTATION TEST ON REDUCED MODEL ##########
#Partial linear model exploring effect of Gonadal Somatic Index(GSI) on EOG response.
lm.GSI <- lm(EOG.norm ~ GSI, data=axodata)
summary(lm.GSI)

lm.GSI.fstat  = summary(lm.GSI)$fstat # just pulls out the F statistic
lm.GSI.fstat   
lm.GSI.fstat1 = lm.GSI.fstat[1]
anova(lm.GSI)

### Permutation Function ###
PermuteFunction <- function(y=axodata$EOG.norm, x=axodata$GSI){
	model.resample = lm(sample(y, replace=F) ~ x) 
	#permutes response, then runs model
	fstats = summary(model.resample)$fstat[1] # place all of the fstats in this vector 
	coef(model.resample)
	return(fstats)}

### Distribution of p values under permutation ###
PermutePvals <- function(){
	GSI.resample = lm(sample(EOG.norm, replace=F) ~ GSI, data=axodata) 
	anova(GSI.resample)[1,5]}
	
#Permutation test run 50x
GSI.permute.50 <- replicate(50, PermuteFunction())
GSI.permute.p.50 <- replicate(50, PermutePvals())

#Permutation run 5000x
GSI.permute.5000 <- replicate(5000, PermuteFunction())
GSI.permute.p.5000 <- replicate(5000, PermutePvals())

#To observe F-stats (we don't really want a table with 5000 F stats do we?!)
t(t(GSI.permute.50[1:20]))
t(t(GSI.permute.5000[1:20]))

#Hist of F-statistics from permutation resampling
par(mfrow=c(2,2))
hist(GSI.permute.50, xlab="F Statistic", freq=F, main="F values from F-distribution N=50 (DF=1,901)")
curve(df(x,1,901), add=T, col="red")
arrows(lm.GSI.fstat1, 0.2, lm.GSI.fstat1, 0, lwd = 2)
hist(GSI.permute.5000, xlab="F Statistic", freq=F, main="F values from F-distribution N=5000 (DF=1,901)")
curve(df(x,1,901), add=T, col="red")
arrows(lm.GSI.fstat1, 0.2, lm.GSI.fstat1, 0, lwd = 2)
hist(GSI.permute.p.50, xlab="P-value", freq=F, main="Distribution of P-values, N=50")
hist(GSI.permute.p.5000, xlab="P-value", freq=F, main="Distribution of P-values, N=50")



hist(f.theor, xlim = c(0,65),ylim = c(0,1), xlab = "F Statistic", freq=F,
main = " F values from F-distribution (1,901)")
curve(df(x, 1,901), add=T)


############ BOOTSTRAP ANALYSIS ##############

summary(lm.full) #the full model we ran earlier

#### Pairs resampling function####
PairedBootstrap <- function(dat=axodata, mod.formula=formula(EOG.norm ~ GSI*Neuropeptide)){
  dat.boot <- dat[sample(x = nrow(dat), size = nrow(dat), replace=T),] # samples along index  
  boot.lm <- lm(mod.formula, data=dat.boot)
  coef(boot.lm)}

# 5000 bootstrap iterations
paired.bootstrap <- t(replicate(5000, PairedBootstrap()))

# standard error of the estimates via bootstrap
apply(paired.bootstrap, MARGIN = 2, sd) 

apply(paired.bootstrap, MARGIN = 2, coef) 

# Percentile CIs (transposed so it is oriented like the output of confint)
t(apply(paired.bootstrap, MARGIN = 2, quantile, probs=c(0.025, 0.975)))

# asymptotic normal CIs
confint(lm.full)

### Function to output multiple histograms ###
par(mfrow=c(3,2))
MultipleHistograms <- function(X=paired.bootstrap){
    for (i in 1:ncol(X)) {
	    hist(X[,i], freq=F,
	        main = colnames(X)[i],
	        xlab = colnames(X)[i])}}

MultipleHistograms()	    


#### Residual resampling function ####

lm.full.resid <- resid(lm.full)
plot(density(lm.full.resid, bw=0.5)) 


BootstrapFromResiduals <- function(mod.object = lm.full, dat = axodata) {
	resids = mod.object$resid # extracts residuals from model
	fittedValues = mod.object$fitted # extracts fitted values
	matr <- model.matrix(mod.object)
	# generating new values for each y[i], by adding  bootstrapped resids to fitted values.
	Y <- fittedValues + sample(resids, length(resids), replace=T) 
	
	# Using model.matrix for the predictors (not pretty, I know)
	model.boot <- lm( Y ~ 0 + matr, data=dat ) # refit model with new Y values
	
	coef(model.boot) # Extract the co-efficients
	}

lm.full.bootstrap <- t(replicate(5000, BootstrapFromResiduals()))
t(apply(lm.full.bootstrap, MARGIN = 2, quantile, probs=c(0.025, 0.975)))
confint(lm.full)


# Plot
par(mfrow=c(3,2))
MultipleHistograms(X=lm.full.bootstrap)

############### Problem 1F: Redundancy ################
#original model
coef(lm.full)
confint(lm.full)

# Monte Carlo Simulation (5000x)



################################### QUESTION 2 #####################################

# Using covariance to make inferences about effect size
lm.GSI <- lm(EOG.norm ~ GSI, data=axodata)

library(MASS)
n <- nrow(axodata)

ModelSimParameters <- function(model=lm.GSI) {
  
	# Here we are inputting the estimated co-efficients from the model in mu, and 
	# Sigma is the variance covariance matrix for the parameter estimates.

	par.mat <- mvrnorm(n=1, # Number of simulated values
	    mu=c( coef(model)[1], coef(model)[2]), # vector of "means", in this case mean and intercept
	    Sigma=vcov(model))   # matrix of variances and co-variances
	return(par.mat)    		
}

simulated_slope_intercept <-replicate(5000, ModelSimParameters())
simulated_slope_intercept <- t(simulated_slope_intercept)

#confidence intervals (simulated and asymptotic)
quantile(simulated_slope_intercept[,1], probs=c(0.025, 0.975))
quantile(simulated_slope_intercept[,2], probs=c(0.025, 0.975))

confint(lm.GSI)[,2]


################################### QUESTION 3 ######################################

lm.GSI <- lm(EOG.norm ~ GSI, data=axodata)
summary(lm.GSI)

###### POWER ANALYSIS based on SMEE Code #####

N=1000  # Number of simulations for inner loop. You generally want this to be >1000. 

p = rep(NA, N) # initializing the vector to store the p values in the inner for loop. 

#Global Parameter values
a=0.5 # intercept
b <- seq(from=0, to=1.0,by=0.1)

sample_size <- seq(from=5,to=50,by=1)  # Incremently increasing sample size from 10 to 100 by 10 observations at a time.

power.size <- numeric(length(sample_size)) # initializing the vector to store the power at each sample size for the outer for loop.

### initialize the matrix to store all of the power estimates
power.b <- matrix(NA,length(sample_size),length(b))

## Now the actual for loop
system.time(
for (k in 1:length(b))  # across the different effect sizes
 {
  
  b_b <- b[k]
  
   for (j in 1:length(sample_size))  # looping through the different sample_sizes

    {
   
      s_s = sample_size[j]
      for (i in 1:N)
      {
       x <- rnorm(s_s, mean=8, sd=2)  # simulate values of predictor
       y_det <- a + b_b*x             # deterministic part of model
       y_sim <- rnorm(s_s, mean=y_det,sd=2)  # Simulate y|x values
       lm1 <- lm(y_sim~x)                    # fit model given simulation 
       p[i] <- coef(summary(lm1))[2,4] # You may want to extract a different p-value from the model.
	  
     }
    
      power.size[j] <- length(p[p<0.05])/N   # How many p-values are less than alpha (0.05)
   }
   
    power.b[,k] <- power.size 
}
)


# Now we graph it.
par(mfrow=c(1,1))

#3D perspective plot
persp(y=b, x=sample_size, z=power.b, col="blue", theta=-65, 
    shade=0.75, ltheta=45, ylab="slope", xlab="Sample Size", main="Power Analysis", 
    lphi=30, zlim=c(0,1.25), ticktype = "detailed")

# contour plot
contour(z=power.b, x=sample_size, y=b, col="blue",  ylab="slope", xlab="Sample Size", main="Power Analysis")

#fancy contour
filled.contour(z=power.b, x=sample_size, y=b, 
    ylim=c(min(b),max(b)), xlim=c(min(sample_size), max(sample_size)), 
    xlab="Sample Size", ylab="slope", main="Power Analysis", color = topo.colors)


# Modifying alphas:
alpha <- seq(0.001,0.1,0.01) #initialize a vector of alpha values for (k in 1:length(b)){ 
b_b <- b[k] for (j in 1:length(sample_size)){ 
##remove s_s 
for (i in 1:N){ x <- rnorm(50, mean=8, sd=2) y_det <- a + b_b*x y_sim <- rnorm(50, mean=y_det,sd=2) lm1 <- lm(y_sim~x) p[i] <- coef(summary(lm1))[2,4] 
} power.size[j] <- length(p[p < alpha[j]])/N ##change to 0.05 
alpha[j] 
} 
power.b[,k] <- power.size } 







############################### EXTRAS ###############################
with(axolotl_data, plot(log(Afterpotential), log(EOG), main="Correlation between EOG and After-potential Responses"))

OutlierFinder <- function(variable){
	for (i in variable)
		if (i>(10*sd(variable)))
			print(i)
}

OutlierFinder(axolotl_data$EOG)


#GSI (gonadal somatic index) is my only continuous explanatory variable, so I ran one model with just this variable. May use it later on (perhaps for residual resampling methods where categorical variables cannot be used)
axodata$GSI.cent <- axodata$GSI - mean(axodata$GSI)
lm.GSI <- lm(EOG ~ GSI.cent, data=axodata)


#I subset the data to hold sex and nutritional state constant (male and well-fed, respectively), and only included the Ringers (control) and Neuropeptide (treatment) that we are interested in comparing.
lm.stage <- lm(EOG ~ Stage, data=axodata)
lm.neuropeptide <- lm(EOG ~ Neuropeptide, data=axodata)
lm.both <- lm(EOG ~ Stage + Neuropeptide, data=axodata)
lm.interaction <- lm(EOG ~ Stage + Neuropeptide + Stage * Neuropeptide, data=axodata)

summary()

LinearModelOutput(lm.stage)
LinearModelOutput(lm.neuropeptide)
LinearModelOutput(lm.both)
LinearModelOutput(lm.interaction)
#coefplot(lm.GnRH, intercept=T, vertical=F, ylim=c(-2, 6))
#coefplot(lm.NPY, intercept=T, vertical=F, ylim=c(-2, 3))

# Small function to combine outputs
LinearModelOutput <- function(linear_model){
	print(confint(linear_model))
	print(coef(linear_model))
	#print(cov2cor(vcov(linear_model)))
}


#### 
SimReg <- function(mod.input = lm.interaction){
	a = coef(mod.input)[1] 
	b = coef(mod.input)[2]
	c = coef(mod.input)[3]
	d = coef(mod.input)[4]
	df.sim <- mod.input$df # residual degrees of freedom 
	rse = summary(mod.input)$sigma  # residual standard "error" from model
	rse.sim <- rse*sqrt(df.sim/rchisq(1, df = df.sim)) # incorporates UNCERTAINTY in RSE.
    y.sim <- rnorm(n = length(axodata$EOG), mean = a + b*axodata$Stage=="Treatment" + c*axodata$Stage=="Wash 1" + d*axodata$Stage=="Wash 2" +, sd=rse.sim) 
	lm.sim <- lm(y.sim ~ x) 
	coef(lm.sim)}
	
###########################

require(arm)
model.sim <- sim(lm.GnRH, n.sims=50)
sim.slopes <- model.sim@coef # extracts the simulated slopes
quantile(sim.slopes,c(0.025,0.975)) # simulated confidence intervals
confint(lm.GnRH)


par(mfrow=c(2,1))
plot(EOG ~ Neuropeptide, data=axolotl_data_GnRH, ylim=c(-4, 8))
plot(simulated.EOG ~ Neuropeptide, data=axolotl_data_GnRH)

lm.GnRH.sim <- lm(simulated.EOG ~ Stage + Neuropeptide + Odorant, data=axolotl_data_GnRH)
LinearModelOutput(lm.GnRH.sim)
LinearModelOutput(lm.GnRH)


#### Simulations under Null ####

SimulationUnderNull <- function(model, covariate){
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


lm.GnRH.null <- replicate(N, SimulationUnderNull(lm.GnRH, axolotl_data_GnRH$Neuropeptide))
sd(apply(X = lm.GnRH.null, MARGIN = 2, FUN = mean))
#### Code stolen by Patric Vaelli from Ian Dworkin for performing bootstrap resampling analyses

require(car)
setwd("/Users/pvaelli/Desktop") 
dll.data = read.csv("dll.csv", header=TRUE)   #data frame 
dll.data = na.omit(dll.data)
dll.data$temp <- factor(dll.data$temp)
dll.data$replicate <- factor(dll.data$replicate)
dll.data$genotype <-relevel(dll.data$genotype, ref="wt")

observ.lm <- lm(SCT ~ genotype*temp, data=dll.data)
summary(observ.lm)
confint(observ.lm)



##################################################################

N = 5000  # of Resampling iterations

# Pairs/ random variable resampling
BootstrapRandomX <- function(dat=dll.data, mod.formula=formula(SCT ~ genotype*temp)){
  dat.boot <- dat[sample(x = nrow(dat), size = nrow(dat), replace=TRUE),] # samples along index  
  boot.lm <- lm(mod.formula, data=dat.boot)
  coef(boot.lm)}
  
  # standard error of the estimates via bootstrap
apply(vector.boot, MARGIN = 2, sd) 

# Percentile CIs (transposed so it is oriented like the output of confint)
t(apply(vector.boot, MARGIN = 2, quantile, probs=c(0.025, 0.975)))

# asymptotic normal CIs
confint(observ.lm

##################################################################

# Write a function to perform the bootstrapping
# To make this easy we are extracting everything from the model call to lm
BootstrapFromResiduals <- function(mod.object = observ.lm, dat = dll.data) {
	
	resids = mod.object$resid # extracts residuals from model
	fittedValues = mod.object$fitted # extracts fitted values
	matr <- model.matrix(mod.object)
	# generating new values for each y[i], by adding  bootstrapped resids to fitted values.
	Y <- fittedValues + sample(resids, length(resids), replace=TRUE) 
	
	# Using model.matrix for the predictors (not pretty, I know)
	model.boot <- lm( Y ~ 0 + matr, data=dat ) # refit model with new Y values
	
	coef(model.boot) # Extract the co-efficients
	}

residual.boot.N <- t(replicate(N, BootstrapFromResiduals()))

pairs(residual.boot.N)

t(apply(residual.boot.N, MARGIN = 2, quantile, probs=c(0.025, 0.975)))
confint(observ.lm)

##################################################################

require(boot)



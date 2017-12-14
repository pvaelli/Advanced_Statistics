

# bias in estimates, variance in estimates, and bias-variance tradeoffs.


 # let's say we have a the following data
 
 
 x1 <- rnorm(30, 0, 1)
 
 y <- rnorm(length(x1), 5 + 2*x1, sd = 1.5)
 
 plot(y~x1)
 
 (summary(lm(y~x1)))
 
 
# Now let us say that we have additional co-variates x2-x5, that we think are important (but are actually unrelated to y)

x2 <- rnorm(30, 0, 1)
x3 <- rnorm(30, 0, 1)
x4 <- rnorm(30, 0, 1)
x5 <- rnorm(30, 0, 1)

# What happens to the estimates for the intercept and x1
 (summary(lm(y~x1))[[4]][1:2,1:2])
 (summary(lm(y~ x1  + x2))[[4]][1:2,1:2])
 (summary(lm(y~ x1 + x2 + x3))[[4]][1:2,1:2])
 (summary(lm(y~ x1 + x2 + x3 + x4))[[4]][1:2,1:2])
 (summary(lm(y~ x1 + x2 + x3 + x4 + x5))[[4]][1:2,1:2])
 
 
 # Let's look at this overall a whole bunch of simulations (only worrying about the SE for the coefficient associated with x1)
 
 # Let's write a general function to do this.
 BiasVarianceSimulator <- function(sample_size = 30, b1 = 2, b2 = 0, b3 =0, b4 = 0, b5=0, Model_SE = 1.5){
 	# the b's determine how much influence each of the covariates have on y. Defaults to only x1 influencing it.
 	# Our potential covariates
 	x1 <- rnorm(sample_size, 0, 1)
    x2 <- rnorm(sample_size, 0, 1)
    x3 <- rnorm(sample_size, 0, 1)
    x4 <- rnorm(sample_size, 0, 1)
    x5 <- rnorm(sample_size, 0, 1)
    
    deterministicPartModel <- 5 + b1*x1 + b2*x2 + b3*x3 + b4*x4 + b5*x5
    y <- rnorm(length(x1), deterministicPartModel, sd = Model_SE)
    
    # We are only extracting the SE associated with the parameter estimate for b1 for this simulation
    SE.1 <- (summary(lm(y~ x1))[[4]][2,2])
    SE.2 <- (summary(lm(y~ x1 + x2 ))[[4]][2,2])
    SE.3 <- (summary(lm(y~ x1 + x2 + x3 ))[[4]][2,2])
    SE.4 <- (summary(lm(y~ x1 + x2 + x3 + x4))[[4]][2,2])
    SE.5 <- (summary(lm(y~ x1 + x2 + x3 + x4 + x5))[[4]][2,2])
    
    return(c(SE.1= SE.1, SE.2  = SE.2, SE.3 = SE.3, SE.4 = SE.4, SE.5 = SE.5))
 }
 
 
BiasVarianceSimOutput <- t(replicate(1000,  BiasVarianceSimulator()))
colMeans(BiasVarianceSimOutput)
apply(BiasVarianceSimOutput, 2, sd)
plot(colMeans(BiasVarianceSimOutput), xlab = "number of parameters fit", ylab = "SE of estimate - mean from simulations")

# So all else being equal. When you estimate additional parameters (that contain no information) it decreases your ability to estimate yuor other parameters with as high a degree of certainty.


# Of course life is often much more complicated than this!!!!



# Let's do the same thing, but with several co-variates influencing y (but x1 still being the most important one)

BiasVarianceSimOutput.2 <- t(replicate(1000,  BiasVarianceSimulator(b2=0.25, b3=0.25)))
colMeans(BiasVarianceSimOutput.2)
apply(BiasVarianceSimOutput.2, 2, sd)
plot(colMeans(BiasVarianceSimOutput.2), xlab = "number of parameters fit", ylab = "SE of estimate - mean from simulations")


# And as we increase the effect sizes of those other covariates
BiasVarianceSimOutput.3 <- t(replicate(1000,  BiasVarianceSimulator(b2=0.75, b3=0.5)))
colMeans(BiasVarianceSimOutput.3)
apply(BiasVarianceSimOutput.3, 2, sd)
plot(colMeans(BiasVarianceSimOutput.3), xlab = "number of parameters fit", ylab = "SE of estimate - mean from simulations")



# Finally let us try with all 5 covariates being of equal importance

BiasVarianceSimOutput.4 <- t(replicate(1000,  BiasVarianceSimulator(b1=1, b2=1,b3=1, b4=1, b5=1)))
colMeans(BiasVarianceSimOutput.4)
apply(BiasVarianceSimOutput.4, 2, sd)
plot(colMeans(BiasVarianceSimOutput.4), xlab = "number of parameters fit", ylab = "SE of estimate - mean from simulations")

# So you can see the influence of additional variables on your ability to estimate the variables can vary alot.

# It would also be worth examining the bias in estimation, or variance accounted for, for these simulations, (and also thinking about including some degree of colinearity), which I will leave for an exercise.


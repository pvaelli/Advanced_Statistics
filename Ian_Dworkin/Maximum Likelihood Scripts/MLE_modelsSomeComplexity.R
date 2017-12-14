

# November 1st 2011


# MLE for more complicated models.


# So we have shown how to do mle for a simple regression model. This is pretty easily extended to more complicated models.

### libraries
require(bbmle)

# Input and check data as normal

setwd("/Users/ian/R/R scripts/Dll data/") 
dll.data = read.csv("dll.txt", header=TRUE)   #data frame input
dll.data$temp <- factor(dll.data$temp)
dll.data$replicate <- factor(dll.data$replicate)
dll.data$genotype <- relevel(dll.data$genotype, "wt") 
dll.data <- na.omit(dll.data)
dll.data$tarsus.centered <- scale(dll.data$tarsus, scale=F)



# Let's first start with a simple ANOVA style model, with one factor (with two levels)

Model_1_genotype = lm(SCT ~ genotype, data=dll.data) # same format as aov
summary(Model_1_genotype)$sigma  # RSE for the model
coef(Model_1_genotype)
confint(Model_1_genotype)
logLik(Model_1_genotype)

# Let's remind ourselves of the structure of this simple design matrix
head(model.matrix(Model_1_genotype))
tail(model.matrix(Model_1_genotype))


################MLE
# so how do we code this for our MLE support function (our negative log likelihood calculator)

### We can use the design matrix for temp
x <- with(data=dll.data, model.matrix(~ genotype))  # The design matrix for our simple model

# Let's remind ourselves about is in this simple matrix...
head(x)
tail(x)

# So the first column is all 1s (representing the intercept...). Remember that by default R uses treatment contrasts, so the intercept is the mean for the base level for this factor (in this case "wt")

# So for this particular model we are only interested in the the second column (but this may not always be the case)

# neg. log lik calc
 glm.mle = function(b0,b1, sigma=0.1) { 
 Y.pred = b0 + b1*x[,2]
 -sum(dnorm(dll.data$SCT, mean = Y.pred, sd = sigma, log = TRUE)) 
 } 

# Always worth checking that our nll calculator works (in as much as there are no errors)

glm.mle(1,1,1)

# Let's pick some starting values
start.b0 <- mean(dll.data$SCT)
start.sigma <- sd(dll.data$SCT)

mle.model <- mle2(glm.mle, start=list(b0=start.b0, b1=0, sigma=start.sigma))
# Some warnings. We could double check by using method = 'L-BFGS-B'

summary(mle.model) # gives same parameters and -2logLik

# Let's set up for profiling...
par(mfrow=c(1,3))
prof <- profile(mle.model)  # Store the profiles
plot(prof,conf = c(99, 95, 90, 80, 50)/100, absVal=T)

confint(prof) # uses the profile to get CI's
vcov(mle.model)
logLik(mle.model)


#### Using this for a likelihood ratio test


## How do we calculate the reduced model? IN this case it it very easy!

glm.mle.reduced = function(b0, sigma=0.1) { 
 Y.pred = b0 
 -sum(dnorm(dll.data$SCT, mean = Y.pred, sd = sigma, log = TRUE)) 
 } 

 mle.model.reduced <- mle(glm.mle.reduced, list(b0=start.b0, sigma=start.sigma), method="BFGS")
summary(mle.model.reduced)

 
Likelihood_ratio <- -2*(logLik(mle.model.reduced) - logLik(mle.model))

pchisq(Likelihood_ratio[1], df = 1,lower.tail=F)

# you can use this approach with the anova function using test="Chisq"
anova(Model_1_genotype, test="Chisq")

# Or with a null
Model_1_null = lm(SCT ~ 1, data=dll.data) # same format as aov
anova(Model_1_null, Model_1_genotype, test="Chisq")
 

# As a group let's write out the support function for a cell means model instead!

cell.means.model <- lm(SCT ~ 0 + genotype, data= dll.data)
head(model.matrix(cell.means.model))
tail(model.matrix(cell.means.model))
summary(cell.means.model)

x.cell <- model.matrix(cell.means.model)

cell.means.MLE <- function(b0,b1, sigma=0.1) { 
 Y.pred = b0*x.cell[,1] + b1*x.cell[,2]
 -sum(dnorm(dll.data$SCT, mean = Y.pred, sd = sigma, log = TRUE)) 
 } 

mle.model.cell.means <- mle2(cell.means.MLE, start=list(b0=10, b1=11, sigma=1.6))
profile(mle.model.cell.means)
plot(profile(mle.model.cell.means))
confint(mle.model.cell.means)
### Ok, we have done this for a very simple model. What should we do for a more complicated model



# let's examine try to do the following model as an MLE problem

linear_model_2 <- lm(SCT ~ tarsus.centered + genotype + temp + genotype:temp, data=dll.data)
coef(linear_model_2)
summary(linear_model_2)
-logLik(linear_model_2)
AIC(linear_model_2)
confint(linear_model_2)
summary(linear_model_2)$sigma

# What are the steps we need to take?
# 1 pull out the design matrix for the model
# 2 - write support function/ NLL calculator function
# 3 decide on some starting values for numerical optimization
# 4 - fit the model!!!! (Decide on our "black box" for optimization method)
# 5 - check the model fit, and whether it converged.. (profile)
# IF the model has not converged... then what?
   #- change starting values
    # try different optimization method ("SANN" or Nelder-Mead)


# step 1 get the design matrix

design_matrix <- model.matrix( ~ tarsus.centered + genotype + temp + genotype:temp, data=dll.data)
head(design_matrix)


# 2 - write support function/ NLL calculator function
glm_MLE_support <- function(b0 = mean(dll.data$SCT), b1=0, b2=0, b3=0,b4=0,s0= 1.5, s1=0.25, X=design_matrix, 
 y = dll.data$SCT, sigma = sd(dll.data$SCT) ){
 deterministic_part <- b0*X[,1] + b1*X[, 2] + b2*X[, 3] + b3*X[,4] + b4*X[,5]
 #sigma <- s0 + s1*X[,2]  # This is the line we added to model the additional components of variance
 -sum(dnorm( y, mean = deterministic_part, sd = sigma, log=T)) 
}


# 3 decide on some starting values for numerical optimization

# let's use the standard error of the mean for sigma
SE <- function(x){ 
	sd(x)/sqrt(length(x))}

#4 try a fit	
glm.MLE.fit <- mle2(glm_MLE_support, 
            start=list(b0 = mean(dll.data$SCT), b1 = 0, b2 = 0, 
            b3 =0, b4=0, sigma = SE(dll.data$SCT)), method="BFGS")
#5 look at fit

summary(glm.MLE.fit)
-logLik(glm.MLE.fit)
prof.model <- profile(glm.MLE.fit)
plot(prof.model)


# Let's try a refit

glm.MLE.fit.try2 <- mle2(glm_MLE_support, 
            start=list(b0 = 10.58, b1 = 22.1, b2 = 0.605, 
            b3 =-0.4, b4=1, sigma = 1.51), method="BFGS")
summary(glm.MLE.fit.try2)


glm.MLE.fit.try3 <- mle2(glm_MLE_support, 
            start=list(b0 = 10.58, b1 = 22.1, b2 = 0.605, 
            b3 =-0.4, b4=1, sigma = 1.51), method="SANN")
summary(glm.MLE.fit.try3)
coef(glm.MLE.fit.try3)
-logLik(glm.MLE.fit.try3)

# We have taken the estimates from SANN optimization, and are going back to derivative based optimization
glm.MLE.fit.try4 <- mle2(glm_MLE_support, 
            start=list(b0 = 10.803, b1 = 18.935, b2 = 1.256, 
            b3 =0.137, b4=-1.29, sigma=1.51), method="BFGS")  # Remember to add back s0 & s1 starting values if you want to play with the variance!
prof.model.4 <- profile(glm.MLE.fit.try4)
plot(prof.model.4)
-logLik(glm.MLE.fit.try4)            
-logLik(linear_model_2)           
coef(linear_model_2)
coef(glm.MLE.fit.try4)   
summary(glm.MLE.fit.try4)         
            
# compare sigma from mle to RSE from linear model
summary(linear_model_2)$sigma            
            
            
            
############### IGNORE THIS
# model_line <- lm(SCT ~ line, data=dll.data) 
# head(model.matrix(model_line))
# summary(model_line)
# anova(model_line) 
# logLik(model_line)



# # We can essentially use the same tools we used above, but with a few programming tricks to make our lives easier!

# x <- with(data = dll.data, model.matrix(~ line))

# # First we need to know how many parameters we need to estimate
# N = ncol(x) # extracts # of columns of x, so we know how many parameters we need to estimate.


# # right now you can copy and paste this from b (excluding the quotes and the final "+" at the end). 
# b <- paste( ' b', 0:(N-1), '*x[,', 1:N, '] ', sep='', collapse='+') # creates a vector for the names of the parameters

# b <- noquote(b)  # removes the quotes

# arguments <- paste('b', 0:(N-1), '=1', sep='', collapse= ', ')
# arguments

 # glm.mle = function(b0=0, b1=0, b2=0, b3=0, b4=0, b5=0, b6=0, b7=0, b8=0, b9=0, b10=0, b11=0, b12=0, b13=0, b14=0, b15=0, b16=0, b17=0, b18=0, b19=0, b20=0, b21=0, b22=0, b23=0, b24=0, b25=0, b26=0 , sigma=0.1) { 
 	 # y.pred <- b0*x[,1] + b1*x[,2] + b2*x[,3] + b3*x[,4] + b4*x[,5] + b5*x[,6] + b6*x[,7] + b7*x[,8] 
 	   # + b8*x[,9] + b9*x[,10] + b10*x[,11] + b11*x[,12] + b12*x[,13] + b13*x[,14] + b14*x[,15] + b15*x[,16] 
 	   # + b16*x[,17] + b17*x[,18] + b18*x[,19] + b19*x[,20] + b20*x[,21] + b21*x[,22] + b22*x[,23] 
 	   # + b23*x[,24] + b24*x[,25] + b25*x[,26] + b26*x[,27]   # remember to copy b in
 	 # -sum(dnorm(dll.data$SCT, mean = y.pred, sd = sigma, log = TRUE)) 
 # } 



# mle.model <- mle2(glm.mle, start= list(b0=mean(dll.data$SCT), b1=0, b2=0, b3=0, b4=0, b5=1, b6=1, b7=1,
                    # b8=1, b9=1, b10=1, b11=1, b12=1, b13=1, b14=1, b15=1, b16=1, b17=1, b18=1, b19=1, 
                    # b20=1, b21=1, b22=1, b23=1, b24=1, b25=1, b26=1, sigma=sd(dll.data$SCT)), 
                    # method="Nelder-Mead") 
# logLik(mle.model)
# summary(mle.model)
# plot(profile(mle.model))



# mle.model <- mle2(glm.mle, start = list(b0 = 10.8, b1 = 1.48, 
    # b2 = -0.453352199132526, b3 = 1.20302968899609, b4 = 0.503395906537996, 
    # b5 = 0.370928543508496, b6 = -0.336213973514549, b7 = -0.906087404259375, 
    # b8 = 1.10485134884686, b9 = 1.24117495111627, b10 = 1.25536701657221, 
    # b11 = 1.41499668128165, b12 = 1.04721934146329, b13 = 1.94656185391255, 
    # b14 = 1.21913668390102, b15 = 0.864066893777478, b16 = 1.18163467037058, 
    # b17 = 1.3582523875637, b18 = 1.630373685413, b19 = 1.62022533860867, 
    # b20 = 1.4546946107973, b21 = 1.52906469599927, b22 = 0.319117887518249, 
    # b23 = 2.01857767704545, b24 = 0.868904844927269, b25 = 1.64516987192319, 
    # b26 = 1.88071076620015, sigma = 1.55530819451678), 
                    # method="Nelder-Mead") 


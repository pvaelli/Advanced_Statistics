#As a motivating example, we measured the wingspan of a number of males and females 
#peregrine falcons.  We are interested in whether 
#the sexes differ in size. For Western Europe, Monneret (2006) gives a range of 
#male wingspan as 70-85 cm and that for females as 95-115 cm. 
#Assuming Normal distributions, this implies means and standard deviations of 
#77.5 and 2.5 cm for males and of 105 and 3 for females. 


#We will first look at a model where we assume equal variances with the sexes.

#Simulate the data
n1 <- 60  				# Number of females
n2 <- 40					# Number of males
mu1 <- 105					# Population mean of females
mu2 <- 77.5					# Population mean of males
sigma <- 2.75				# Average population SD of both

n <- n1+n2					# Total sample size
y1 <- rnorm(n1, mu1, sigma)		# Data for females separately
y2 <- rnorm(n2, mu2, sigma)		# Date for males separately
data1 <- data.frame(y=c(y1, y2), sex=c(rep("f",n1),rep("m",n2)))			# Aggregate both data sets
boxplot(data1$y ~ data1$sex, col = "grey", xlab = "Sex", ylab = "Wingspan (cm)", las = 1)

#This is a "means parameterization" approach to generating the data.
#An "effects parameterization" approach to generating data would look like the following:

n <- n1+n2  				# Total sample size
alpha <- mu1				# Mean for females serves as the intercept
beta <- mu2-mu1				# Beta is the difference male-female
E.y <- alpha + beta*x			# expectation
y.obs <- rnorm(n = n, mean = E.y, sd = sigma)	# Add random variation
x <- rep(c(0,1), c(n1, n2))		# Indicator for male
boxplot(y.obs ~ x, col = "grey", xlab = "Male", ylab = "Wingspan (cm)", las = 1)

#You may want repeatedly generate the data to get an idea for sampling variance

#######################################################

#Analysis in R (MLE)

fit1 <- lm(data1$y ~ data1$sex)  		# Analysis the data with an effects parameterization
fit2 <- lm(data1$y ~ data1$sex-1) 		# Analysis the data with a means parameterization
summary(fit1)
summary(fit2)

#Why is there a difference between fit1 and fit2?

#Take a look at the design matrices for the two models (they are the same):
model.matrix(fit1)
model.matrix(fit2)

#What is the interpretation of the coefficients?
fit1$coefficients
fit2$coefficients

#Is there a difference in wingspan between males and females?
#It's easiest to figure this out using fit1.  Why?
confint(fit1)
#But we can also tell with the means parameterization
confint(fit2)
#How?

#Pull out the residuals
residual = fit1$residuals
predicted= fit1$fitted.values

#Is there a pattern in the residuals?
plot(predicted, residual, main = "Residuals vs. predicted values", 
     las = 1, xlab = "Predicted values", ylab = "Residuals")
abline(h = 0)
#What are we looking for here?

#Another way to do this analysis is with the function t.test
t.test(data1$y ~ data1$sex, var.equal=TRUE)
#What is the default parameterization for t.test?
#Suggestion: Use t.test for the homework


#Save your workspace if you want to come back to your work and not have to rerun everything.
save.image()




#Pretend we take a Swiss survey of the Wallcreeper, a little cliff-inhabiting bird 
#that seems to have declined greatly in Switzerland in recent years. Assume that we 
#had data on the proportion of sample quadrats in which the species was observed 
#in Switzerland for the years 1990-2005 and that we were willing to assume that the 
#random deviations about a linear time trend were normally distributed. 

#This is for illustration only, usually, we would use logistic regression or a 
#site-occupancy model to make inference about such data that have to do with the 
#distribution of a species.

#Simulate the data
n <- 16  				# Number of years
a = 40					# Intercept
b = -1.5					# Slope
sigma2 = 25					# Residual variance



#Run this part a few times.  What do we see about the data?

#######################################################

#Analysis in R (MLE)

linreg.fit=

#Inspect the design matrix of that model and print out results


  

#What are some of the important things we should look at in our linear regression?
#R-squared or adjusted R-squared if sample size is small.  Adjusted by the # of dependent data points and # of variables. 
#When the # of variables is small and the # of dependent data points is very large then adjusted R-squared = R-squared.

#Degrees of freedom - Do we have enough power for our covariates, in this case the slope term?

#The effect size of the parameters (i.e., the coefficients)?  What are the parameter values and Standard errors?  
#Does the slope overlap zero?  What is the p-value of the slope?

#Draw a line with the model fit
abline()

#Another way that we could draw the line would be to pull out the coefficients of the model
#alpha=intercept
#beta=slope
alpha=linreg.fit$coefficients[1]
beta=linreg.fit$coefficients[2]
predicted <- alpha+beta*years
lines(years,predicted, col="red", lwd=1)

#Look at residuals 

#Plot the observed vs the predicted


#Add a 1-1 line for reference
abline(0,1)

#Pull out the residuals


#Is there a pattern in the residuals?
plot(predicted, residual, main = "Residuals vs. predicted values", 
     las = 1, xlab = "Predicted values", ylab = "Residuals")
abline(h = 0)

#Determine whether there is a trend.
#What is the slope coefficient?

#Why is your estimate different from b=-1.5, the value we used to generate the data?

#What is the confidence interval for the parameter?  Standard error?


#Is the wallcreeper declining?

predictions <- predict(linreg.fit, interval="confidence")
matplot(years, predictions, lty=c(1,2,2), type="l", lwd=c(4,2,2), col=c("black", "red", "red"))
points(years,y, pch=15, col="grey")

######

#Does this matter if we transform the data by adding 1989?  Check by using x instead of years as the independent variable.
#What changes?
lm(y ~ I(x))













---
title: "Standard_error"
author: "Patric Vaelli"
date: "September 17, 2015"
output: html_document
---

n <- 500

x <- rnorm(n = n, mean = 10, sd=1)

# histogram
hist(x)

#this is same data as histogram above, but it's more accurate in that it's not clumping data into groups
kde <- density(x)
plot(kde, xlim=c(min(x), max(x)), main="Distribution of sample of 500 observations")

# basic 

mean(x)
sd(x)

# this is the forumula for Standard Error:
#sd divided by the square root of the sample size
sd(x)/sqrt(n)

# Writing a function to simulate a distribution of means
Dist1 <- function(n, mean=10, sd=1){
	x <- rnorm(n = n, mean = mean, sd=sd)
	y <- mean(x)
	return(y)	
}

# testing the function
Dist1(n=500)

# Now generating a sample of means, replicating our function 1000 times
sample_means <- replicate(1000, Dist1(n=500))

# Now we plot our simulated means
plot(density(sample_means), xlim=c(min(x), max(x)), main="sampling distribution of the means")

# Whats the standard deviation of this sample of means? It's the standard error from before!!
# SE is essentially the variance in recalculating the mean over and over again.
sd(sample_means)
sd(x)/sqrt(500)


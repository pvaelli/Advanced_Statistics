---
title: "Standard_error"
author: "Patric Vaelli"
date: "September 17, 2015"
output: html_document
---

new_sample <- rnorm(n=150, mean=72, sd=4)

plot(density(new_sample, bw=1.5), ylim=c(0,0.15), xlim=c(55,90), lwd=2, xlab="new sample")

#adds normal dist curve to above plot
curve(dnorm(x, mean=72, sd=4), 55, 90, add=T, col="red", lty=3, lwd=2)

#vertical lines depicting each standard deviation (1, then 2)
abline(v=mean(new_sample)-sd(new_sample), col="grey", lty=5, lwd=2)
abline(v=mean(new_sample)+sd(new_sample), col="grey", lty=5, lwd=2)
abline(v=mean(new_sample)-2*sd(new_sample), col="grey", lty=6, lwd=2)
abline(v=mean(new_sample)+2*sd(new_sample), col="grey", lty=6, lwd=2)

sample_means <- replicate(1000, mean(rnorm(50, 20, 5)))
plot(density(sample_means))
abline(v=mean(sample_means)-2*sd(sample_means), col="grey", lty=6, lwd=2)
abline(v=mean(sample_means)+2*sd(sample_means), col="grey", lty=6, lwd=2)

# SE function. uses length of X vector to determine sample size "n"
Standard_Error <- function(x){
	return(sd(x)/sqrt(length(x)))
}

# Fuction to generate a mean and standard error
Sampling_Function <- function(n=50, mean=20, sd=5){
	one_sample <- rnorm(n, mean, sd)
	return(c(mean=mean(one_sample), SE=Standard_Error(one_sample)))
}

# replicate sampling function with default inputs
replicate(100, Sampling_Function())

# want to transpose this vector for analysis to have two columns, one with mean and other with SE
samples_for_CI <- t(replicate(100, Sampling_Function()))

# 95% CI
ci_95 <- 1.96 * samples_for_CI[,2]

lower <- samples_for_CI[,1] - ci_95
upper <- samples_for_CI[,1] + ci_95

replicates <- 1:length(samples_for_CI[,1])

plot(samples_for_CI[,1] ~ replicates, ylab="sample values", ylim=c(min(lower), max(upper)))

# adding CI to the means from the previous plot!
for (i in 1:100){
	lines(x=c(replicates[i], replicates[i]), y=c(lower[i], upper[i]), lwd=2.75, col="grey")
}
abline(h=20, col="red", lwd=3, lty=3)
# see!! the CI includes the mean 95% of the time! notice the CI's that do not include the mean

# Despite the fact that you can not interpret CIs in the manner you would like, they do approximate that interpretation pretty well..

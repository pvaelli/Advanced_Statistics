
#This computes a likelihood surface over values sd and mean
setwd("/Users/ian/R/R scripts/Dll data/") 
dll.data = read.csv("dll.txt", header=TRUE)   #data frame input
dll.data$temp <- factor(dll.data$temp)
dll.data$replicate <- factor(dll.data$replicate)
dll.data$genotype <- relevel(dll.data$genotype, "wt") 
dll.data <- na.omit(dll.data)
dll.data$tarsus.centered <- scale(dll.data$tarsus, scale=F) 
  
range(dll.data$SCT) # 6  to 22
sd(dll.data$SCT)    # 1.63
#Parameter range for computing likelihood.
mean.x <- seq(10, 12, by=0.01)
sd.x <- seq(1.3, 1.8, by=0.0025)


lik.mean <- numeric(length(mean.x)) # initializing the vector to store the log likelihood for the inner for() loop.

lik.sd <- matrix(NA,length(mean.x),length(sd.x)) # initializing the matrix to store the log likelihoods

for (k in 1:length(sd.x))
 {
  
  b_b <- sd.x[k]
  
   for (j in 1:length(mean.x))

    {
   
      s_s = mean.x[j]
      
      sample.prob2 <- dnorm(x=dll.data$SCT, mean=s_s, sd=b_b, log=T) 
      lik.mean[j] <- -sum(sample.prob2) 
 
     }
   
    lik.sd[,k] <- lik.mean 
   }

###Graphics

par(mfrow=c(1,2))

#3D
persp(x=mean.x,y=sd.x, z=lik.sd, 
    col="pink", theta=0, phi=0, shade=1, lphi=40, ltheta=135,
    xlab="mean", ylab="sd", zlab="-Log Likelihood",
    ticktype="detailed")
 
 
 #Simple contour plot
#contour plot
#contour(x=mean.x,y=sd.x,z=lik.sd, xlab="Mean",ylab="SD", levels=c(25,24,23,22,21,20,19,18,17))


#Contour plot With color
contour(x=mean.x,y=sd.x,z=lik.sd, xlab="Mean",ylab="SD", 
    col=heat.colors(17), levels=seq(3600, 4400, by=50), lwd=2) 
    # change the levels based on the likelihood, or just remove them entirely.



filled.contour(x=mean.x,y=sd.x,z=lik.sd, 
    xlab="Mean",ylab="SD", levels=pretty(c(3650,4400),25))

#Using Trellis style graphics 
require(lattice)
wireframe(lik.sd, xlab="mean", ylab="sd", zlab="Log Likelihood")


levelplot(lik.sd)
contourplot(lik.sd)


# Work in groups to import all of the femur data. Apply this approach to subsets of the dataset for each genotype, and compare the results for each. 



# Using expand.grid & mapply to generate the likelihood surface
grid.lik <- expand.grid(mean.x, sd.x) # generates all possible combinations of the two vectors
dim(grid.lik) # 7381 rows, by 2 columns 
head(grid.lik, n=90) # This shows a small amount of the data


lik.sd.mean <- function(mean, sd, data=femur.sample){
	sample.prob2 <- dnorm(x=data, mean=mean, sd=sd, log=T) 
      lik.mean <- sum(sample.prob2) }

likelihood.values <- mapply(lik.sd.mean, grid.lik[,1], grid.lik[,2]) # We are inputting values from each row of all of the possible combinations of the sd and mean from the grid.lik object.
# This vector (likelihood.values) contains all of the values for the likelihoods we have calculated, but AS A VECTOR (not a matrix like before)
   
grid.lik.2 <- cbind(grid.lik, likelihood.values) 
  
grid.lik.2[which.max(grid.lik.2[,3]),] # One quick way of getting out the approximate maximum likelihood

# The only problem with this approach is that you need to make the matrix of the likelihoods, which are otherwise just stored as a single vector
grid.matrix <- matrix(data=likelihood.values, nrow=length(mean.x), ncol=length(sd.x))

persp(x=mean.x,y=sd.x,z=grid.matrix, col="blue", theta=45, phi=15,shade=0.5, lphi=40,xlab="mean", ylab="sd", zlab="Log Likelihood",ticktype="detailed")

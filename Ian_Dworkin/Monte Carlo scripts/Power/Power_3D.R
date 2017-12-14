# Power Analysis in 3D: Ian Dworkin, last modified October 12th 2014.

# One of the additional advantages to this Monte Carlo Approach is that we can modify several variables, and examine how they influence the power. 

# Plotting multiple series of power simulations on a single graph.
# Power analysis in the linear regression framework

# Below I perform this using a set of nested for loops.  I also have an alternative script (Power3D_ExpandGrid.R), which utilizes a more `R` approach to this using some vectorized functions (expand.grid and mapply). Both are reasonably similar speed wise (although the for loop is better for really large simulations)

N=200  # Number of simulations for inner loop. You generally want this to be >1000. 

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
par(mfrow=c(1,2))

#3D perspective plot
persp(y=b, x=sample_size, z=power.b, col="blue", theta=-65, 
    shade=0.75, ltheta=45, ylab="slope", xlab="Sample Size", 
    lphi=30, zlim=c(0,1.25), ticktype = "detailed")

# contour plot
contour(z=power.b, x=sample_size, y=b, col="blue",  ylab="slope", xlab="Sample Size")

#fancy contour
filled.contour(z=power.b, x=sample_size, y=b, 
    ylim=c(min(b),max(b)), xlim=c(min(sample_size), max(sample_size)), 
    xlab="Sample Size", ylab="slope", color = topo.colors)



# Scree plots using lattice
### This requires a little bit of set-up.
require(lattice)

power.b.t  <- t(power.b)  # Transpose the matrix to make it easier to work with.
dim(power.b.t)


# from a numeric vector of our slopes to factor for the legend (required for xyplot)
lev.pow <- as.factor(b) 


###  We need to write out a long formula for each component of this model for the graphics to grab each row.

# generate the left side of the formula
variable.matrix <- paste("power.b.t[", 1:nrow(power.b.t), ",]", sep="", collapse=" + ") 

# This would have been the ugly thing we had to write out!
variable.matrix 

# Now we Generate the whole formula
formula.1 <- as.formula(paste(variable.matrix, "sample_size", sep="~")) 
# combines variable.matrix with the right side of the formula.

formula.1

# xyplot is in the lattice library.
xyplot(formula.1, type="b", ylab="power", key=simpleKey(levels(lev.pow), space="right", title="slope")) 

### If you want to output
# write.csv(power.b.t,file="power_matrix.csv")
# q(save="no")

# Power Simulation Example using expand.grid and mapply instead of a for loop
# Ian Dworkin October 12 2014.
# This performs the same simulation as the "Pwer_3D.r" script, but using a more `R` like way, without the for loop, but using several vectorized functions.


# Not much evidence for a speed gain as compared to a for loop (may be a touch slower), but the differences seem negligible up to N=200.

#Global Parameter values
a = 0.5 # intercept
b <- seq(from=0, to=1.0,by=0.1)
std_dev = 2
sample_size <- seq(from=5,to=50,by=1)


# use expand grid to get all combinations of b and sample_size
b_N <- expand.grid(b, sample_size)  # may be worth looking at b_N to see what is being stored.
dim(b_N)
colnames(b_N) <- c("b", "sample_size") 

# Here is the function to generate the simulation and fit the model given the simulated data.
SimulatePower <- function(sample_size, b_b, a, std_dev){
	x <- rnorm(sample_size, mean=0, sd=1)
	y_det <- a + b_b*x
    y_sim <- rnorm(sample_size, mean=y_det, sd=std_dev)
    lm1 <- lm(y_sim~x)
    pval <- coef(summary(lm1))[2,4]
 }

# We Can use this for one sample_size and slope
check_it_works <- replicate(1000, SimulatePower(sample_size=15, b_b=0.4, a=0, std_dev=3))
hist(check_it_works, freq=T)

# The basic approach works like this. This goes through all combinations of sample_size and b (in b_N) and runs the SimulationPower().
p_values <- mapply(SimulatePower, 
    sample_size  = b_N$sample_size, b_b = b_N$b, 
    MoreArgs=list(a=0, std_dev=1)) 

system.time(
# And if we want to repeat this, we can do it easily with replicate    
rep_p <- replicate(100, mapply(SimulatePower, 
    sample_size  = b_N$sample_size, b_b = b_N$b, 
    MoreArgs=list(a=0, std_dev=1)) ) 
)

# Each row represents a distinct combination of sample size and slope. Each column an iteration of that simulation
dim(rep_p)

# Now we can compute the power. We use the apply like funcion to get determine the fraction of p-values less than 0.05
power_lev <- apply(rep_p, MARGIN=1, 
    function(x) length(x[x<=0.05])/length(x)) # how many p-values are less than 0.05


# The only problem with this approach is that you need to make the matrix of p-values, which are otherwise just stored as a single vector
grid_matrix <- matrix(data=power_lev, nrow=length(b), ncol=length(sample_size))

persp(x=b ,y=sample_size, z=grid_matrix, col="blue", 
    theta=-10, shade=0.75, phi=15, d=2, r=0.1,
    xlab="slope", ylab="sample size", zlab="power", 
    ticktype="detailed")

contour(z=grid_matrix, y=sample_size, x=b, col="blue",  xlab="slope", ylab="Sample Size")

filled.contour(z=grid_matrix, y=sample_size, x=b, 
    xlim=c(min(b),max(b)), ylim=c(min(sample_size), max(sample_size)), 
    ylab="Sample Size", xlab="slope", color = topo.colors)

# See the other script for additional plotting examples.  
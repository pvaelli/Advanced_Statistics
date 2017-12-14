# As an exercise write a for loop that will run this same model 200 times and plot the lines from each of the simulations (you do not need to plot the simulated data each time). you can use the type="n" option in the plot to suppress the points. I will post the "answer" Monday night!

# Here is one way to do it!

a=5
b=0.7
x <- seq(2,20) # our predictor variable
y_fixed <- a + b*x  # the deterministic part of our model
par(mfrow=c(1,1))
plot(y_fixed~x, type="n", ylab= "y", xlab = "x")  # The deterministic relationship.


N <- 200
for (i in 1:N){
	y.sim.loop <- rnorm(length(x), mean=y_fixed, sd=2.5)
	y.sim.loop.lm <- lm(y.sim.loop ~ x)
	abline(y.sim.loop.lm, lwd=1, col="grey", lty=2) 
}
abline(a=5, b=0.7, lwd=2) # fit based on "known" values (in black)


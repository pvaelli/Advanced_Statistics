# Overfitting!!


x <- 1:10
y <- rnorm(length(x), x, 3)
plot(y~x)

lm.1 <- lm(y~x)
abline(lm.1)

lm.2 <- lm(y ~ x + I(x^2))
curve(coef(lm.2)[1] + x*coef(lm.2)[2] + (x^2)*coef(lm.2)[3], col="red", add=T)

lm.5 <- lm(y ~ x + I(x^2) + I(x^3) + I(x^4) + I(x^5)  )

curve(coef(lm.5)[1] + x*coef(lm.5)[2] + (x^2)*coef(lm.5)[3] + (x^3)*coef(lm.5)[4] + (x^4)*coef(lm.5)[5]
      + (x^5)*coef(lm.5)[6], col="blue", add=T)
      
# Perfect fit...but at what cost... Can you generalize from this (will it be useful for prediction)      
lm.10 <- lm(y ~ x + I(x^2) + I(x^3) + I(x^4) + I(x^5) + I(x^6) + I(x^7) + I(x^8) + I(x^9) )

curve(coef(lm.10)[1] + x*coef(lm.10)[2] + (x^2)*coef(lm.10)[3] + (x^3)*coef(lm.10)[4] + (x^4)*coef(lm.10)[5]
      + (x^5)*coef(lm.10)[6] + (x^6)*coef(lm.10)[7] + (x^7)*coef(lm.10)[8] + (x^8)*coef(lm.10)[9]+ (x^9)*coef(lm.10)[10], col="purple", add=T)
      
      
# Let's examine how predictive it may be from another set of data that comes from the same underlying "process" (in otherwords except for random variation, the mechanistic process for the two samples, y & y.2 are the same)
y.2 <- rnorm(length(x), x, 3)
plot(y.2 ~x)  
abline(lm.1)  # Not too bad
curve(coef(lm.2)[1] + x*coef(lm.2)[2] + (x^2)*coef(lm.2)[3], col="red", add=T)
curve(coef(lm.5)[1] + x*coef(lm.5)[2] + (x^2)*coef(lm.5)[3] + (x^3)*coef(lm.5)[4] + (x^4)*coef(lm.5)[5]
      + (x^5)*coef(lm.5)[6], col="blue", add=T)
curve(coef(lm.10)[1] + x*coef(lm.10)[2] + (x^2)*coef(lm.10)[3] + (x^3)*coef(lm.10)[4] + (x^4)*coef(lm.10)[5]
      + (x^5)*coef(lm.10)[6] + (x^6)*coef(lm.10)[7] + (x^7)*coef(lm.10)[8] + (x^8)*coef(lm.10)[9]+ (x^9)*coef(lm.10)[10], col="purple", add=T)
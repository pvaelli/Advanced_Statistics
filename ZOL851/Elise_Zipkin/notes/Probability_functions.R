# wrote a function to compare simulations of the probabilty vs the actual probability

binomial_p_est <- function(n=50, size=10, prob=1/6){
  x <- rbinom(n=n, prob=prob, size=size)
  prob_est <- mean(x)/ size
  prob_diff <- prob - prob_est
  return(prob_diff) 
}

binomial_p_est()
prob_1000 <- replicate(1000, binomial_p_est())

#plots
hist(prob_1000)
plot(density(prob_1000))


x <- rgamma(n=1000, shape=1, scale=1)
hist(x, 20)

x <- rgamma(n=1000, shape=10, scale=10)
hist(x, 20)

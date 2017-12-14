#ONE WAY ANOVA

#Assume five groups with 10 snakes measured in each with SVL averages of 
#50,40,45,55, and 60. This corresponds to a base line population mean of 
#50 and effects of populations 2-5 of -10, -5, 5, and 10. We choose a 
#residual standard deviation of SVL of 3 and assemble everything

#Simulate the data
ngroups <-5 # Number of populations
nsample <-10 # Number of snakes in each
pop.means <-c(50,40,45,55,60) # Population mean SVL
sigma <-3 # Residual sd

  # Total number of datapoints
  # Residuals
  # Indicator for population

  #Create design matrix

 # Inspect the data

#Plot the data



#Fit a fixed effect ANOVA in R


#Look at the design matrix for the ANOVA
#Remember R fits an effects parameterization


#What does the anova function do?  



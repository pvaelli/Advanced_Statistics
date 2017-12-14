#Homework

## Due Thurs 11/19 at 11:59 PM

#Use the basic multi-species occupancy modeling code to create a script file to do the following:

#1. Add in a covariate on species specific occupancy for understory foliage (ufc).  Include a linear and square term.
#Don't forget to standardize ufc.
#The covariate effects should be species specific but included as a random effect linked at the community level.
#Turn in the R code and a print out of the mean estimate of the linear and squared effects (on the logit scale).
#(5 pts)

#2. Plot the relationship between ufc and occupancy probability for the hooded warbler (HOWA) over the range of observed
#values of ufc.  (I suggest doing this using the standardized values- for simplicity.)  Use all values of the
#posterior distribution so that you can look at the complete range of uncertainity in the relationship (a thin black lines
#for each iterations of the MCMC). Then plot the mean estimate in red (thicker line).  Label the axes and give 
#the figure a title. Give a few line of text to interpret the results. (5 pts)

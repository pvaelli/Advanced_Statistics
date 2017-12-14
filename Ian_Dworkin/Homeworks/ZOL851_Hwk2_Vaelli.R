### Written by PATRIC VAELLI
#ZOL851 Assignment 2
#Due Oct 7, 2014

#Note: calls or analyses I attempted but did not work are ##'d out 


#Calling libraries
require(arm)
require(car)
require(lattice)
require(plotrix)

#Reading the data into R
setwd("/Users/pvaelli/Desktop")
axolotl_data <- read.csv("all_animals.csv", header=TRUE)
summary(axolotl_data)
str(axolotl_data)

# Minor edits to the data frame
axolotl_data$EOG.centered <- axolotl_data$EOG - mean(axolotl_data$EOG)
axolotl_data$Stage <- factor(axolotl_data$Stage)
axolotl_data$Axolotl_Number <- factor(axolotl_data$Axolotl_Number)
axolotl_data$Neuropeptide <- relevel(axolotl_data$Neuropeptide, "Ringers") 
axolotl_data$Odorant <- relevel(axolotl_data$Odorant, "IAA")
axolotl_data$Odorant <- factor(axolotl_data$Odorant, c("IAA", "Food", "Male", "Female"))

### Examining variation in EOG response across different explanatory variables

par(mfrow=c(3,2)) #add this to visualize all plots simultaneously. However, legends do not scale down...
#par(mfrow=c(1,1)) #run this to cancel the 6 panel visualization

boxplot(EOG ~ Axolotl_Number, data=axolotl_data,
	xlab="Animal Number",
	ylab="EOG Response",
	main="Variation in EOG Responses Across Experiments (Animals)")
	
boxplot(EOG ~ Neuropeptide, data=axolotl_data,
	xlab="Neuropeptide Modulation",
	ylab="EOG Response",
	main="Variation in EOG Responses by Neuropeptide")
	
boxplot(EOG ~ Odorant, data=axolotl_data,
	xlab="Odorant",
	ylab="EOG Response",
	main="Variation in EOG Responses by Odorant")
	
boxplot(EOG ~ Sex, data=axolotl_data,
	xlab="Sex",
	ylab="EOG Response",
	main="Variation in EOG Responses by Sex")
	
boxplot(EOG ~ Nutritional_State, data=axolotl_data,
	xlab="Nutritional State",
	ylab="EOG Response",
	main="Variation in EOG Responses by Nutritional State")
	
boxplot(EOG ~ GSI, data=axolotl_data,
	xlab="GSI",
	ylab="EOG Response",
	main="Variation in EOG Responses by Gonadal Somatic Index")

## Since we are interested in the effects of neuropeptides on odorant responses, let's look at Odorant by Neuropeptide treatments
par(mfrow=c(1,1)) #run this to cancel the 6 panel visualization
boxplot(EOG ~ Neuropeptide:Odorant, data=axolotl_data,
	xlab="Odorant by Peptide",
	ylab="EOG Response",
	main="Variation in EOG Responses by Odorant")
	
#On centered data	
par(mfrow=c(1,1)) #run this to cancel the 6 panel visualization
boxplot(EOG.centered ~ Neuropeptide:Odorant, data=axolotl_data,
	xlab="Odorant by Peptide",
	ylab="EOG Response",
	main="Variation in EOG Responses by Odorant")

	

### QUESTION 2 -- Starting to build some models..

#lm on neuropeptide
lm.neuropeptide <- lm(EOG ~ Neuropeptide, data=axolotl_data)
display(lm.neuropeptide)
coef(lm.neuropeptide)
confint(lm.neuropeptide)
summary(lm.neuropeptide)
plot(EOG ~ Neuropeptide, data=axolotl_data)
plot(lm.neuropeptide)
#dispersion(lm.neuropeptide)
#error.bars(lm.neuropeptide, add=TRUE)
#plotCI(lm.neuropeptide)
sum(lm.neuropeptide$residuals^2) #to get residual SS


#lm on odorant
lm.odorant <- lm(EOG ~ Odorant, data=axolotl_data)
plot(EOG ~ Odorant, data=axolotl_data)
display(lm.odorant)
confint(lm.odorant)
summary(lm.odorant)
plot(lm.odorant)
sum(lm.odorant$residuals^2)

#lm including both
lm.neuropeptide.2way <- lm(EOG ~ Neuropeptide + Odorant, data=axolotl_data)
display(lm.neuropeptide.2way)
confint(lm.neuropeptide.2way)
summary(lm.neuropeptide.2way)
plot(lm.neuropeptide.2way)
vcov(lm.neuropeptide.2way)
sum(lm.neuropeptide.2way$residuals^2)

#lm including both + interaction
lm.neuropeptide.interaction <- lm(EOG.centered ~ 1 + Neuropeptide + Odorant + Odorant:Neuropeptide, data=axolotl_data)
display(lm.neuropeptide.interaction)
confint(lm.neuropeptide.interaction)
summary(lm.neuropeptide.interaction)
anova(lm.neuropeptide.interaction)
plot(lm.neuropeptide.interaction)
vcov(lm.neuropeptide.interaction)
sum(lm.neuropeptide.interaction$residuals^2)
avPlots(lm.neuropeptide.interaction)

#Anova on Neuropeptide
anova.neuropeptide <- aov(EOG ~ Neuropeptide, data=axolotl_data)
anova(anova.neuropeptide)
TukeyHSD(anova.neuropeptide)
confint(anova.neuropeptide)
print(anova.neuropeptide)

#Anova on Odorant
anova.odorant <- aov(EOG ~ Odorant, data=axolotl_data)
anova(anova.odorant)
TukeyHSD(anova.odorant)
confint(anova.odorant)
print(anova.odorant)

#Anova on Neuropeptide, Odorant, and interaction
anova.neuropeptidexodorant <- aov(EOG ~ Odorant*Neuropeptide, data=axolotl_data)
anova(anova.neuropeptidexodorant)
TukeyHSD(anova.neuropeptidexodorant)
confint(anova.neuropeptidexodorant)
print(anova.neuropeptidexodorant)
vcov(anova.neuropeptidexodorant)

# Type III SS ANOVA 
anova.neuropeptide.interaction <- Anova(lm.neuropeptide.interaction, data=axolotl_data)
summary(anova.neuropeptide.interaction)
#display(anova.neuropeptide.interaction)
print(anova.neuropeptide.interaction)

#Trying a within factor model to look at effects within each animal
nested.neuropeptide <- aov(EOG ~ Neuropeptide+Error(Axolotl_Number/Neuropeptide), data=axolotl_data)
summary(nested.neuropeptide)

nested.odorant <- aov(EOG ~ Odorant+Error(Axolotl_Number/Odorant), data=axolotl_data)
summary(nested.odorant)

nested.neuropeptidexodorant <- aov(EOG ~ Neuropeptide*Odorant+Error(Axolotl_Number/Neuropeptide*Odorant), data=axolotl_data)
summary(nested.neuropeptidexodorant)

### QUESTION 6

nested.all_interactions <- aov(EOG ~ Neuropeptide*Odorant*Nutritional_State*GSI+Error(Axolotl_Number/Neuropeptide*Odorant*Nutritional_State*GSI), data=axolotl_data)
summary(nested.all_interactions)
print(nested.all_interactions)
#confint(nested.all_interactions)
#vcov(nested.all_interactions)
#TukeyHSD(nested.all_interactions, conf.level = 0.95)

nested.typeII <- Anova(nested.all_interactions)

avPlots(nested.all_interactions)





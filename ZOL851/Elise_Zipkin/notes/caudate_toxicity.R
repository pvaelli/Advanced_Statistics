library(ggplot2)
library(car)

setwd("~/Desktop")
ttx_data <- read.csv("caudate_toxicity.csv")
str(ttx_data)

xaxis <- as.vector(ttx_data$Taxa)
par(mfrow = c(1,1), oma = c(1,1,1,1), cex.axis = 1.2)
barplot(ttx_data$High_TTX, las=2, col = "black", pos = c(0,0), lwd = 2, ylab="TTX in micrograms (µg)", main="TTX Toxicity Compared Among Caudates")
text(seq(1, 15, by=1.2), par("usr")[3] - 0.2, labels = xaxis, srt = 45, adj=5, pos = 2, xpd = TRUE)

help(text)
help(par)
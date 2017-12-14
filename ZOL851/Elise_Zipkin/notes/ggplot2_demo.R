install.packages("ggplot2")
library("ggplot2")

diamonds # sample dataset present in ggplot2 library
summary(diamonds)
help(diamonds) # can only do this because it's a sample dataset

ggplot(diamonds, aes(x=carat, y=price)) + geom_point() #(dataset, aesthetics, ) + layers for plot()

ggplot(diamonds, aes(x=carat, y=price, color=clarity)) + geom_point() #(dataset, aesthetics, adding color coding) + layers for plot()

ggplot(diamonds, aes(x=carat, y=price, color=cut)) + geom_point() #(dataset, aesthetics, adding color coding) + layers for plot()

ggplot(diamonds, aes(x=carat, y=price, color=clarity, size=cut)) + geom_point() #(dataset, aesthetics, adding color coding) + layers for plot()

ggplot(diamonds, aes(x=carat, y=price, color=clarity, size=cut)) + geom_point(alpha=0.3) #(dataset, aesthetics, adding color coding) + layers for plot()

ggplot(diamonds, aes(x=carat, y=price)) + geom_point() + geom_smooth() #(dataset, aesthetics, adding color coding) + layers for plot()

ggplot(diamonds, aes(x=carat, y=price)) + geom_point() + geom_smooth(se=FALSE) #(dataset, aesthetics, adding color coding) + layers for plot()

ggplot(diamonds, aes(x=carat, y=price, color=clarity)) + geom_point() + geom_smooth(se=FALSE) #(dataset, aesthetics, adding color coding) + layers for plot()

ggplot(diamonds, aes(x=carat, y=price, color=clarity)) + geom_smooth(se=FALSE) #(dataset, aesthetics, adding color coding) + layers for plot()

ggplot(diamonds, aes(x=carat, y=price, color=clarity)) + geom_point() + geom_smooth(se=FALSE) + ggtitle("Colors so many !!") #(dataset, aesthetics, adding color coding) + layers for plot()

ggplot(diamonds, aes(x=carat, y=price, color=clarity)) + geom_point() + geom_smooth(se=FALSE) + ggtitle("Colors so many !!") + xlab("Weight (carats)") + ylab("Price (dollars)")#(dataset, aesthetics, adding color coding) + layers for plot()

# these lines are getting really long. here's a trick:

plot <- ggplot(diamonds, aes(x=carat, y=price, color=cut))
plot <- plot + geom_point()
plot <- plot + ggtitle("Colors so many !!!") 
plot <- plot + xlab("Weight (carats)") + ylab("Price (dollars)")
plot <- plot + xlim(0,2) + ylim(0, 10000)
plot <- plot + scale_y_log10()
plot

# Histograms

ggplot(diamonds, aes(x=price)) + geom_histogram() # by default each bin is 30

ggplot(diamonds, aes(x=price)) + geom_histogram(binwidth=2000) 

ggplot(diamonds, aes(x=price)) + geom_histogram(binwidth=200) 

ggplot(diamonds, aes(x=price)) + geom_histogram(binwidth=200) +facet_wrap(~clarity) 

ggplot(diamonds, aes(x=price)) + geom_histogram(binwidth=200) +facet_wrap(~clarity, scale="free_y") # can add this last part to change the y scaling. Good for looking at shape of distribution; otherwise all y axes are identical across plots

ggplot(diamonds, aes(x=price, fill=cut)) + geom_histogram() # this could be useful for microbial abundance bar plots!

# Density plots

ggplot(diamonds, aes(x=price)) + geom_density()

ggplot(diamonds, aes(x=price, color=cut)) + geom_density()

ggplot(diamonds, aes(x=price, fill=cut)) + geom_density()

ggplot(diamonds, aes(x=price, fill=cut)) + geom_density(alpha=0.3)



# Boxplots

ggplot(diamonds, aes(x=color, y=price)) + geom_boxplot()

ggplot(diamonds, aes(x=color, y=price)) + geom_boxplot() + scale_y_log10()


# Violin plots

ggplot(diamonds, aes(x=color, y=price)) + geom_violin()

ggplot(diamonds, aes(x=color, y=price)) + geom_violin() + scale_y_log10()


# bar plots

ggplot(diamonds, aes(x=cut)) + geom_bar()

plot <- ggplot(diamonds, aes(x=color, fill=cut)) + geom_bar() # stacked vs adjacent (next code)

ggplot(diamonds, aes(x=color, fill=cut)) + geom_bar(position = "dodge")

x <- 0:5
y <- x * 3000
model <- data.frame(x,y)

ggplot(diamonds, aes(x=carat, y=price)) + geom_point() + geom_line(data=model, x=x, y=y)
ggsave(filename="pretty.pdf", plot) # can do .png or .jpeg etc.




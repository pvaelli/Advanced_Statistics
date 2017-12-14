date <- 8
day_of_week <- "Tuesday"

x <- 1:10
y <- exp(1:10)

plot(y~x)

even <- c(0,2,4,6)
log(even)

log_of_even <- log(even)

add_one <- function(x){
  x+1
}

add_one(9)

divide_by_two <- function(x){
  x/2
}

divide_by_two(100)

add_together <- function(x, y) {
  x+y
}

print("this machine kills puppies")


print_both <- function(first, second){
  print(first); print(second)
}
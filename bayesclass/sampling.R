## Inverse CDF method

## Simulate draws from a Normal(5,1) distribution 

n.sim <- 10000
u <- runif(n.sim, 0, 1) #uniform draws
x.draws <- qnorm(u, mean=5, sd=1) #inverse CDF

## Compare with actual Normal(5,1) density

par(mfrow=c(1,2))
x.axis <- seq(2,8,by=.1)
plot(x=x.axis, y=dnorm(x.axis, mean=5), type="l")
lines(density(x.draws), col=2)


## Simulate draws from the following density

f.x <- function(x) {
  if (x >= 0 && x < 0.25)
    8 * x
  else if (x >= 0.25 && x <= 1)
    8/3 - 8 * x/3
  else 0
  }
  
plot(x = seq(0,1,.01), y = unlist(lapply(seq(0,1,.01),
f.x)), col = "blue", xlab = "x", ylab = "Density", main = "", type = "l")


## inverse CDF 

invcdf.func <- function(u){
  if (u >= 0 && u < .25)
  sqrt(u)/2
  else if(u >= .25 && u <= 1)
  1 - sqrt(3*(1-u))/2
}

## sampling

n.sim <- 10000
u <- runif(n.sim, 0, 1) #uniform draws
x <- unlist(lapply(u, invcdf.func))

plot(x = seq(0,1,.01), y = unlist(lapply(seq(0,1,.01),
triangle.function)), col = "blue", xlab = "x", ylab = "Density", main = "", type = "l")
lines(density(x))
legend(x = "topright", legend = c("Inverse CDF Draws", "Target Density"), lty = c(1,1), col = c(1, "blue"))

## why does it work?

par(mfrow=c(1,2))
x.val <- seq(0,1,.01)
pdf.val <- unlist(lapply(x.val, f.x))
plot(x = x.val, y = pdf.val, type = "l", xlab = "x", ylab = "Density",
     main = "PDF") 
abline(v = .2, col = 2)
abline(v = .4, col = 2)
cdf.func <- function(x) {
  if (x < 0)
    0
  else if (x >= 0 && x < 0.25)
    4 * x^2
  else if (x >= 0.25 && x <= 1)
    (8/3) * x - (4/3) * x^2 - (1/3)
  else 1
}
cdf.val <- unlist(lapply(x.val, cdf.func))
plot(x = x.val, y = cdf.val, type = "l", xlab = "x", ylab =
  "Probability (also u)",
     main = "CDF") 
abline(v = .2, col = 2)
abline(v = .4, col = 2)


############################################

## Rejection Sampling

## Simulate draws from the triangle density
## Use Uniform(0,1)*M as envelope density

g.x <- function(x){
  if(x >= 0 && x <= 1)
    1
  else 0
}
M <- 3

n.sim <- 10000
n.draws <- 0
draws <- c()

while(n.draws < n.sim){
  x.c <- runif(1,0,1)
  accept.prob <- f.x(x.c)/(M*g.x(x.c))
  if(runif(1) <= accept.prob){
    draws <- c(draws, x.c)
    n.draws <- n.draws + 1
  }
}

plot(x = x.val, y = unlist(lapply(x.val, f.x)), xlab = "x", ylab =
  "Density", type = "l", col = "blue")
lines(density(draws))
legend(x = "topright", legend = c("Rejection Sampling Draws", "Target Density"), col = c(1,"blue"), lty = c(1,1))
## Binomial model 

## Visualizing Prior and Posterior

bin.beta.plot <- function(y,N,a,b){
  x.axis <- seq(0,1,.01)
  mle <- y/N
  prior <- dbeta(x.axis, shape1=a, shape2=b)
  posterior <- dbeta(x.axis, y+a, N-y+b)
  plot(x=x.axis, y=posterior, xlab=expression(pi), ylab="", type="l", col="blue", axes=F, ylim=c(min(c(posterior,prior)), max(c(posterior,prior))))
  axis(1)
  box()
  lines(x=x.axis, y=prior, col="red")
  abline(v=mle)
  legend("topright", legend=c("MLE", "Prior", "Posterior"), lty=c(1,1,1), col=c(1,"red","blue"))
}


## Uninformative prior
bin.beta.plot(y=41, N=54, a=1, b=1)

## Informative Beta(1,20) prior
bin.beta.plot(y=41, N=54, a=1, b=20)

## Informative Beta(1,20) prior with more n
bin.beta.plot(y=410, N=540, a=1, b=20)


## Using Posterior as a Prior in Bayesian updating

## suppose we observe new observations of size N1 generated from the true pi

bin.beta.update <- function(y,N,a0,b0,true.pi,N1,n.iter){
  x.axis <- seq(0,1,.01)
  prior <- dbeta(x.axis, shape1=a0, shape2=b0)
  posterior <- dbeta(x.axis, y+a0, N-y+b0)
  post.mean <- (y+a0)/(y+a0+N-y+b0)


  par(mfrow=c(1,2))

  for(i in 1:n.iter){

    plot(x=x.axis, y=prior, xlab=expression(pi), ylab="", type="l", col="red", axes=F, ylim=c(min(c(posterior,prior)), max(c(posterior,prior))))
    axis(1)
    box()
    lines(x=x.axis, y=posterior, col="blue")
    abline(v=true.pi)
    legend("topright", legend=c("Truth", "Prior", "Posterior"), lty=c(1,1,1), col=c(1,"red","blue"))

    plot(x=1:i, y=post.mean, col="blue", ylab=expression(pi), xlab="Iteration", xlim=c(0,n.iter), ylim=c(0,1), type="l")
    abline(h=true.pi)
    legend("topright", legend=c("Truth", "Posterior Mean"), lty=c(1,1), col=c(1,"blue"))
    
    if(i == 1) readline("\nPress <return> to start simulation: ")
    new.y <- rbinom(1, size=N1, prob=true.pi)
    N <- N+N1
    y <- y + new.y
    prior <- posterior
    posterior <- dbeta(x.axis, y+a0, N-y+b0)
    post.mean <- c(post.mean, (y+a0)/(y+a0+N-y+b0))
  }
} 


bin.beta.update(y=41, N=54, a=1, b=1, true.pi=.1, N1=1, n.iter=100)
bin.beta.update(y=41, N=54, a=1, b=1, true.pi=.1, N1=1, n.iter=10)
bin.beta.update(y=41, N=54, a=100, b=1, true.pi=.1, N1=10, n.iter=100)



## Posterior summary

y <- 41
N <- 54
a <- 1
b <- 1

# simulation

n.draws <- 10000
post.draws <- rbeta(n.draws, y+a, N-y+b)

## simulation versus analytical

# mean
analytical.mean <- (y+a)/(a+N+b); analytical.mean
sum(post.draws)/n.draws # MC formula
mean(post.draws)

# var
analytical.var <- ((y+a)*(N-y+b))/((a+N+b)^2 * (a+N+b+1)); analytical.var
sum((post.draws - mean(post.draws))^2)/n.draws # MC formula
var(post.draws)

# sd
analytical.sd <- sqrt(analytical.var); analytical.sd
sd(post.draws)


# central credible interval
cent.cred.95 <- quantile(post.draws, probs = c(.025,.975))

# highest posterior density region
library(hdrcde)
hpd.95 <- hdr(post.draws, prob=95)$hdr
hdr.den(post.draws, prob=95)

funky.dist <- c(rbeta(100,1,10), rbeta(100,10,1)) 
hdr(funky.dist, prob=95)$hdr
hdr.den(funky.dist, prob=95)


# probability statements

# p(.6 > pi > .7)
mean(post.draws > .6 & post.draws < .7)

# p(pi > .75)
mean(post.draws > .75)

# p(lakers > celtics)
lakers <- rbeta(n.draws, 41+a, 54-41+b)
celtics <- rbeta(n.draws, 32+a, 50-32+b)
mean(lakers > celtics)


## Prediction

# prior predictive distribution

y <- 41
N <- 54
a <- 1
b <- 1


prior.draws <- rbeta(n.draws,a,b)
prior.pred.dist <- rbinom(n.draws, size=N, prob=prior.draws)
hist(prior.pred.dist)
summary(prior.pred.dist)


# posterior predictive distribution
post.draws <- rbeta(n.draws, y+a, N-y+b)
post.pred.dist <- rbinom(n.draws, size=N, prob=post.draws)
hist(post.pred.dist)

# posterior prediction of # of wins over next 29 games
post.pred.29 <- rbinom(n.draws, size=29, prob=post.draws)
hist(post.pred.29)



## Non-conjugate binomial model 

bin.nc.plot <- function(y,N){
  x.axis <- seq(0,1,.01)

  prior.func <- function(p){
    if(p >=0 & p <= 0.2) 0
    else if(p > 0.2 & p <= 0.5) 1/3    
    else if(p > 0.5 & p <= 0.9) 2 
    else if(p > 0.9 & p <= 1) 1
    else 0
  }
  

  mle <- y/N
  prior <- sapply(x.axis, FUN = prior.func)
  unnormal.post.ord <- prior * x.axis^y * (1-x.axis)^(N-y)

  ## grid method to normalize posterior
  posterior <- unnormal.post.ord/sum(0.01 * unnormal.post.ord) # .01 is size of each grid

  plot(x=x.axis, y=posterior, xlab=expression(pi), ylab="", type="l", col="blue", axes=F, ylim=c(min(c(posterior,prior)), max(c(posterior,prior))))
  axis(1)
  box()
  lines(x=x.axis, y=prior, col="red")
  abline(v=mle)
  legend("topright", legend=c("MLE", "Prior", "Posterior"), lty=c(1,1,1), col=c(1,"red","blue"))
}


bin.nc.plot(y=41,N=54)
bin.nc.plot(y=20,N=54)
bin.nc.plot(y=200,N=540)
bin.nc.plot(y=5,N=54)

## calculating quantities with grid method

unnormal.post.func <- function(p,y,N){
	if(p >=0 & p <= 0.2) 0 
    else if(p > 0.2 & p <= 0.5) 1/3 * p^y * (1-p)^(N-y)    
    else if(p > 0.5 & p <= 0.9) 2 * p^y * (1-p)^(N-y) 
    else if(p > 0.9 & p <= 1) 1 * p^y * (1-p)^(N-y) 
    else 0
	}

k <- .01
grids <- seq(0,1,by=k); grids
unnormal.post.ord <- sapply(grids, FUN=unnormal.post.func, y=41, N=54)
p.y <- sum(k * unnormal.post.ord)

# posterior mean \int pi p(pi | y) dpi
post.mean <- sum(k * grids * (unnormal.post.ord/p.y)); post.mean

# posterior var \int (pi - post.mean)^2 p(pi | y) dpi
post.var <- sum(k * (grids - post.mean)^2 * unnormal.post.ord/p.y); post.var


## via sampling

post.draws2 <- sample(grids, size=n.draws, prob=unnormal.post.ord, replace=T)
mean(post.draws2)
var(post.draws2)


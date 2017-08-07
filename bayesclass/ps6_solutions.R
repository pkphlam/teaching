## Week 6 Problems Solutions

## 1.

## Linear Regression Gibbs Sampler

library(Zelig)
library(MCMCpack)
data(macro)

gibbs.lm <- function(y,X,m,V,nu,delta,sigma2.start,n.sims=5000,burnin=500){
	
	library(mvtnorm)
	
	n <- nrow(X)
	tX <- t(X)
	V.inv <- solve(V)
	
	sigma2.cur <- sigma2.start
	
	beta.update <- function(sigma2,...){
		s2I.inv <- solve(diag(sigma2, nrow=n))
		tX.s2I.inv <- tX %*% s2I.inv
		V.star <- solve(tX.s2I.inv %*% X + V.inv)
		m.star <- V.star %*% (tX.s2I.inv %*% y + V.inv %*% m)
		rmvnorm(1,mean=m.star,sigma=V.star)
		}
	
	sigma2.update <- function(beta,...){
		xb <- X%*%t(beta)
		shape <- (n+nu)/2
		scale <- (t(y-xb)%*%(y-xb) + delta)/2
		rinvgamma(1,shape,scale)
		}
	
	beta.draws <- matrix(NA, nrow=n.sims+burnin, ncol=ncol(X))
	sigma2.draws <- c()
	
	for(i in 1:(n.sims+burnin)){
		beta.draws[i,] <- beta.cur <- beta.update(sigma2.cur)
		sigma2.draws[i] <- sigma2.cur <- sigma2.update(beta.cur)
		print(i)
		}
	
	beta.res <- beta.draws[(burnin+1):nrow(beta.draws),]
	sigma2.res <- sigma2.draws[(burnin+1):length(sigma2.draws)]
	res <- mcmc(cbind(beta.res, sigma2=sigma2.res))
	return(res)
	}
	

sigma2.start <- 1
y <- macro$unem
X <- cbind(1,macro$gdp,macro$trade)
V <- diag(10000,ncol(X))
m <- rep(0,ncol(X))
nu <- delta <- 1

posterior.lm <- gibbs.lm(y,X,m,V,nu,delta,sigma2.start)
lm.check <- MCMCregress(unem~gdp+trade, data=macro)


## Convergence Diagnostics

plot(posterior.lm)
geweke.diag(posterior.lm)
raftery.diag(posterior.lm)
heidel.diag(posterior.lm)

# gelman-rubin
sigma2.starts <- runif(5,0,100)
posterior.lm1 <- gibbs.lm(y,X,m,V,nu,delta,sigma2.starts[1])
posterior.lm2 <- gibbs.lm(y,X,m,V,nu,delta,sigma2.starts[2])
posterior.lm3 <- gibbs.lm(y,X,m,V,nu,delta,sigma2.starts[3])
posterior.lm4 <- gibbs.lm(y,X,m,V,nu,delta,sigma2.starts[4])
posterior.lm5 <- gibbs.lm(y,X,m,V,nu,delta,sigma2.starts[5])
posterior.lm.list <- mcmc.list(posterior.lm1, posterior.lm2, posterior.lm3, posterior.lm4, posterior.lm5)
gelman.diag(posterior.lm.list)



## 2. 

## Probit Regression with Data Augmentation 

data(turnout)

gibbs.probit <- function(y,X,m,V,beta.start,n.sims=5000,burnin=500){
	
	library(msm)
	beta.cur <- t(beta.start)
	
	tX <- t(X)
	V.inv <- solve(V)
	V.star <- solve(tX %*% X + V.inv)
	
	beta.update <- function(ystar,...){
		m.star <- V.star %*% (tX %*% ystar + V.inv %*% m)
		rmvnorm(1,mean=m.star,sigma=V.star)		
		}
	
	n1 <- sum(y==1)
	n0 <- sum(y==0)
	index.1 <- which(y==1)
	index.0 <- which(y==0)
	
	ystar.update <- function(beta,...){
		ystar <- c()
		xb <- X %*% t(beta)
		ystar[index.0] <- rtnorm(n0,mean=xb[index.0],sd=1,upper=0)
		ystar[index.1] <- rtnorm(n1,mean=xb[index.1],sd=1,lower=0)
		return(ystar)
		}
	
	beta.draws <- matrix(NA, nrow=n.sims+burnin, ncol=ncol(X))
	for(i in 1:(n.sims+burnin)){
		ystar.cur <- ystar.update(beta.cur)
		beta.draws[i,] <- beta.cur <- beta.update(ystar.cur)
		print(i)
		}
	
	res <- mcmc(beta.draws[(burnin+1):nrow(beta.draws),])
	return(res)	
	}

y <- turnout$vote
X <- cbind(1, turnout$age, turnout$income)
V <- diag(10000,ncol(X))
m <- rep(0,ncol(X))
beta.start <- rep(0,ncol(X))

posterior.prob <- gibbs.probit(y,X,m,V,beta.start)
probit.check <- MCMCprobit(vote~age+income, data=turnout)


## Convergence Diagnostics
plot(posterior.prob)
geweke.diag(posterior.prob)
raftery.diag(posterior.prob)
heidel.diag(posterior.prob)

# gelman-rubin

beta.starts <- rmvnorm(5, mean=c(0,0,0), sigma=diag(5,3))
posterior.prob1 <- gibbs.probit(y,X,m,V,beta.starts[1,])
posterior.prob2 <- gibbs.probit(y,X,m,V,beta.starts[2,])
posterior.prob3 <- gibbs.probit(y,X,m,V,beta.starts[3,])
posterior.prob4 <- gibbs.probit(y,X,m,V,beta.starts[4,])
posterior.prob5 <- gibbs.probit(y,X,m,V,beta.starts[5,])
posterior.prob.list <- mcmc.list(posterior.prob1, posterior.prob2, posterior.prob3, posterior.prob4, posterior.prob5)
gelman.diag(posterior.prob.list)

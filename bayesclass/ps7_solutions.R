## Week 7 Solutions

library(MCMCpack)
library(msm)
library(mvtnorm)
data(SupremeCourt)


irt.gibbs <- function(y,T,t,V,v, n.sims=1000, burnin=100){
	
	## y has n rows and m columns
		
	index.0 <- which(y == 0)
	index.1 <- which(y == 1)
	index.na <- which(is.na(y))
	n0 <- length(index.0)
	n1 <- length(index.1)
	n.na <- length(index.na)
	n <- nrow(y)
	m <- ncol(y)
	
	ystar.update <- function(x,alpha.beta,...){
		ystar <- c()
		beta <- replicate(n, alpha.beta[,2])
		alpha <- replicate(n, alpha.beta[,1])
		x <- rep(x,m)
		xb <- x*c(t(beta)) - c(t(alpha))
		ystar[index.0] <- rtnorm(n0, mean=xb[index.0], upper=0)
		ystar[index.1] <- rtnorm(n1, mean=xb[index.1], lower=0)
		ystar[index.na] <- rnorm(n.na, mean=xb[index.na])
		ystar <- matrix(ystar,nrow=n,ncol=m)
		return(ystar)
		}
	
	T.inv <- solve(T)
	T.inv.t <- T.inv %*% t
	alpha.beta.update <- function(x,ystar,...){	
		X <- cbind(-1,x)
		tX <- t(X)
		Sigma <- solve(tX %*% X + T.inv)
		reg.func <- function(ystar,...){
			m <- Sigma %*% (tX %*% ystar + T.inv.t)
			rmvnorm(1, mean=m, sigma=Sigma)
			}
		ab <- t(apply(ystar,MARGIN=2,FUN=reg.func))
		return(ab)
		}
	
	V.inv <- solve(V)
	V.inv.v <- V.inv %*% v
	restrict.low <- which(rownames(y) == "Breyer")
	restrict.hi <- which(rownames(y) == "Scalia")
	restrict <- c(restrict.low,restrict.hi)
	no.restrict <- (1:n)[-restrict]
	x.update <- function(ystar,alpha.beta,...){
		x.res <- c()
		alpha <- alpha.beta[,1]
		w <- ystar + matrix(alpha,nrow=n,ncol=m,byrow=T)
		B <- as.matrix(alpha.beta[,2])
		tB.B <- t(B) %*% B
		Sigma <- solve(tB.B + V.inv)
		reg.func <- function(w,...){
			m <- Sigma %*% (t(B) %*% w + V.inv.v)
			rnorm(1,mean=m, sd=sqrt(Sigma))
			}
		x <- apply(w[no.restrict,],FUN=reg.func,MARGIN=1)
		x.res[no.restrict] <- x
		V.inv <- solve(.001)
		Sigma.res <- solve(tB.B + V.inv) 
		m.low <- Sigma.res %*% (t(B) %*% w[restrict.low,] + V.inv %*% -2)		
		m.high <- Sigma.res %*% (t(B) %*% w[restrict.hi,] + V.inv %*% 2)		
		x.res[c(restrict.low,restrict.hi)] <- rnorm(2,mean=c(m.low,m.high),sd=sqrt(Sigma.res))
		return(x.res)
		}
	
	x.cur <- rep(0,n)
	alpha.beta.cur <- cbind(rep(0,m), rep(0,m))
	
	alpha.draws <- beta.draws <- matrix(NA, nrow=n.sims+burnin, ncol=m)
	x.draws <- matrix(NA, nrow=n.sims+burnin, ncol=n)
	for(i in 1:(n.sims+burnin)){
		ystar.cur <- ystar.update(x.cur,alpha.beta.cur)
		alpha.beta.cur <- alpha.beta.update(x.cur,ystar.cur)
		alpha.draws[i,] <- alpha.beta.cur[,1]
		beta.draws[i,] <- alpha.beta.cur[,2]
		x.draws[i,] <- x.cur <- x.update(ystar.cur,alpha.beta.cur)
		print(i)
		}
	
	res <- list(ideal=mcmc(x.draws[(burnin+1):(n.sims+burnin),]), alphas=mcmc(alpha.draws[(burnin+1):(n.sims+burnin),]), betas=mcmc(beta.draws[(burnin+1):(n.sims+burnin),]))
	return(res)
	}
	
	
data <- t(SupremeCourt)
n <- nrow(data)
m <- ncol(data)

t <- rep(0,2)
T <- diag(100, nrow=2)
v <- 0
V <- 1

posterior <- irt.gibbs(data,T,t,V,v,n.sims=100000, burnin=5000)
check <- MCMCirt1d(t(SupremeCourt), theta.constraints=list(Scalia=2, Breyer=-2), burnin=500,mcmc=10000,store.item=T, store.ability=T)


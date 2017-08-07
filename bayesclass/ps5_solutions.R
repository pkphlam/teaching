## Week 5 Problems Solution

## Logistic Regression with Metropolis Hastings

library(Zelig)
data(turnout)


logit.mh <- function(y,X,beta.start,prior.mean=0,prior.var=10000,jump.var,n.sims=10000,burnin=0){
	
	library(mvtnorm)
	library(coda)

	k <- ncol(X)
	beta.cur <- t(beta.start)
	
	log.post.func <- function(beta,...){
		pi.i <- 1/(1+exp(-X %*% t(beta)))
		log.like <- sum(dbinom(y, size=1, prob=pi.i, log=T))
		log.prior <- dmvnorm(beta, mean=rep(prior.mean,k), sigma=diag(prior.var, k), log=T)
		log.post <- log.like + log.prior
		return(log.post)		
		}
	
	beta.update <- function(beta,...){
 		beta.cand <- rmvnorm(1, mean=beta, sigma=jump.var)		
 		r <- exp(log.post.func(beta.cand) - log.post.func(beta))
 		if(runif(1) <= r)
 			beta.cand
 		else beta
 		}
 		
 		
	beta.update2 <- function(beta,...){
		
		for(i in 1:k){
			beta.cand <- rnorm(1, mean=beta[i], sd=sqrt(jump.var[i]))		
			beta.old <- beta
			beta[i] <- beta.cand
			r <- exp(log.post.func(beta) - log.post.func(beta.old))
			if(runif(1) <= r)
				beta <- beta
			else beta <- beta.old
		}
		return(beta)
		}
		
	
	draws <- matrix(NA, nrow=n.sims+burnin, ncol=k)
	for(i in 1:(n.sims+burnin)){
		draws[i,] <- beta.cur <- beta.update(beta.cur)
		#draws[i,] <- beta.cur <- beta.update2(beta.cur)
		print(i)
		}
	res <- mcmc(draws[(burnin+1):(n.sims+burnin),])
	cat("Acceptance Rate:", 1-rejectionRate(res), "\n")
	return(res)
	}
	
	
y <- turnout$vote
X <- cbind(1,turnout$age, turnout$income)
mle <- glm(vote~age+income, data=turnout, family=binomial)
start.val <- c(0,0,0)

## part 1

jump.var <- vcov(mle)
system.time(posterior <- logit.mh(y=y,X=X,beta.start=start.val, jump.var=jump.var, n.sims=5000, burnin=500))
plot(posterior)

## part 2
## change beta.update to beta.update2

posterior2 <- logit.mh(y=y,X=X,beta.start=start.val, jump.var=c(.01,.0001,.001), n.sims=5000, burnin=500)
plot(posterior2)



check <- MCMClogit(vote~age+income, data=turnout)



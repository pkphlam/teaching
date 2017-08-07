## Metropolis-Hastings algorithm

## Draw from a Gamma distribution

library(coda)
mh.gamma <- function(n.sims, start, burnin, cand.sd, shape, rate){
  theta.cur <- start
  draws <- c()
  
  theta.update <- function(theta.cur, shape, rate){
    theta.can <- rnorm(1, mean=theta.cur, sd=cand.sd)
    accept.prob <- dgamma(theta.can, shape=shape, rate=rate)/dgamma(theta.cur, shape=shape, rate=rate)
    if(runif(1) <= accept.prob)
      theta.can
    else
      theta.cur
  }

  for(i in 1:n.sims){
    draws[i] <- theta.cur <- theta.update(theta.cur, shape = shape, rate = rate)
  }
  
  res <- mcmc(draws[(burnin + 1):n.sims])
  cat("Acceptance Rate:", 1-rejectionRate(res), "\n")
  return(res)
}

mh.draws <- mh.gamma(10000, start = 1, burnin = 1000, cand.sd = 1, shape = 1.7, rate = 4.4)
gamma.draws <- rgamma(10000, shape = 1.7, rate = 4.4)


plot(mh.draws)
autocorr.plot(mh.draws)



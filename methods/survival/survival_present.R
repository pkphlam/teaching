# normal transform plot

expo <- rexp(1000,rate=1)
par(mfrow=c(1,2))
plot(density(expo), xlab = "Y", main = "Y")
plot(density(log(expo)), xlab = "log(Y)", main = "log(Y)")

#censor plot

x <- c(1,2,3,4)
plot(x, type="n", xlim = c(0,4), ylim = c(0,4), xlab = "time", ylab =
"observations", main = "Observation 3 is Censored")
points(2,1,pch=19)
points(3,2,pch=19)
points(3.5,3,pch=21)
points(3,3,pch=19)
segments(x0=0,x1=2,y0=1,y1=1)
segments(x0=0,x1=3,y0=2,y1=2)
segments(x0=0,x1=3.5,y0=3,y1=3)
abline(v=3, lty = "dashed")

library(Zelig)
data(coalition)
library(rgenoud)
X <- coalition[,c("invest", "fract", "polar", "numst2", "crisis")]
X <- as.matrix(cbind(1,X))
T <- coalition$duration
C <- coalition$ciep12
expo.lik <- function(par,t,X,c){
  beta <- par
  lambda <- exp(-(X %*% beta))
  log.lik <- sum(c*(log(lambda) - lambda * t) + (1-c)*(-lambda * t))
  return(log.lik)
}

expo.lik <- function(par,T,X,C){
  beta <- as.matrix(par)
  lambda <- exp(-(X %*% beta))
  log.lik <- sum(C*(-(X%*%beta)) - (lambda * T))
  return(log.lik)
}


my.expo <- optim(par = c(1,1,0,0,0,0), fn = expo.lik, T=T, X=X, C=C,
                 control=list(fnscale=-1, trace=0), hessian=TRUE, method = "BFGS")

my.expo.constant <- optim(par = 1, fn = expo.lik, T=T, X=X[,1], C=C,
                          control=list(fnscale=-1, trace=0),
                          hessian=TRUE, method = "BFGS")

gen.out <- genoud(fn = expo.lik, nvars=6, X=X, T=T, C=C, max=TRUE,
                  print.level = 1, project.path = "heh.txt",
                  default.domains = 1)

library(Zelig)
data(coalition)
z <- zelig(Surv(duration, ciep12)~invest+fract+polar+numst2+crisis,
           data=coalition, model="exp")

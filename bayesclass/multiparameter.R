## Multinomial model

## suicide terrorist target type 1981-2006
## 1 - civilians
## 2 - tourists and foreigners
## 3 - diplomacy
## 4 - governors and politicians
## 5 - police
## 6 - soldiers and military personnel

suicide <- c(361, 36, 29, 113, 272, 330) 
n <- sum(suicide)

## Uninformative Prior
alpha <- c(1,1,1,1,1,1)

## Posterior
library(MCMCpack)
n.sim <- 10000
post.draws <- rdirichlet(n.sim, alpha=suicide+alpha)


## Suppose we are interested in probability of tourist given soft target (civilians+tourists)
## a = tour/(civ+tour)
## We can simulate from Dirichlet posterior

civ.draws <- post.draws[,1]
tour.draws <- post.draws[,2]
a.draws <- tour.draws/(civ.draws+tour.draws)

## Or we can leverage the fact that posterior of a is a Beta distribution
a.draws2 <- rbeta(n.sim, suicide[1]+1, suicide[2]+1)
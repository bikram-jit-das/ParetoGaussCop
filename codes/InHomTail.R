rm(list=ls())
# Required packages
library(mvtnorm)
library(ggplot2)
require(gridExtra)
#library(evir)


source("./codes/HTailGCop.R")

##################################################################################################
# Function to simulate from a bivariate Gaussian copula with correlation rho (margins are naturally uniform)
##################################################################################################
# n: sample size
# rho: correlation parameter in Gaussian Copula

sim.GC <- function(n, rho){
  R <- rbind(c(1,rho),c(rho,1))
  dat <- rmvnorm(n, mean = c(0,0), sigma = R)
  dat[,1] <- qunif(pnorm(dat[,1]))
  dat[,2] <- qunif(pnorm(dat[,2]))
  return(dat)
}

# now transform to margins qmarg1, qmarg2

GC.margin <- function(dat, qmarg1, qmarg2){
  mar<-array(0,dim(dat))
  mar[,1] <- qmarg1((dat[,1]))
  mar[,2] <- qmarg2((dat[,2]))
  return(mar)
}



qn1 <- function(p) qnorm(p)
qn2 <- function(p) qnorm(p)

alphaX=2
alphaY=3
rho =  -0.4



qp1 <- function(p) qpareto(p,a=alphaX)
qp2 <- function(p) qpareto(p,a=alphaY)


qf1 <- function(p) qfrechet(p,a=alphaX)
qf2 <- function(p) qfrechet(p,a=alphaY)


qg1 <- function(p) qgpd(p,xi=1/alphaX)
qg2 <- function(p) qgpd(p,xi=1/alphaY)

b1<-alphaX
b2<-alphaY
g<-(b1+b2-2*rho*sqrt(b1*b2))/(1-rho^2)
r<-rho

# Generate the copula
set.seed(Sys.time())
ns <- 10000  # samples
sim <- sim.GC(n = ns,rho = rho)
#simpar<- GC.margin(sim,qf1,qf2)
simpar<- GC.margin(sim,qp1,qp2)
#simpar<- GC.margin(sim,qg1,qg2)



#pdf("./simplots/Paretoa2a3pl03.pdf",width = 20, height = 7)
#pdf("./simplots/Paretoa2a3pl09.pdf",width = 20, height = 7)
pdf("./simplots/Paretoa2a3neg04.pdf",width = 20, height = 7)
plotdataparsim(simpar,b1,b2,g,r)
dev.off()



#pdf("./simplots/Frecheta2a3pl08.pdf",width = 20, height = 7)
#pdf("./simplots/Frecheta2a3neg03.pdf",width = 20, height = 7)
#pdf("./simplots/Frecheta2a3neg08.pdf",width = 20, height = 7)
#plotdataparsimwithCI(simpar,pct=0.95,b1,b2,g,r)
#dev.off()


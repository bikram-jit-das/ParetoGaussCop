rm(list=ls())
# Required packages
library(mvtnorm)
library(ggplot2)
require(gridExtra)
library(fitdistrplus)
library(evd)
library(igraph)


source("./codes/HTailGCop.R")


x<-read.table("soc-Ep.txt",header=T)
socdata<-graph_from_data_frame(x)
degin<-degree(socdata,mode="in")
degout<-degree(socdata,mode="out")
deg<-cbind(degin,degout)

st1<- "Users trusted (in-degree)"
st2<-  "Users trusted by (out-degree)"


pdf("./socnet/degreesocEP.pdf",width = 15, height = 5.5)
plotdataparwithCI(deg,0.95,lt=5,os=1500,st1,st2)
dev.off()



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

alphaX= 1.75
alphaY= 2
rho =  0.8



qp1 <- function(p) qpareto(p,a=alphaX)
qp2 <- function(p) qpareto(p,a=alphaY)


qf1 <- function(p) qfrechet(p,a=alphaX)
qf2 <- function(p) qfrechet(p,a=alphaY)

b1<-alphaX
b2<-alphaY
g<-(b1+b2-2*rho*sqrt(b1*b2))/(1-rho^2)
r<-rho

# Generate the copula
set.seed(Sys.time())
ns <- length(deg[,1])  # samples
sim <- sim.GC(n = ns,rho = rho)
simpar<- GC.margin(sim,qp1,qp2)
a=2
plot(simpar,xlab=st1,ylab=st2,col="darkblue",cex.lab=a)

pdf("./socnet/Epinion.pdf",width = 20, height = 5.5)
par(mfrow=c(1,4))
a=1.7 # size of labels on the axes

QQplotrvwithCI(deg[,1],kpct=0.025,pct=0.9,xi=1/alphaX, xl=a)
QQplotrvwithCI(deg[,2],kpct=0.025,pct=0.9,xi=1/alphaY, xl=a)
plot(deg,xlab=st1,ylab=st2,cex.lab=a)

plot(simpar,xlab=expression("Variable "~X^1),ylab=expression("Variable "~X^2),col="darkblue",cex.lab=a)
dev.off()

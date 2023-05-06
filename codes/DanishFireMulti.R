rm(list=ls())
# Required packages
library(mvtnorm)
library(ggplot2)
require(gridExtra)
library(fitdistrplus)
library(evd)


setwd("/Users/bikram/Dropbox/Projects/InhomogeneousTails/Rcode/")

source("HTailGCop.R")


data(danishmulti)

x<-danishmulti

st1<-"Building"
st2<-"Content"
st3<-"Profit"

#####################

temp<-x[,c(2,3)]
dbc<-temp[(temp[,1]>0)&(temp[,2]>0),]

temp<-x[,c(3,4)]
dcp<-temp[(temp[,1]>0)&(temp[,2]>0),]

temp<-x[,c(2,4)]
dbp<-temp[(temp[,1]>0)&(temp[,2]>0),]

pdf("./danishplots/Danishbuilcont.pdf",width = 15, height = 5.5)
plotdataparwithCI(dbc,0.95,lt=5,os=1000,st1,st2)
dev.off()

pdf("./danishplots/Danishbuildprofit.pdf",width = 15, height = 5.5)
plotdataparwithCI(dbp,0.95,lt=5,os=1000,st1,st3)
dev.off()

pdf("./danishplots/Danishcontprofit.pdf",width = 15, height = 5.5)
plotdataparwithCI(dcp,0.95,lt=5,os=1000,st2,st3)
dev.off()

##################################################################################################
# Function to simulate from a bivariate Gaussian copula with correlation rho (margins are naturally uniform)
##################################################################################################
# n: sample size
# rho: correlation parameter in Gaussian Copula

# sim.GC <- function(n, rho){
#   R <- rbind(c(1,rho),c(rho,1))
#   dat <- rmvnorm(n, mean = c(0,0), sigma = R)
#   dat[,1] <- qunif(pnorm(dat[,1]))
#   dat[,2] <- qunif(pnorm(dat[,2]))
#   return(dat)
# }
# 
# # now transform to margins qmarg1, qmarg2
# 
# GC.margin <- function(dat, qmarg1, qmarg2){
#   mar<-array(0,dim(dat))
#   mar[,1] <- qmarg1((dat[,1]))
#   mar[,2] <- qmarg2((dat[,2]))
#   return(mar)
# }
# 
# 
# 
# qn1 <- function(p) qnorm(p)
# qn2 <- function(p) qnorm(p)
# 
# alphaX=1.55
# alphaY=1.3
# rho =  0.8
# 
# 
# 
# qp1 <- function(p) qpareto(p,a=alphaX)
# qp2 <- function(p) qpareto(p,a=alphaY)
# 
# 
# qf1 <- function(p) qfrechet(p,a=alphaX)
# qf2 <- function(p) qfrechet(p,a=alphaY)
# 
# b1<-alphaX
# b2<-alphaY
# g<-(b1+b2-2*rho*sqrt(b1*b2))/(1-rho^2)
# r<-rho
# 
# # Generate the copula
# set.seed(Sys.time())
# ns <- 604  # samples
# sim <- sim.GC(n = ns,rho = rho)
# simpar<- GC.margin(sim,qp1,qp2)
# a=2
# plot(simpar,xlab=expression("Variable "~X^1),ylab=expression("Variable "~X^2),col="darkblue",cex.lab=a)
# 
# pdf("./danishplots/Danish1.pdf",width = 20, height = 5.5)
# par(mfrow=c(1,4))
# a=2 # size of labels on the axes
# parfit(temp[,1],k=100,xl=a)
# parfit(temp[,2],k=100,xl=a)
# plot(temp,cex.lab=a)
# 
# plot(simpar,xlab=expression("Variable "~X^1),ylab=expression("Variable "~X^2),col="darkblue",cex.lab=a)
# dev.off()
# 
# 
# pdf("./danishplots/Danish90.pdf",width = 20, height = 5.5)
# par(mfrow=c(1,4))
# a=2 # size of labels on the axes
# QQplotrvwithCI(temp[,1],kpct=0.15,pct=0.9,xi=1/1.55, xl=a)
# QQplotrvwithCI(temp[,2],kpct=0.15,pct=0.9,xi=1/1.3, xl=a)
# plot(temp,cex.lab=a)
# 
# plot(simpar,xlab=expression("Variable "~X^1),ylab=expression("Variable "~X^2),col="darkblue",cex.lab=a)
# dev.off()

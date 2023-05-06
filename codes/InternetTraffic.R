rm(list=ls())
# Required packages
library(mvtnorm)
library(ggplot2)
require(gridExtra)
library(fitdistrplus)
library(evd)
library(igraph)

#setwd("/Users/bikram/Github/ParetoGaussCop/")

source("./codes/HTailGCop.R")

#### Another data Duration Rate
e2eses<-read.table("e2eses.txt",header=F)
m<-dim(e2eses)[1]

data<-cbind(e2eses[,1],e2eses[,3])
str1="Duration"
str2="Rate"


pdf("durationrate.pdf",width = 15, height = 5.5)
plotdataparwithCI(data,0.95,lt=5,os=1500,str1,str2)
dev.off()


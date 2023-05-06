rm(list=ls())
# Required packages
library(mvtnorm)
library(ggplot2)
require(gridExtra)
library(fitdistrplus)
library(evd)
library(igraph)

setwd("/Users/bikram/Dropbox/Projects/InhomogeneousTails/Rcode/")

source("HTailGCop.R")


x<-read.table("cit-Hep.txt",header=T)
socdata<-graph_from_data_frame(x)
degin<-degree(socdata,mode="in")
degout<-degree(socdata,mode="out")
deg<-cbind(degin,degout)

st1<- "Articles cited (in-degree)"
st2<-  "Articles cited by (out-degree)"


pdf("./socnet/degreesCiteHEP.pdf",width = 15, height = 5.5)
plotdataparwithCI(deg,0.95,lt=5,os=1500,st1,st2)
dev.off()



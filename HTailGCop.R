
Hill<-function(x)
{
  ordered <- rev(sort(x))
  ordered <- ordered[ordered[] > 0.]
  n <- length(x)
  loggs <- log(ordered)
  hill <- cumsum(loggs[1.:(n - 1.)])/(1.:(n - 1.)) - loggs[2.:n]
  alpha <- 1./hill
  return(alpha)
}

HillwCI<-function(x,pct=0.95)
{
  ordered <- rev(sort(x))
  ordered <- ordered[ordered[] > 0.]
  n <- length(x)
  loggs <- log(ordered)
  hill <- cumsum(loggs[1.:(n - 1.)])/(1.:(n - 1.)) - loggs[2.:n]
  alpha <- 1./hill
  p<-1-(1-pct)/2
  lci<-alpha-qnorm(p)*alpha/sqrt(1:n)
  uci<-alpha+qnorm(p)*alpha/sqrt(1:n)
  return(list(alpha=alpha,lci=lci,uci=uci))
}



Hillalphaplot<-function(x)
{
  ordered <- rev(sort(x))
  ordered <- ordered[ordered[] > 0.]
  n <- length(x)
  loggs <- log(ordered)
  hill <- cumsum(loggs[1.:(n - 1.)])/(1.:(n - 1.)) - loggs[2.:n]
  hill <- 1./hill
  plot(1.:length(hill), hill, 
       type = "l", xlab = "number of order statistics", 
       ylab = "Hill estimate of alpha", main="Hill plot") 
}



## Pareto random variables: density,  quantiles, cdf, generation.  

dpareto <- function(x, a=2, b=1) a*b^a/x^(a+1)
ppareto <- function(x, a=2, b=1) (x > b)*(1-(b/x)^a)
qpareto <- function(u, a=2, b=1) b/(1-u)^(1/a)
rpareto <- function(n, a=2, b=1) qpareto(runif(n),a,b)

## Frechet random variables:   quantiles, cdf, generation.  

pfrechet <- function(x, a=2) exp(-(x)^(-a))
qfrechet <- function(u, a=2) (-log(u))^(-1/a)
rfrechet <- function(n, a=2) qfrechet(runif(n),a)

## GPD random variables:   quantiles, cdf, generation.  

pgpd <- function(x, xi=1/2, mu=0, sigma=1) (x > mu)*1-(1+(xi-mu)/sigma)^(-1/xi)
qgpd <- function(u, xi=1/2,mu=0, sigma=1) ((1-u)^(-xi)-1)/xi
rgpd <- function(n, xi=1/2,mu=0, sigma=1) qfrechet(runif(n),xi, mu, sigma)

Estrho<-function(simpar, pctci=0.95) ### Needs to be checked and modified but works moderately
{
  xe1<-simpar[,1]
  xe2<-simpar[,2]
  xe12<- apply(simpar,1,min)
  
  n<-length(xe1)

  
  a1hat<-Hill(xe1)
  a2hat<-Hill(xe2)
  a12hat<-Hill(xe12)
  
  rho1<- sqrt(a1hat*a2hat)/a12hat-sqrt(pmax(1+(a1hat*a2hat)/a12hat^2-(a1hat+a2hat)/a12hat,0))
  
  rho<- pmin(pmax(rho1,-1),1)
  
  s1<-a1hat^2
  s2<-a2hat^2
  s12<-a12hat^2
  
  bigsqrt<-sqrt(pmax(a1hat*a2hat+a12hat^2-a12hat*(a1hat+a2hat),0))
  delal1<-(sqrt(a1hat/a2hat) - ((a2hat-a12hat)/bigsqrt))*(1-(bigsqrt==0))/(2*a12hat)
  delal2<-(sqrt(a2hat/a1hat) - ((a1hat-a12hat)/bigsqrt))*(1-(bigsqrt==0))/(2*a12hat)
  delg<-(bigsqrt-sqrt(a1hat*a1hat))/a12hat^2 - ((2*a12hat-a1hat-a2hat)/(2*bigsqrt*a12hat))*(1-(bigsqrt==0))
  
  nu<-sqrt(s1*delal1^2+s2*delal2^2+s12*delg^2)

  p<-1 -(1-pctci)/2
  
  lci<-pmax(rho-qnorm(p)*nu/(sqrt(1:length(rho))),-1)
  uci<-pmin(rho+qnorm(p)*nu/(sqrt(1:length(rho))),1)

  return(list(rho=rho,lci=lci,uci=uci))
}

  



plotdatapar<-function(simpar)
{
xe1<-simpar[,1]
xe2<-simpar[,2]
xe12<- apply(simpar,1,min)

n<-length(xe1)
lim<-20:n/20

a1hat<-Hill(xe1)[lim]
a2hat<-Hill(xe2)[lim]
a12hat<-Hill(xe12)[lim]

#mx<-max(xe1,xe2)
n <- length(c(xe1,xe2))
mx<-sort(c(xe1,xe2),partial=n-5)[n-5]
l<-ceiling(max(a1hat[n/20],a2hat[n/20],a12hat[n/20]))
rho<- sqrt(a1hat*a2hat)/a12hat-sqrt(pmax(1+(a1hat*a2hat)/a12hat^2-(a1hat+a2hat)/a12hat,0))

#rho2<- sqrt(a1hat*a2hat)/a12hat-sqrt(1+(a1hat*a2hat)/a12hat^2-(a1hat+a2hat)/a12hat)

par(mfrow=c(1,3))

plot(simpar,xlim=c(0,mx),ylim=c(0,mx),xlab=expression("Variable "~X^1), ylab=expression("Variable "~X^2),cex=1.5, cex.axis=1.3,cex.lab=1.3)
plot(lim,a1hat, ylim=c(0,l+4),type="l", lty=1,lwd=2, col="red", xlab="order statistics", ylab="Hill estimates",cex=1.5, cex.axis=1.3,cex.lab=1.3)
lines(lim,a2hat, ylim=c(0,l+4),type="l", lty=1,lwd=2, col="blue")
lines(lim,a12hat, ylim=c(0,l+4),type="l", lty=1,lwd=2, col="darkgreen")
legend(30, l+3.5, legend=c(expression("Hill estimate of"~alpha[1]),expression("Hill estimate of"~alpha[2]),expression("Hill estimate of"~gamma)), col=c("red","blue","darkgreen"), lty=1, lwd=2, cex=1.5)
plot(lim,rho, ylim=c(-1,1),type="l", lty=1,lwd=2, col="red", xlab="order statistics", ylab=expression("Estimate of "~rho),cex=1.5, cex.axis=1.3,cex.lab=1.3)
#lines(lim,rho2, ylim=c(-1,1),type="l", lty=1,lwd=2, col="blue")
}




plotdataparonlyrho<-function(simpar)
{
  xe1<-simpar[,1]
  xe2<-simpar[,2]
  xe12<- apply(simpar,1,min)
  
  n<-length(xe1)
  lim<-20:n/20
  
  a1hat<-Hill(xe1)[lim]
  a2hat<-Hill(xe2)[lim]
  a12hat<-Hill(xe12)[lim]
  
  #mx<-max(xe1,xe2)
  n <- length(c(xe1,xe2))
  mx<-sort(c(xe1,xe2),partial=n-5)[n-5]
  l<-ceiling(max(a1hat[n/20],a2hat[n/20],a12hat[n/20]))
  rho<- sqrt(a1hat*a2hat)/a12hat#-sqrt(pmax(1+(a1hat*a2hat)/a12hat^2-(a1hat+a2hat)/a12hat,0))
  
  #rho2<- sqrt(a1hat*a2hat)/a12hat-sqrt(1+(a1hat*a2hat)/a12hat^2-(a1hat+a2hat)/a12hat)

  plot(lim,rho, ylim=c(-1,1),type="l", lty=1,lwd=2, col="red", xlab="order statistics", ylab=expression("Estimate of "~rho))
  #lines(lim,rho2, ylim=c(-1,1),type="l", lty=1,lwd=2, col="blue")
}





plotdataparwithCI<-function(simpar,pct,lt=20,os=1000,str1=expression("Variable "~X^1),str2=expression("Variable "~X^2) )
{
  xe1<-simpar[,1]
  xe2<-simpar[,2]
  xe12<- apply(simpar,1,min)
  
  n<-length(xe1)
  lim<-20:min(n/5,os)
  
  
  a1hat<-HillwCI(xe1,pct)$alpha[lim]
  a1hatlc<-HillwCI(xe1,pct)$lci[lim]
  a1hatuc<-HillwCI(xe1,pct)$uci[lim]
  
  a2hat<-HillwCI(xe2,pct)$alpha[lim]
  a2hatlc<-HillwCI(xe2,pct)$lci[lim]
  a2hatuc<-HillwCI(xe2,pct)$uci[lim]
  
  
  a12hat<-HillwCI(xe12,pct)$alpha[lim]
  a12hatlc<-HillwCI(xe12,pct)$lci[lim]
  a12hatuc<-HillwCI(xe12,pct)$uci[lim]
  
  rho<-Estrho(simpar,pct)$rho[lim]
  rholc<-Estrho(simpar,pct)$lci[lim]
  rhouc<-Estrho(simpar,pct)$uci[lim]
  
  
  n1 <- length(c(xe1,xe2))
  mx<-sort(c(xe1,xe2),partial=n1-5)[n1-5]
  l1<-max(min(ceiling(max(a1hat)),lt),5)
  l2<-max(min(ceiling(max(a2hat)),lt),5)
  l12<-max(min(ceiling(max(a12hat)),lt),5)
  
  
  
  layout(matrix(c(1,2,3,1,4,5), 2, 3, byrow = T))
  plot(simpar,xlim=c(0,mx),ylim=c(0,mx),xlab=str1, ylab=str2,cex=1.8, cex.axis=1.3,cex.lab=1.8)
  plot(lim,a1hat, ylim=c(0,l1),type="l", lty=1,lwd=2, col="red",xlab="Order statistics", ylab="",cex=1.8, cex.axis=1.3,cex.lab=1.8)
  lines(lim,a1hatlc,type="l", lty=1,lwd=2, col="blue")
  lines(lim,a1hatuc,type="l", lty=1,lwd=2, col="blue")
  legend(n/7, l1+3, legend=c(expression("Hill estimate of"~alpha[1]),expression("95% CI")), col=c("red","blue"), lty=1, lwd=2, cex=1.5)
  plot(lim,a12hat, ylim=c(0,l12),type="l", lty=1,lwd=2, col="red", xlab="Order statistics", ylab="",cex=1.8, cex.axis=1.3,cex.lab=1.8)
  lines(lim,a12hatlc,type="l", lty=1,lwd=2, col="blue")
  lines(lim,a12hatuc,type="l", lty=1,lwd=2, col="blue")
  legend(n/7, l12+3, legend=c(expression("Hill estimate of"~gamma),expression("95% CI")), col=c("red","blue"), lty=1, lwd=2, cex=1.5)
  plot(lim,a2hat, ylim=c(0,l2),type="l", lty=1,lwd=2, col="red", xlab="Order statistics", ylab="", cex=1.8, cex.axis=1.3,cex.lab=1.8)
  lines(lim,a2hatlc,type="l", lty=1,lwd=2, col="blue")
  lines(lim,a2hatuc,type="l", lty=1,lwd=2, col="blue")
  legend(n/7, l2+3, legend=c(expression("Hill estimate of"~alpha[2]),expression("95% CI")), col=c("red","blue"), lty=1, lwd=2, cex=1.5)
  plot(lim,rho, ylim=c(-1,1),type="l", lty=1,lwd=2, col="red", xlab="Order statistics", ylab=expression("Estimate of "~rho),cex=1.8, cex.axis=1.3,cex.lab=1.8)
  lines(lim,rholc, ylim=c(-1,1),type="l", lty=1,lwd=2, col="blue")
  lines(lim,rhouc, ylim=c(-1,1),type="l", lty=1,lwd=2, col="blue")
 # legend(n/7, 1, legend=c(expression("Hill estimate of"~rho),expression("95% CI ")), col=c("red","blue"), lty=1, lwd=2, cex=1.5)
}



plotdataparonlyrho<-function(simpar)
{
  xe1<-simpar[,1]
  xe2<-simpar[,2]
  xe12<- apply(simpar,1,min)
  
  n<-length(xe1)
  lim<-20:n/20
  
  a1hat<-Hill(xe1)[lim]
  a2hat<-Hill(xe2)[lim]
  a12hat<-Hill(xe12)[lim]
  
  #mx<-max(xe1,xe2)
  n <- length(c(xe1,xe2))
  mx<-sort(c(xe1,xe2),partial=n-5)[n-5]
  l<-ceiling(max(a1hat[n/20],a2hat[n/20],a12hat[n/20]))
  rho<- sqrt(a1hat*a2hat)/a12hat#-sqrt(pmax(1+(a1hat*a2hat)/a12hat^2-(a1hat+a2hat)/a12hat,0))
  
  #rho2<- sqrt(a1hat*a2hat)/a12hat-sqrt(1+(a1hat*a2hat)/a12hat^2-(a1hat+a2hat)/a12hat)
  
  plot(lim,rho, ylim=c(-1,1),type="l", lty=1,lwd=2, col="red", xlab="order statistics", ylab=expression("Estimate of "~rho))
  #lines(lim,rho2, ylim=c(-1,1),type="l", lty=1,lwd=2, col="blue")
}







plotdataparsim<-function(simpar,b1,b2,g,r)
{
  xe1<-simpar[,1]
  xe2<-simpar[,2]
  xe12<- apply(simpar,1,min)
  
  n<-length(xe1)
  lim<-20:n/5
  
  l<-ceiling(max(b1,b2,g))
  
  a1hat<-Hill(xe1)[lim]
  a2hat<-Hill(xe2)[lim]
  a12hat<-Hill(xe12)[lim]
  
  #mx<-max(xe1,xe2)
  n <- length(c(xe1,xe2))
  mx<-sort(c(xe1,xe2),partial=n-5)[n-5]
  rho<- sqrt(a1hat*a2hat)/a12hat-sqrt(pmax(1+(a1hat*a2hat)/a12hat^2-(a1hat+a2hat)/a12hat,0))
  
  #rho2<- sqrt(a1hat*a2hat)/a12hat-sqrt(1+(a1hat*a2hat)/a12hat^2-(a1hat+a2hat)/a12hat)
  
  par(mfrow=c(1,3))
  
  plot(simpar,xlim=c(0,mx),ylim=c(0,mx),xlab=expression("Variable "~X^1), ylab=expression("Variable "~X^2),cex=1.8, cex.axis=1.3,cex.lab=1.8)
  plot(lim,a1hat, ylim=c(0,l+3),type="l", lty=1,lwd=2, col="red", xlab="order statistics", ylab="Hill estimates",cex=1.8, cex.axis=1.3,cex.lab=1.8)
  lines(lim,a2hat, ylim=c(0,l+3),type="l", lty=1,lwd=2, col="blue")
  lines(lim,a12hat, ylim=c(0,l+3),type="l", lty=1,lwd=2, col="darkgreen")
  abline(h=b1,col="darkred",lty=2,lwd=2)
  abline(h=b2,col="darkblue",lty=2,lwd=2)
  abline(h=g,col="darkgreen",lty=2,lwd=2)
  legend(ceiling(n/20), l+3, legend=c(expression("Hill estimate of"~alpha[1]),expression("Hill estimate of"~alpha[2]),expression("Hill estimate of"~gamma)), col=c("red","blue","green"), lty=1, lwd=2, cex=1.5)
  plot(lim,rho, ylim=c(-1,1),type="l", lty=1,lwd=2, col="red", xlab="order statistics", ylab=expression("Estimate of"~rho),cex=1.8, cex.axis=1.3,cex.lab=1.8)
  abline(h=r,col="darkred",lty= 2,lwd=2)
  #lines(lim,rho2, ylim=c(-1,1),type="l", lty=1,lwd=2, col="blue")
}



plotdataparsimwithCI<-function(simpar,pct,b1,b2,g,r)
{
  xe1<-simpar[,1]
  xe2<-simpar[,2]
  xe12<- apply(simpar,1,min)
  
  n<-length(xe1)
  lim<-20:n/5
  
  
  a1hat<-HillwCI(xe1,pct)$alpha[lim]
  a1hatlc<-HillwCI(xe1,pct)$lci[lim]
  a1hatuc<-HillwCI(xe1,pct)$uci[lim]
  
  a2hat<-HillwCI(xe2,pct)$alpha[lim]
  a2hatlc<-HillwCI(xe2,pct)$lci[lim]
  a2hatuc<-HillwCI(xe2,pct)$uci[lim]
  
  
  a12hat<-HillwCI(xe12,pct)$alpha[lim]
  a12hatlc<-HillwCI(xe12,pct)$lci[lim]
  a12hatuc<-HillwCI(xe12,pct)$uci[lim]
  
  rho<-Estrho(simpar,pct)$rho[lim]
  rholc<-Estrho(simpar,pct)$lci[lim]
  rhouc<-Estrho(simpar,pct)$uci[lim]
  
  
  n1 <- length(c(xe1,xe2))
  mx<-sort(c(xe1,xe2),partial=n1-5)[n1-5]
  l1<-ceiling(max(a1hat))
  l2<-ceiling(max(a2hat))
  l12<-ceiling(max(a12hat))

  
  
  layout(matrix(c(1,2,3,1,4,5), 2, 3, byrow = T))
  plot(simpar,xlim=c(0,mx),ylim=c(0,mx),xlab=expression("Variable "~X^1), ylab=expression("Variable "~X^2),cex=1.5, cex.axis=1.3,cex.lab=1.3)
  plot(lim,a1hat, ylim=c(0,l1+3),type="l", lty=1,lwd=2, col="red",xlab="", ylab="",cex=1.5, cex.axis=1.3,cex.lab=1.3,)
  lines(lim,a1hatlc,type="l", lty=1,lwd=2, col="blue")
  lines(lim,a1hatuc,type="l", lty=1,lwd=2, col="blue")
  abline(h=b1,col="darkred",lty=2,lwd=2)
  legend(n/10, l1+3, legend=c(expression("Hill estimate of"~alpha[1]),expression("95% CI"),expression(alpha[1])), col=c("red","blue","darkred"), lty=c(1,1,2), lwd=2, cex=1.5)
  plot(lim,a12hat, ylim=c(0,l12+3),type="l", lty=1,lwd=2, col="red", xlab="Order statistics", ylab="",cex=1.5, cex.axis=1.3,cex.lab=1.3)
  lines(lim,a12hatlc,type="l", lty=1,lwd=2, col="blue")
  lines(lim,a12hatuc,type="l", lty=1,lwd=2, col="blue")
  abline(h=g,col="darkred",lty=2,lwd=2)
  legend(n/10, l12+3, legend=c(expression("Hill estimate of"~gamma),expression("95% CI"),expression(gamma)), col=c("red","blue","darkred"), lty=c(1,1,2), lwd=2, cex=1.5)
  plot(lim,a2hat, ylim=c(0,l2+3),type="l", lty=1,lwd=2, col="red", xlab="", ylab="", ex=1.5, cex.axis=1.3,cex.lab=1.3)
  lines(lim,a2hatlc,type="l", lty=1,lwd=2, col="blue")
  lines(lim,a2hatuc,type="l", lty=1,lwd=2, col="blue")
  abline(h=b2,col="darkred",lty=2,lwd=2)
  legend(n/10, l2+3, legend=c(expression("Hill estimate of"~alpha[2]),expression("95% CI"),expression(alpha[2])), col=c("red","blue","darkred"), lty=c(1,1,2), lwd=2, cex=1.5)
  plot(lim,rho, ylim=c(-1,1),type="l", lty=1,lwd=2, col="red", xlab="Order statistics", ylab=expression("Estimate of "~rho),cex=1.5, cex.axis=1.3,cex.lab=1.3)
  lines(lim,rholc, ylim=c(-1,1),type="l", lty=1,lwd=2, col="blue")
  lines(lim,rhouc, ylim=c(-1,1),type="l", lty=1,lwd=2, col="blue")
  abline(h=r,col="darkred",lty=2,lwd=2)
  legend(n/10, 1, legend=c(expression("Hill estimate of"~rho),expression("95% CI"),expression(rho)), col=c("red","blue","darkred"), lty=c(1,1,2), lwd=2, cex=1.5)
 
  #lines(lim,a2hat, ylim=c(0,l+4),type="l", lty=1,lwd=2, col="blue")
  #lines(lim,a12hat, ylim=c(0,l+4),type="l", lty=1,lwd=2, col="darkgreen")
  #legend(30, l+3.5, legend=c(expression("Hill estimate of"~alpha[1]),expression("Hill estimate of"~alpha[2]),expression("Hill estimate of"~gamma)), col=c("red","blue","darkgreen"), lty=1, lwd=2, cex=1.5)
  # plot(lim,rho, ylim=c(-1,1),type="l", lty=1,lwd=2, col="red", xlab="order statistics", ylab=expression("Estimate of "~rho),cex=1.5, cex.axis=1.3,cex.lab=1.3)
  #lines(lim,rho2, ylim=c(-1,1),type="l", lty=1,lwd=2, col="blue")
}



parfit<-function(x, k, xl=1)
{
  l <- length(x)
  s <- seq(1./(l + 1.), l/(l + 1.), length = l)
  y <-  - log(1. - s)
  plot(y[(l - k + 1.):l], log(sort(x)[(l - k + 1.):l]), pch =
         21, xlab =  "quantiles of exponential", ylab = "log sorted data",cex.lab=xl) 
  coeffs <- lsfit(y[(l - k + 1.):l], log(sort(x)[(l - k + 1.):l]))$coef
  abline(coeffs[1.], coeffs[2.])
  names(coeffs) <- NULL
  list(logxl = coeffs[1.], alpha = 1./coeffs[2.], bn =
         exp(coeffs[1.] + log(l) * coeffs[	2.]))
}






QQplotrvwithCI<-function(x,kpct=0.25,pct=0.95,xi=1,xl=1)
{
  xs=sort(x,decreasing=TRUE)
  n=length(x)
  k=kpct*n
  del=0.01
  y=xs[1:k]
  z=-log((1:k)/k)
  v=log(y/y[k])
  cilen=qnorm(1-(1-pct)/2)*mean(v)/(del*k)^(1/2)
   #   cilen95=2.2414*mean(v)/(del*k)^(1/2)
   #   cilen99=2.80703*mean(v)/(del*k)^(1/2)
      plot(z,v,t="l",cex.lab=xl,xlab="(scaled) exponential quantiles", ylab="(scaled) log-data quantiles")
      xx=c(-log(((del*k):k)/k),rev(-log(((del*k):k)/k)))
     # yy95=c(v[(del*k):k]+cilen95,rev(v[(del*k):k])-cilen95)
      yy=c(v[(del*k):k]+cilen,rev(v[(del*k):k])-cilen)
    #  yy99=c(v[(del*k):k]+cilen99,rev(v[(del*k):k])-cilen99)
     # polygon(xx,yy99,col="skyblue1",border=NA)
     # polygon(xx,yy95,col="skyblue2",border=NA)
      polygon(xx,yy,col="skyblue3",border=NA)
      oldpar<-par(new=TRUE)
      plot(z,v,t="l",lwd=2,xlab="",ylab="",cex.lab=xl,)
      lines(z,z*xi,col="brown",lty=2,lwd=2)

}



QQPlotSG<-function(x,kpct,delpct,flname,xi)
{
  xs=sort(x,decreasing=TRUE)
  n=length(x)
  pdf(file=paste(flname,".pdf",sep=""),height=4.5,width=8)
  par(mfrow=c(2,3),mar=c(2,2,2,1))
  for (j in 1:2){
    for (i in 1:3){
      k=kpct[i]*n
      del=delpct[j]
      y=xs[1:k]
      z=-log((1:k)/k)
      v=log(y/y[k])
      cilen90=1.95996*mean(v)/(del*k)^(1/2)
      cilen95=2.2414*mean(v)/(del*k)^(1/2)
      cilen99=2.80703*mean(v)/(del*k)^(1/2)
      plot(z,v,t="l")
      xx=c(-log(((del*k):k)/k),rev(-log(((del*k):k)/k)))
      yy95=c(v[(del*k):k]+cilen95,rev(v[(del*k):k])-cilen95)
      yy90=c(v[(del*k):k]+cilen90,rev(v[(del*k):k])-cilen90)
      yy99=c(v[(del*k):k]+cilen99,rev(v[(del*k):k])-cilen99)
      polygon(xx,yy99,col="skyblue1",border=NA)
      polygon(xx,yy95,col="skyblue2",border=NA)
      polygon(xx,yy90,col="skyblue3",border=NA)
      oldpar<-par(new=TRUE)
      plot(z,v,t="l",lwd=2)
      lines(z,z*xi,col="brown",lty=2,lwd=2)
      title(main=paste("k=",k,", delta=",del))
    }
  }
  dev.off()
}

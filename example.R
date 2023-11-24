rm(list=ls())
setwd("D:/RICD")
source("RICD.R")

##MDP(Ro et al.2015)
MDP<-function(X, m1=NULL){
  if(is.null(m1)){
    m1 <- 100
  }
  p=ncol(X) #number of variables in data
  N=nrow(X) #number of observations in data
  h=floor(N/2)+1
  H0=matrix(0,2,m1)
  H0=apply(H0,2,function(X){sort(sample(N,2))}) #initial subset for MDP
  H0=t(H0) #dimension m1*k
  
  H0_LTS=matrix(0,m1,h )  #save the best MDP subset 
  detD_mdp=matrix(0,m1,1) #save the best MDP detD
  trR2_mdp=matrix(0,m1,1) #save the best MDP trR2
  dis_mdp=matrix(0,m1,N)
  
  for(i in 1:m1){
    Y=X[H0[i,],]
    Ybar=apply(Y,2,mean)
    S1=cov(Y)
    D=diag(S1)
    detD=prod(D)
    dis=matrix(0,N,1)
    for (j in 1:N){
      temp2=as.matrix(X[j,]-Ybar)
      dis[j]=t(temp2/D)%*%temp2
    }
    nn=sort(order(dis)[1:h]) 
    crit=100 
    
    k=1
    while(crit!=0 & k<10){
      Y=X[nn,]
      Ybar=apply(Y,2,mean)
      S2=cov(Y)
      D1=diag(S2)
      
      detD=prod(D1)
      dis=matrix(0,N,1)
      for (j in 1:N){
        temp2=as.matrix(X[j,]-Ybar)
        dis[j]=t(temp2/D1)%*%temp2
      }
      nn2=sort(order(dis)[1:h])
      crit=sum(abs(nn2-nn))
      nn=nn2
      k=k+1
    }
    ER=cor(X[nn,])
    trR2=sum(diag(ER%*%ER))-p^2/h
    
    H0_LTS[i,]=nn
    detD_mdp[i]=detD
    trR2_mdp[i]=trR2
    dis_mdp[i,]=dis
  }
  loc_mdp=which.min(detD_mdp)
  list(Hmdp=H0_LTS[loc_mdp,],dis=dis_mdp[loc_mdp,],trR2=trR2_mdp[loc_mdp])
}#####end function MDP

#Refined MDP(RMDP Ro et al.2015)
ReMDP<-function(MDP_result,alpha,X){
  p=ncol(X) #number of variables in data
  N=nrow(X) #number of observations in data
  h=floor(N/2)+1
  delta=alpha/2
  cpn=1+(MDP_result$trR2+p^2/h)/p^1.5	 
  und=p+qnorm(1-delta)*sqrt(2*cpn*MDP_result$trR2)
  scale=median(MDP_result$dis)/p
  od=which(MDP_result$dis>(und*scale))
  nod=length(od)
  
  Y=X[-od,]
  if(nod==0){Y=X}
  newN=N-nod
  Ybar=apply(Y,2,mean)
  
  S1=cov(Y)#
  D=diag(S1)
  detD=prod(D)
  
  dis=matrix(0,N,1)
  for (i in 1:N){
    temp2=as.matrix(X[i,]-Ybar)
    dis[i]=t(temp2/D)%*%temp2
  }
  
  ER=cor(Y)
  trRAW=sum(diag(ER%*%ER))
  cpn=1+trRAW/p^1.5
  trR2=trRAW-p^2/newN
  und=p+qnorm(1-alpha)*sqrt(2*cpn*trR2)
  scale=1+sqrt(2*trR2)/(1-delta)/p/sqrt(2*pi)*exp(-qnorm(1-delta)^2/2)
  cutoff = und*scale
  ou_RMDP <- rep(0, N)
  for(dd in 1:N){
    if(dis[dd] > und*scale ){
      ou_RMDP[dd] <- 1
    }
  }
  list(ou_RMDP=ou_RMDP, dist_RMDP=dis, cutoff=cutoff)
}#end function ReMDP
# Real data example
#---------------------------------------------------------------------------------------------------------
#installing required packages and libraries
library(rrcov)
library(rlist)
library(ggplot2)
library(ggrepel)
library(cowplot)
library(mvoutlier)
# Example 1 
# stockdata <- as.data.frame(read.table('stock2017o.csv', header=TRUE, sep=',', row.names=1, fileEncoding = 'UTF-8'))
stockdata <- as.data.frame(read.table('stock2017c.csv', header=TRUE, sep=',', row.names=1, fileEncoding = 'UTF-8'))
# X <- as.matrix(stockdata)+1
X <- as.matrix(stockdata)
# X <- log(stockdata)[ ,15:39]
X <- log(stockdata)
# X <- (stockdata)[, 43:67]
X <- t(X)
# Display data(Figure S1)
par(mfrow=c(1,2))
matplot(X, type='l', col='gray', lty=1, lwd=1.5, ylim=c(-0.5, 0.5), cex.axis=0.8, cex.lab=1.5, mgp=c(1.8, 0.5, 0),
        xlab='Week', ylab='Return')
matplot(X, type='l', col='gray', lty=1, lwd=1.5, ylim=c(-0.5, 3), cex.axis=0.8, cex.lab=1.5, mgp=c(1.8, 0.5, 0),
        xlab='Week', ylab='Return')

RICD_output <- RICD(X, m=300, alpha=0.01)
n <- nrow(X)
p <- ncol(X)
dist_RICD <- RICD_output$dist_RICD
cutoff_RICD <- RICD_output$cutoff
ou_RICD <- RICD_output$ou_RICD
which(ou_RICD==1)
#[1] 2  8  9 13 22 23
RICD_output <- RICD(X, m=300, alpha=0.05)
ou_RICD <- RICD_output$ou_RICD
which(ou_RICD==1)
#[1]  2  7  8  9 10 13 19 22 23
RICD_output <- RICD(X, m=300, alpha=0.1)
ou_RICD <- RICD_output$ou_RICD
which(ou_RICD==1)
#[1]  2  7  8  9 10 13 19 22 23
MDP_output <- MDP(X, m1=300)
RMDP_output <- ReMDP(MDP_output, 0.01, X)
ou_RMDP <- RMDP_output$ou_RMDP
which(ou_RMDP==1)
# [1]  2  8 13 19 22
MDP_output <- MDP(X, m1=300)
RMDP_output <- ReMDP(MDP_output, 0.05, X)
ou_RMDP <- RMDP_output$ou_RMDP
which(ou_RMDP==1)
# [1]  2  7  8 13 19 21 22
MDP_output <- MDP(X, m1=300)
RMDP_output <- ReMDP(MDP_output, 0.1, X)
ou_RMDP <- RMDP_output$ou_RMDP
which(ou_RMDP==1)
# [1]  2  7  8 13 19 21 22
PCout <- pcout(X, makeplot = FALSE)
#Error in pcout(X, makeplot = FALSE) : 
#More than 50% equal values in one or more variables!

# display data
round(log(dist_RICD), 2)#alpha=0.01
#V15   V16   V17   V18   V19   V20   V21   V22   V23   V24   V25   V26   V27   V28   V29   V30   V31 
#-2.11  0.33 -2.09 -2.29 -2.17 -2.49 -1.79  0.48  0.44 -1.98 -2.12 -2.51  0.48 -2.30 -2.23 -2.14 -2.34 
#V32   V33   V34   V35   V36   V37   V38   V39 
#-2.28 -1.96 -2.28 -2.26  0.67  0.57 -2.20 -2.30
#start:"2017-05-22"
#end:""2017-12-19"
#[1] "2017-06-01" "2017-07-21" "2017-07-31" "2017-09-01" "2017-11-23" "2017-12-01"
log(cutoff_RICD)
#[1] -1.761582
df1 <- as.data.frame(t(rbind(c(1:n), log(dist_RICD))))
df1.highlight <- df1[c(2,  8,  9, 13, 22, 23), ]
df1.text <- df1.highlight
df1.text$label <- c("06-01:0.33", "07-21:0.48", "07-31:0.44", "09-01:0.48", "11-23:0.67", "12-01:0.57")
pic1 <- ggplot(df1, aes(V1, V2)) + geom_point(size = 3, shape = 15)+ 
  geom_point(data = df1.highlight, aes(V1, V2), size = 5, col="black")+
  geom_text_repel(data = df1.text, aes(V1, V2, label = label), color="black", seed = 1, size = 5)+
  geom_hline(yintercept = log(cutoff_RICD), linetype = 'dotted', col = 'darkcyan', size = 1)+
  annotate("text", x = 1, y = log(cutoff_RICD), label = "cutoff:-1.76", vjust = -0.5, size = 5)+
  theme(axis.line = element_line(color="black", size = 0.2), 
        axis.title = element_text(size = 18), plot.title = element_text(size = 18))+
  ylim(-3, 3) + xlim(0, 25) + labs(x = "Index", y = "Modified distance", title = "Method RICD")
pic1
# Test the same data set with RMDP(Ro et al.2015) procedure
MDP_output <- MDP(X, m1=300)
RMDP_output <- ReMDP(MDP_output, 0.01, X)
ou_RMDP <- RMDP_output$ou_RMDP
which(ou_RMDP==1)
# [1]  2  8 13 19 22
dist_RMDP <- RMDP_output$dist_RMDP
t(round(log(dist_RMDP), 2))
#[,1] [,2] [,3] [,4] [,5] [,6] [,7] [,8] [,9] [,10] [,11] [,12] [,13] [,14] [,15] [,16] [,17] [,18]
#[1,] 5.77 6.79 5.92 5.31 5.73 5.19 6.06 6.57 5.75  5.65  5.84  5.32  6.67  5.36  5.43  5.62  5.57  5.76
#[,19] [,20] [,21] [,22] [,23] [,24] [,25]
#[1,]  6.72   5.7  6.24  6.89  5.79   5.5  5.62
# t(round(log(dist_RMDP), 2))[c(2,  8,  9, 13, 22, 23)]
#[1] 6.79 6.57 5.75 6.67 6.89 5.79
cutoff_RMDP <- RMDP_output$cutoff
log(cutoff_RMDP)
#[1] 6.273443
ou_RMDP <- RMDP_output$ou_RMDP
RMDP1 <- as.data.frame(t(dist_RMDP))
df2 <- as.data.frame(t(rbind(c(1:n), log(RMDP1))))
df2.highlight <- df2[c(2,  8,  9, 13, 22, 23), ]
df2.text <- df2.highlight
df2.text$label <- c("06-01:6.79", "07-21:6.57", "07-31:5.75", "09-01:6.67", "11-23:6.89", "12-01:5.79")
pic2 <-
  ggplot(df2, aes(V1, V2)) + geom_point(size = 3, shape = 15)+ 
  geom_point(data = df2.highlight, aes(V1, V2), size = 5, col="black")+
  geom_text_repel(data = df2.text, aes(V1, V2, label = label), color = "black", seed = 2, size = 5)+
  geom_hline(yintercept = log(cutoff_RMDP), linetype='dotted', col = 'darkcyan', size = 1)+
  annotate("text", x = 1, y = log(cutoff_RMDP), label = "cutoff:6.27", vjust = -0.5, size = 5)+
  theme(axis.line = element_line(color="black", size = 0.2), 
        axis.title = element_text(size = 18), plot.title = element_text(size = 18))+
  ylim(5, 11) + xlim(0, 25) + labs(x = "Index", y = "Modified distance", title = "Method RMDP")
pic2
pic <- cowplot::plot_grid(pic1, pic2, ncol=2)#Figure S2: Plots of the modified distances based on the RICD and the RMDP
pic

#Example 2
data(octane)
X <- as.matrix(octane)
RICD_output <- RICD(X, h=32, alpha=0.01)
n <- nrow(X)
p <- ncol(X)
dist_RICD <- RICD_output$dist_RICD
cutoff_RICD <- RICD_output$cutoff
ou_RICD <- RICD_output$ou_RICD
which(ou_RICD==1)
robpca = PcaHubert(octane, k=2, alpha=0.75, mcd=FALSE)
(outl.robpca = which(robpca@flag==FALSE)) #check the 'true' outliers
#---------------------------------------------------------------------------------------------------------
# Q-Q plot of RICD distance measure
#---------------------------------------------------------------------------------------------------------
Theta1 <- RICD_output$Theta1
Theta2 <- RICD_output$Theta2
RICD1 <- (dist_RICD/p-Theta1)*sqrt(p)/sqrt(2*Theta2)
nor.points <- qnorm(ppoints(n))
qqplot(nor.points, RICD1, cex=1.1, cex.axis=0.9, cex.lab=1.1, xlab='Normal quantile', ylab='Distance measure', mgp=c(2,1,0))
abline(0, 1, col = 'gray30')
abline(h=(cutoff_RICD/p-Theta1)*sqrt(p)/sqrt(2*Theta2), lwd=1, col='darkcyan', lty=2)
points(sort(nor.points)[34:39], sort(RICD1)[34:39], pch = 16, col = 'black')
# Test the same dataset with RMDP(Ro et al.2015) procedure
MDP_output <- MDP(X)
RMDP_output <- ReMDP(MDP_output, 0.01, X)
ou_RMDP <- RMDP_output$ou_RMDP
which(ou_RMDP==1)

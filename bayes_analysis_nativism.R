### Bayes chain for nativism data

library(GGally)


## inverse code Q4
data_transformed$trans_Q9_4 <-  6 - data_transformed$trans_Q9_4
####
### data_bayes

data_bayes = data_transformed %>%
  mutate(nativism = ((trans_Q9_1+trans_Q9_2+trans_Q9_3+trans_Q9_4+trans_Q9_5)/5)) %>%
  select(country, nativism)

data_us1 = data_US %>%
  select(trans_Q9_1, trans_Q9_2, trans_Q9_3, trans_Q9_4, trans_Q9_5)

data_test = data_transformed %>%
  select(trans_Q9_1, trans_Q9_2, trans_Q9_3, trans_Q9_4, trans_Q9_5)

colMeans(data_test)

apply(data_test, 2, sd)



empirical_mean = data_bayes %>%
  group_by(country) %>%
  summarise(empirical = mean(nativism))
hist(empirical_mean$empirical)


hist(data_bayes$nativism)


nu0<-1  ; s20<-1
eta0<-1 ; t20<-1
mu0<-3 ; g20<-2
###
Y = data_bayes
Y = as.data.frame(Y)
### starting values
m<-(unique(Y[, 1], drop = TRUE))
m = as_vector(m)
n<-sv<-ybar<-rep(NA,count(m)) 

Y$country = droplevels(Y$country)

for(j in 1:length(m)){
  ybar[j] = colMeans(Y[Y[,1] == m[j],2, drop = FALSE])
  sv[j]<-var(Y[Y[,1] == m[j], 2])
  n[j]<-sum(Y[,1] == m[j]) 
}

ybar = ybar[1:22]
sv = sv[1:22]
n = n[1:22]

######colMeans(Y[Y[,1] == "US",2], na.rm = TRUE)
##ybar[3]
##var(Y[Y[,1] == "US",2], na.rm = TRUE)
##sv[3]
##sum(Y[,1] == "US")
##n[3]
##

###
theta<-ybar
sigma2<-mean(sv)
mu<-mean(theta)
tau2<-var(theta)
###

### setup MCMC
set.seed(1)
S<-5000
THETA<-matrix( nrow=S,ncol=length(m))
MST<-matrix( nrow=S,ncol=3)
###

### MCMC algorithm
for(s in 1:S) 
{
  
  # sample new values of the thetas
  for(j in 1:length(m)) 
  {
    vtheta<-1/(n[j]/sigma2+1/tau2)
    etheta<-vtheta*(ybar[j]*n[j]/sigma2+mu/tau2)
    theta[j]<-rnorm(1,etheta,sqrt(vtheta))
  }
  
  #sample new value of sigma2
  nun<-nu0+sum(n)
  ss<-nu0*s20
  for(j in length(m)){ss<-ss+sum((Y[Y[,1] == m[j], 2]-theta[j])^2)}
  sigma2<-1/rgamma(1,nun/2,ss/2)
  
  #sample a new value of mu
  vmu<- 1/(length(m)/tau2+1/g20)
  emu<- vmu*(length(m)*mean(theta)/tau2 + mu0/g20)
  mu<-rnorm(1,emu,sqrt(vmu)) 
  
  # sample a new value of tau2
  etam<-eta0+length(m)
  ss<- eta0*t20 + sum( (theta-mu)^2 )
  tau2<-1/rgamma(1,etam/2,ss/2)
  
  #store results
  THETA[s,]<-theta
  MST[s,]<-c(mu,sigma2,tau2)
} 

plot(density(MST[, 1]))

#######
stationarity.plot<-function(x,...){
  S<-length(x)
  scan<-1:S
  ng<-min(round(S/100),10)
  group<-S*ceiling( ng*scan/S) /ng
  boxplot(x~group,...)
}

##########

stationarity.plot(MST[,1],xlab="iteration",ylab=expression(mu))
stationarity.plot(MST[,2],xlab="iteration",ylab=expression(sigma^2))
stationarity.plot(MST[,3],xlab="iteration",ylab=expression(tau^2))

##########

pdf("nativism_bayes.pdf",family="Times",height=1.75,width=5)
par(mfrow=c(1,3),mar=c(2.75,2.75,.5,.5),mgp=c(1.7,.7,0))
plot(density(MST[,1],adj=2),xlab=expression(mu),main="",lwd=2,
     ylab=expression(paste(italic("p("),mu,"|",italic(y[1]),"...",italic(y[m]),")")))
abline( v=quantile(MST[,1],c(.025,.5,.975)),col="gray",lty=c(3,2,3) )
plot(density(MST[,2],adj=2),xlab=expression(sigma^2),main="", lwd=2,
     ylab=expression(paste(italic("p("),sigma^2,"|",italic(y[1]),"...",italic(y[m]),")")))
abline( v=quantile(MST[,2],c(.025,.5,.975)),col="gray",lty=c(3,2,3) )
plot(density(MST[,3],adj=2),xlab=expression(tau^2),main="",lwd=2,
     ylab=expression(paste(italic("p("),tau^2,"|",italic(y[1]),"...",italic(y[m]),")")))
abline( v=quantile(MST[,3],c(.025,.5,.975)),col="gray",lty=c(3,2,3) )
dev.off()

####
plot(density(THETA[,3]))

plot(density(THETA[,3]), main = "Posterior distribution of the mean of the US vs India",
     ylim=c(-.05,80), xlim = c(3.3, 3.4),ylab="",yaxt="n")


#points(Y[Y[,1] == m[3],]$nativism,rep(0.01666,1000),col="black",pch=16)

points( ybar[3],-.01666,col="black",pch=16 ,cex=1.5)
abline( h=-.01666,col="black")
lines(density(THETA[,4],adj=2),col="gray",lwd=2)
lines(density(THETA[,6],adj=2),col="red", lwd=2)

mean(THETA[,3] < THETA[,4])



##########
par(mfrow=c(1,1))
plot(density(THETA[,3],adj=2),col="black",xlim=
       range(c(Y[Y[,1] == 3, 2],Y[Y[,1] == 16, 2],THETA[,c(3,16)])),lwd=2,
     main="",xlab="math score",ylim=c(-.05,8),ylab="",yaxt="n")

axis(side=2,at=c(0,0.10,0.20) )
lines(density(THETA[,16],adj=2),col="gray",lwd=2)
abline(h=0)

points( Y[Y[,1] == m[3], 2],rep(-0.01666,n[3]), col="black",pch=16)
points( ybar[3],-.01666,col="black",pch=16 ,cex=1.5)
abline( h=-.01666,col="black")

points( Y[Y[,1] == m[16], 2],rep(-0.0333,n[16]), col="gray",pch=16)
points( ybar[16],-.0333,col="gray",pch=16 ,cex=1.5)
abline( h=-.0333,col="gray")

segments(mean(MST[,1]), 0,mean(MST[,1]),1,lwd=2,lty=2 )

legend(52.5,.15,legend=c("school 46","school 82",
                         expression(paste("E[", mu,"|",italic(y[1]),"...",italic(y[m]),"]"))),
       lwd=c(2,2),lty=c(1,1,2),col=c("black","gray"),bty="n")
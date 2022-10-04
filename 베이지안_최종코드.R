#########################Gruop 3 _ SeoSuMin, KimPuReum, JungYuRim
#HPD.grid function
HPDgrid=function(prob,level=0.95){
  
  prob.sort=sort(prob,decreasing = T)
  
  M=min(which(cumsum(prob.sort)>=level))
  
  height=prob.sort[M]
  
  HPD.index=which(prob>=height)
  
  HPD.level=sum(prob[HPD.index])
  
  res=list(index=HPD.index,level=HPD.level)
  
  return(res)
  
}

######Binomial Distribution######
install.packages("HSAUR")
library(HSAUR)
data("mastectomy")

##prior distribution : theta~Beta(1,1)
a=1;b=1
n=40; x=24
theta<-seq(0,1,length=100)
prior.theta<-dbeta(theta,a,b)

##posterior Distribution : theta|x~Beta(25,17)
posterior.theta<-dbeta(theta,a+x,b+n-x)

x/n #Classic estimate
a+x/(b+n-x)#Bayesian estimate

##Predictive Distribution
m=4; z=c(0:4)
pred.z<-gamma(m+1)/gamma(z+1)/gamma(m-z+1)*beta(a+z+x,b+n-x+m-z)/beta(a+x,b+n-x)
plot(z,pred.z,type="h")

pred.z #P(z=2) = P(z=3) = 0.33373603 

m*(x/n) #Classic predictive value
sum(z*pred.z) #Bayesian predictive value


##Credible Region
#HPD.beta function
HPD.beta=function(a,b,guesses,C,niter){
  con=lgamma(a+b)-lgamma(a)-lgamma(b)-log(C)
  new1=guesses[1]
  new2=guesses[2]
  for(i in 1:niter)
    new1=new1-(con+(a-1)*log(new1)+(b-1)*log(1-new1))/((a-1)/new1-(b-1)/(1-new1))
  for(i in 1:niter)
    new2=new2-(con+(a-1)*log(new2)+(b-1)*log(1-new2))/((a-1)/new2-(b-1)/(1-new2))
  p1=pbeta(new1,a,b)
  p2=1-pbeta(new2,a,b)
  c(new1,new2,p1,p2,p1+p2)
}

#1.Classic Credible Region
p=x/n
l_origin=p-1.96*sqrt((p*(1-p))/n)
u_origin=p+1.96*sqrt((p*(1-p))/n)
c(l_origin,u_origin) #0.4481791 0.7518209

plot(theta,dbeta(theta,a+x,n-x+b),type='l')
abline(v=c(l_origin,u_origin),lty=2,col=4)
#legend(0.2,4,"고전적",lty=2,col=4)

dbeta(c(l_origin,u_origin),a+x,b+n-x) #y (0.8216973,  0.5676645)
pbeta(u_origin,a+x,b+n-x)-pbeta(l_origin,a+x,b+n-x)#Credible Region area = (0.9583567)
u_origin-l_origin #Credible Region length = (0.3036419)

#2.HPDgrid method
ftheta=dbeta(theta,a+x,n-x+b)
prob=ftheta/sum(ftheta)
HPD=HPDgrid(prob,0.95)
HPD1.grid=c(min(theta[HPD$index]),max(theta[HPD$index]))
HPD1.grid#Credible Region (0.4545455 , 0.7373737)

plot(theta,dbeta(theta,a+x,n-x+b),type='l')
abline(v=HPD1.grid,lty=2,col=6)

dbeta(HPD1.grid,25,17)#y (0.9574065 , 0.8810856)
pbeta(HPD1.grid[2],a+x,b+n-x)-pbeta(HPD1.grid[1],a+x,b+n-x)#Credible Region area (0.9423497)
HPD1.grid[2]-HPD1.grid[1]#Credible Region length (0.2828283)

#3.HPD.beta
qu=qbeta(c(0.025,0.975),a+x,b+n-x)
q1=qu[1]; q2=qu[2]

for(k in seq(5,0,length=10000)){
  out=HPD.beta(a+x,b+n-x,c(q1,q2),k,100)
  diff=pbeta(out[2],a+x,b+n-x)-pbeta(out[1],a+x,b+n-x)
  if(diff>=0.95) break
}
c(out[1],out[2])#Credible Region(0.4482229 , 0.7397671)

plot(theta,dbeta(theta,a+x,b+n-x),type='l')
abline(v=c(out[1],out[2]),lty=2,col=2)
abline(h=k,lty=2,col=2)
mtext(round(k,3),side=2,at=k,col=2)

dbeta(c(out[1],out[2]),a+x,b+n-x) #y ( 0.8225823 , 0.8225823) 
pbeta(out[2],a+x,b+n-x)-pbeta(out[1],a+x,b+n-x)#Credible Region area(0.9500073)
out[2]-out[1]#HPD Credible Region length(0.2915442)

#4.quantile method
lc1=qbeta(0.025,a+x,b+n-x)
uc1=qbeta(0.975,a+x,b+n-x)
c(lc1,uc1)#Credible Region(0.4450478 , 0.7368320)

plot(theta,dbeta(theta,a+x,n-x+b),type='l')
abline(v=c(lc1,uc1),lty=2,col=3)
#legend(0.1,4,c("고전적","격자점","사후분위수","HPD"),lty=c(2,2,2,2),col=c(4,6,3,2))

dbeta(c(lc1,uc1),a+x,b+n-x)#y ( 0.7602439 , 0.8946978)
pbeta(uc1,a+x,b+n-x)-pbeta(lc1,a+x,b+n-x)#Credible Region area (0.95)
uc1-lc1#quantile Credible Region length(0.2917842)

##Hypothesis Testing
#H0: theta<=0.4 H1: theta>0.4

#Classical Hypothesis Testing
binom.test(24,40,p=0.4,alternative = "greater")
sum(dbinom(c(24:40),40,0.4))# p-value=0.008342 <0.05. Reject to H0

#Baysian Hypothesis Testing
pie0=pbeta(0.4,a,b)
pie1=1-pie0
p0=pbeta(0.4,a+x,b+n-x)
p1=1-p0

B=(p0/p1)/(pie0/pie1)
B #Bayes factor= 0.008063959 < 1. Reject to H0



######Poisson Distribution######

##prior distribution: theta~Gamma(3,0.22)
x<-c(13,4,8,7,11,18,1,42,12,18,17,11,25,27,5,1,18,8)#prior data
hist(x,breaks = 5)

x1=seq(0,50,length=1000)
plot(x1,dgamma(x1,3,0.22),type='l')

##posterior Distribution: theta|x~Gamma(47,6.22)
x1<-c(4,4,17,7,9,3)
a=3; b=0.22
n1=length(x1); s1=sum(x1)
theta<-seq(0,15,length=100)
plot(theta,dgamma(theta,a+s1,b+n1),type='l',ylab="p(theta|x)")


mean(x1)#Classic estimate
(postmean.theta1=(a+s1)/(b+n1))#Bayesian estimate

##Predictive Distribution
x=seq(0,20)
r<-a+s1
theta100<-(b+n1)/(b+n1+1)
plot(x,dnbinom(x,size=r,prob=theta100),type='h')

mean(x1)#Classic predictive value
(r*(1-theta100))/theta100 #Bayesian predictive value

##Credible Region
#1.HPDgrid
theta=seq(0,15,length=1000)
ftheta=dgamma(theta,a+s1,b+n1)
prob=ftheta/sum(ftheta)
HPD1=HPDgrid(prob,0.95)
HPD1.grid=c(min(theta[HPD1$index]),max(theta[HPD1$index]))
HPD1.grid #Credible Region(5.465465, 9.744745)

plot(theta,dgamma(theta,a+s1,b+n1),type='l',main="격자점 방법을 이용한 95% 신뢰구간")
abline(v=HPD1.grid,lty=2,col=6)

pgamma(HPD1.grid[2],a+s1,b+n1)-pgamma(HPD1.grid[1],a+s1,b+n1)#Credible Region area (0.9493171)
dgamma(HPD1.grid,a+s1,b+n1) #y (0.05426275, 0.05338597)
HPD1.grid[2]-HPD1.grid[1] #Credible Region length(4.279279)

#2&3 HPD & Classic Credible Regionr
x1<-c(4,4,17,7,9,3)
a=3; b=0.22
n1=length(x1); s1=sum(x1)
theta<-seq(0,15,length=100)

HPD.gamma=function(a,b,guesses,C,niter){ 
  con=a*log(b)-lgamma(a)-log(C)
  new1=guesses[1]
  new2=guesses[2]
  for(i in 1:niter){
    new1=new1-(con+(a-1)*log(new1)-b*new1)/(((a-1)/new1)-b)
  }
  for(i in 1:niter){
    new2=new2-(con+(a-1)*log(new2)-b*new2)/(((a-1)/new2)-b)
  }
  p1=pgamma(new1,a,b)
  p2=1-pgamma(new2,a,b)
  c(new1,new2,p1,p2,p1+p2)
}

qu=qgamma(c(0.025,0.975),a+s1,b+n1); dgamma(qu[1],a+s1,b+n1); dgamma(qu[2],a+s1,b+n1);
pgamma(qu[2],a+s1,b+n1)-pgamma(qu[1],a+s1,b+n1)
q1=qu[1]; q2=qu[2]

for(k in seq(0.4, 0, length=10000)){
  out=HPD.gamma(a+s1,b+n1,c(q1,q2),k,4)
  diff=pgamma(out[2],a+s1,b+n1)-pgamma(out[1],a+s1,b+n1)
  if(diff>=0.95) break
}
c(out[1],out[2])#HPD Credible Region (5.455853, 9.748007)

plot(theta,dgamma(theta,a+s1,b+n1),type='l',main="HPD & 고전적 95% 신뢰구간")
abline(v=c(out[1],out[2]),lty=2,col=2)
abline(h=k,lty=2,col=2)
abline(v=c(q1,q2),lty=3,col=3)
mtext(round(k,3),side=2,at=k,col=2)

pgamma(out[2],a+s1,b+n1)-pgamma(out[1],a+s1,b+n1)#HPD Credible Region area(0.9500069)
dgamma(c(out[1],out[2]),a+s1,b+n1)#HPD y (0.05312531, 0.05312531)
out[2]-out[1]#HPD Credible Region length(4.292154)

c(q1,q2)#Classic Credible Region( 5.552063 , 9.864558)
pgamma(q2,a+s1,b+n1)-pgamma(q1,a+s1,b+n1)# Classic Credible Region area(0.95)
dgamma(c(q1,q2),a+s1,b+n1)#Classic y (0.06525707, 0.04445323)
q2-q1#Credible Region length (4.312495)

#4.quantile method 
lc=qgamma(0.025,a+s1,b+n1)
uc=qgamma(0.975,a+s1,b+n1)
c(lc,uc) #Credible Region (5.552063,9.864558)

plot(theta,dgamma(theta,a+s1,b+n1),type='l',main="사후분위수를 이용한 95% 신뢰구간")
abline(v=c(lc,uc),lty=2,col=4)
#legend(1,0.2,c("고전적","격자점","사후분위수","HPD"),lty=c(3,2,2,2),col=c(3,6,4,2))

pgamma(uc,a+s1,b+n1)-pgamma(lc,a+s1,b+n1)#Credible Region area (0.95)
dgamma(c(lc,uc),a+s1,b+n1)#y (0.06525707, 0.04445323)
uc-lc# Credible Region length(4.312495)

##Hypothesis Testing
#H0: theta>=10 H1: theta<10

#Classical Hypothesis Testing
sum(dpois(c(0:44),6*10))
poisson.test(44,T=6,r=10,alternative = "less", conf.level = 0.95) #p-value=0.01897 < 0.05. Reejct to H0

#Baysian Hypothesis Testing
Pi1=pgamma(10,a,b)
Pi0=1-Pi1
P1=pgamma(10,a+s1,b+n1)
P0=1-P1
B01=(P0/P1)/(Pi0/Pi1)
B01#Bayes factor=0.01209929 < 1. Reject to H0

######Normal Distribution######
install.packages("Stat2Data")
library(Stat2Data)
data('Ricci')

set.seed(1209)
k=sample(1:118,5)
pri_Ricci=Ricci[-k,]
predict_Ricci=Ricci[k,]
#Random sampling

##Normal data proof
#H0:Follow the normal distribution. H1:Not follow the normal distribution

qqnorm(pri_Ricci$Written)
qqline(pri_Ricci$Written)
#qqplot

shapiro.test(pri_Ricci$Written)
#shapiro.test result: p-value > 0.05
#Fail to reject H0. Follow the normal distribution

##prior distribution: pie(theta)=1
hist(pri_Ricci$Written)

mean(pri_Ricci$Written) #Mean=71.57522
sd(pri_Ricci$Written) #sd=10.75025

mean(predict_Ricci$Written) #Mean=73.4

##posterior Distribution: theta|x~N(71.57522, 1.011299^2)
n=113; xbar=mean(pri_Ricci$Written); s=sd(pri_Ricci$Written)
mu.post=xbar; sig.post=sqrt(s^2/n) #Classic estimate = Bayesian estimate (71.57522) 

theta=seq(mu.post-5*sig.post,mu.post+5*sig.post,length=100)
plot(theta,dnorm(theta,mu.post,sig.post),type='l',main="posterior of theta")

##Predictive Distribution
mu.new=mu.post; sig.new=sqrt(s^2+sig.post^2)
xnew=seq(mu.new-5*sig.new,mu.new+5*sig.new,length=100)

plot(xnew,dnorm(xnew,mu.new,sig.new),type='l',main="Predictive density of X_new") #Predictive Distribution~N(71.57522, 10.79772^2)

mu.post#Classic predictive value(71.57522)
mu.new #Bayesian predictive value(71.57522)

##Credible Region
#1.Classical Hypothesis Testing
lc2=mean(pri_Ricci$Written)-1.96*sd(pri_Ricci$Written)
uc2=mean(pri_Ricci$Written)+1.96*sd(pri_Ricci$Written)
c(lc2,uc2)#Credible Region (50.50472, 92.64572)

theta=seq(mean(pri_Ricci$Written)-5*sd(pri_Ricci$Written),mean(pri_Ricci$Written)+5*sd(pri_Ricci$Written),length=100)
plot(theta,dnorm(theta,mean(pri_Ricci$Written),sd(pri_Ricci$Written)),type='l',main="고전적 방법을 이용한 95% 신뢰구간")
abline(v=c(lc2,uc2),lty=2,col=4)
#legend(100,0.02,"고전적",lty=2,col=4)

pnorm(uc2,mean(pri_Ricci$Written),sd(pri_Ricci$Written))-pnorm(lc2,mean(pri_Ricci$Written),sd(pri_Ricci$Written))#Credible Region area(0.9500042)
dnorm(c(lc2,uc2),mu.post,sig.post)#y (2.149016e-95, 2.149016e-95)
uc2-lc2 #Credible Region length(42.14099)

#2.HPDgrid method
theta=seq(mu.post-5*sig.post,mu.post+5*sig.post,length=100)
ftheta=dnorm(theta,mu.post,sig.post)
prob=ftheta/sum(ftheta)
HPD1=HPDgrid(prob,0.95)
HPD1.grid=c(min(theta[HPD1$index]),max(theta[HPD1$index]))
HPD1.grid #Credible Region (69.58327, 73.56717)

plot(theta,dnorm(theta,mu.post,sig.post),type='l',main="격자점 방법을 이용한 95% 신뢰구간")
abline(v=HPD1.grid,lty=2,col=1)

pnorm(HPD1.grid[2],mu.post,sig.post)-pnorm(HPD1.grid[1],mu.post,sig.post)#Credible Region area(0.9511269)
dnorm(HPD1.grid,mu.post,sig.post)#y (0.05669739, 0.05669739)
HPD1.grid[2]-HPD1.grid[1]#Credible Region length (3.983904)

#3.quantile method = HPD
lc=qnorm(0.025,mu.post,sig.post)
uc=qnorm(0.975,mu.post,sig.post)
c(lc,uc) #Credible Region (69.59311, 73.55733)

plot(theta,dnorm(theta,mu.post,sig.post),type='l',main="사후분위수를 이용한 95% 신뢰구간")
abline(v=c(lc,uc),lty=2,col=2)
#legend(74,0.2,c("격자점","사후분위수"),lty=c(2,2),col=c(1,2))


pnorm(uc,mu.post,sig.post)-pnorm(lc,mu.post,sig.post)#Credible Region area(0.95)
dnorm(c(lc,uc),mu.post,sig.post) #y (0.05779209, 0.05779209)
uc-lc #Credible Region length (3.964218)

##Hypothesis Testing
#H0:theta<=72  H1:theta>72

#Classical Hypothesis Testing
z=(xbar-72)/(s/sqrt(n))
pvalue=1-pnorm(z)
pvalue #p-value = 0.6627693 > 0.05. Fail to reject H0

#Baysian Hypothesis Testing
p0=pnorm((72-mu.post)/(s/sqrt(n))) #0.6627693
p1=1-p0 #0.3372307

#prior distribution ~N(0,sigma^2) : unknown variance >>> use t-distribution
pi0=pt((xbar-72)/(s/sqrt(n)),df=112) #0.3376331
pi1=1-pi0 #0.6623669

#Bayes factor
B=(p0/p1)/(pi0/pi1)
B #Bayes factor = 3.855573 > 1. Fail to reject H0

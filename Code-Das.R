library( Rcpp)
library(rstanarm)
library(HDInterval)
library(pgdraw)
library(BayesLogit)
library(truncnorm) #for Truncated normal sample 
library(MASS) # for Multivariate normal sample
setwd("C:/Users/apust/Desktop/Bayesian methods/Project")
dat=read.csv("brain_stroke.csv", header=TRUE)
summary(dat) 
str(dat) 
attach(dat)
##########data processing
dat$gender=ifelse(dat$gender=="Female", 0, 1) #female=0, male=1
dat$ever_married=ifelse(dat$ever_married=="No", 0, 1) #no=0, yes=1
dat$residence_type=ifelse(dat$residence_type=="Rural", 0, 1) #rural=0, urban=1
dat$work_type=as.numeric(factor(work_type)) #Private=3, Self-employed=4, Govt_job=2, children=1
dat$smoking_status=as.numeric(factor(dat$smoking_status)) #"formerly smoked"=1 "never smoked"=2 "smokes"=3 "unknown"=4
dat1=dat[which(dat$stroke==1),]
dat2=dat[which(dat$stroke==0 & dat$age>38),]
newdat=rbind(dat1, dat2)
attach(newdat)

#logit model#######
#mle.logit<-glm(stroke~gender+age+hypertension+heart_disease+ever_married+work_type+residence_type+
#avg_glucose_level+bmi+smoking_status,family=binomial(link="logit"))
#stepAIC(mle.logit,direction="backward") #using backward selection method
mle.logit<-glm(stroke~age+hypertension+
                  avg_glucose_level,family=binomial(link="logit"))
hat.Beta.MLE<-as.numeric(summary(mle.logit)$coefficients[,1])
se.Beta.MLE<-as.numeric(summary(mle.logit)$coefficients[,2])
round(hat.Beta.MLE,3)
round(se.Beta.MLE,3)
round(confint(mle.logit),3)

exp(cbind(OR=coef(mle.logit), confint(mle.logit)))

#Bayesian logistic regression
n=dim(newdat)[1]
x<-matrix(1,n,4)
x[,2:4]<-c(age, hypertension, avg_glucose_level)
y=stroke
f.stroke<- factor(stroke)
# preparing the inputs for the model
#m <- model.matrix(f.stroke ~ age+hypertension+heart_disease+ever_married+
                    #avg_glucose_level, data = newdat)
#m=model.matrix(f.stroke~gender+age+hypertension+heart_disease+ever_married+work_type+residence_type+
                 #avg_glucose_level+bmi+smoking_status, data=newdat)

#Will use N(1,1) as prior
posterior3 <- stan_glm(f.stroke~age+hypertension+avg_glucose_level, data=newdat,
                       family = binomial(link = "logit"), 
                      prior = normal(1, 1), QR=TRUE,
                       refresh = 0, chains = 10, iter = 2000)

posterior3$ses #MAD_SD
cof=round(coef(posterior3), 3)
exp(cof)#odds ratios for coefficients
CI=round(posterior_interval(posterior3, prob = 0.95,type = "central"), 3)#Posterior uncertainty Intervals or credible intervals
round(exp(cbind(cof, CI)),3) #OR and interval
#hdi(posterior3, ci =0.95) #95% High Density Intervals (DHI)
plot(posterior3, "areas", prob = 0.95, prob_outer = 1)#Plotting posterior distribution


#linpred <- posterior_linpred(posterior3)
#preds <- posterior_linpred(posterior3, transform=TRUE)
#pred <- colMeans(preds)
#pr <- as.integer(pred >= 0.5)

# posterior classification accuracy by probability
#round(mean(xor(pr,as.integer(y==0))),2)


##Gibbs sampler for logistic regression
#Prior=N(1,1)
post.mode<-mle.logit$coefficients
var.prior=diag(4)
mean.prior=matrix(1, nrow=4, ncol=1)
iter=10000
MC.beta<-matrix(0,iter,4);w<-rep(0,iter)
for (i in 1:iter){
  #w|beta,y
  w.var=x%*%post.mode
  #w=pgdraw(1,w.var) #pgdraw(b,c)
  w=rpg(n,1,w.var) ##confused about 1st parameter (n or 1)
  
  #beta|w,y
  omega=diag(w)
  beta.var=solve(t(x)%*%omega%*%x+var.prior) #conditional variance
  kappa=y-1/2
  beta.mean=beta.var%*%(t(x)%*%kappa+var.prior%*%mean.prior) #conditional variance
  post.mode=mvrnorm(1,beta.mean,beta.var)
  MC.beta[i,]<-post.mode
  
}
MCMC.beta=MC.beta[-(1:5000),]
est.beta=round(apply(MCMC.beta, 2, mean),3) #posterior mean
est.se=round(apply(MCMC.beta, 2, sd),3) #posterior SE
est.CI=round(posterior_interval(MCMC.beta, prob=0.95),3) #95% credible intervals

#trace plot
par(mfrow=c(2,2))
ts.plot(MCMC.beta[,1],ylab="Intercept")
ts.plot(MCMC.beta[,2],ylab="Age")
ts.plot(MCMC.beta[,3],ylab="Hypertension")
ts.plot(MCMC.beta[,4],ylab="Avg glucose level")

#autocorrelation function
par(mfrow=c(2,2))
acf(MCMC.beta[,1],main="Intercept")
acf(MCMC.beta[,2],main="Age")
acf(MCMC.beta[,3], main="Hypertension")
acf(MCMC.beta[,4],main="Avg glucose level")
#CPO and LPML
x.MC.beta=x%*%t(MCMC.beta)
prob.MC = 1/(1+exp(-x.MC.beta))
y=newdat$stroke
CPO.logit=round(1/apply(1/(prob.MC*y+(1-prob.MC)*(1-y)),1,mean),3)
LPML.logit=sum(log(CPO.logit))


########Bayesian Probit model############
neg.log.like<-function(b){
  log.out<-(-1)*y*log(pnorm(x%*%b))+(-1)*(1-y)*log(1-pnorm(x%*%b))
  return(sum(log.out))}
Optim.post<-optim(c(0,0,0,0),neg.log.like,method = c("BFGS"),hessian=TRUE)

post.mode<-Optim.post$par

##########Gibbs sampler for Probit model ############
###########################################
R<-10000
MC.beta<-matrix(0,R,4)
hat.beta<-post.mode #initial values (we can choose any values)
hat.z<-rep(0,n)
for(j in 1:R){
  #Z|beta,y ##########
  for(i in 1:n){
    tx_ibeta<-as.numeric(t(x[i,])%*%hat.beta)
    if(y[i]==1){
      hat.z[i]<-rtruncnorm(1, a=0, b=Inf, mean =tx_ibeta , sd = 1)
    }
    if(y[i]==0){
      hat.z[i]<-rtruncnorm(1, a=-Inf, b=0, mean = tx_ibeta, sd = 1)
    }
  }
  #beta|Z,y##########
  inv.tXX<-solve(t(x)%*%x)
  mean.beta<-as.numeric(inv.tXX%*%t(x)%*%hat.z)
  hat.beta<-mvrnorm(1,mean.beta,inv.tXX)
  MC.beta[j,]<-hat.beta
}

MCMC.beta=MC.beta[-(1:5000),]
est.beta=round(apply(MCMC.beta, 2, mean),3) #posterior mean
est.se=round(apply(MCMC.beta, 2, sd),3) #posterior SE
est.CI=round(posterior_interval(MCMC.beta, prob=0.95),3) #95% CI

#trace plot
par(mfrow=c(2,2))
ts.plot(MCMC.beta[,1],ylab="Intercept")
ts.plot(MCMC.beta[,2],ylab="Age")
ts.plot(MCMC.beta[,3],ylab="Hypertension")
ts.plot(MCMC.beta[,4],ylab="Avg glucose level")

#autocorrelation function
par(mfrow=c(2,2))
acf(MCMC.beta[,1],main="Intercept")
acf(MCMC.beta[,2],main="Age")
acf(MCMC.beta[,3], main="Hypertension")
acf(MCMC.beta[,4],main="Avg glucose level")

###CPO and LPML
x.MC.beta=x%*%t(MCMC.beta)
prob.MC = pnorm(x.MC.beta,0,1)
y=newdat$stroke
CPO.probit=round(1/apply(1/(prob.MC*y+(1-prob.MC)*(1-y)),1,mean),3)
LPML.probit=sum(log(CPO.probit))



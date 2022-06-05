#Program to calculate statistical properties of phos estimates using maximum likelihood theory
#and Monte Carlo Simulation. This code allows inputs from several hatcheries with
#potentially different visual marking (VM) fractions and different parentage-based tagging (PBT) fractions.

#AUTHOR: Richard A. Hinrichsen, 10 July 2013
#CONTACT: rich@hinrichsenenvironmental.com

#Variables and parameters used in the analysis
#inputs
#phosi = true proportions of hatchery origin spawners (hatchery-specific)
#Nsamp = total number of spawners sampled on spawning grounds
#n = total number of spawners tested for PBT
#n1 = number of visually marked spawners tested for PBT (when OPT=FALSE)
#lambda = marking fraction (hatchery-specific)
#ppbt = fraction of fish that are PBT (hatchery-specific)
#OPT = FALSE when n1 is user input, TRUE when program to select an optimal value of n1
#MONTE = FALSE for theoretical results, TRUE for Monte Carlo results
#NSIM = number of Monte Carlo simulations (needed if MONTE=TRUE)
#
#
#Select intermediate variables
#nhatch = number of hatcheries supplying spawners to spawning grounds
#I = Fisher Information Matrix
#x1 = number of visually marked spawners in sample of size Nsamp
#x2 = number of unmarked spawners in sample of size Nsamp
#n2 = number of unmarked spawners tested for PBT
#y = number of visually marked spawners tested that were PBT (hatchery-specific)
#z = number of unmarked spawners tested that were PBT (hatchery-specific)
#Ey = expected value of y
#Ez = expected value of z

#Results
#phos = true proportion of hatchery-origin spawners
#Ex1 = expected number of visually marked spawners (summing over hatcheries)
#Ex2 = expected number of not visually marked spawners (summing over hatcheries)
#SE_MIN.phos = standard error (SE) when all sampled fish are tested for PBT
#CV_MIN.phos = Coefficient of variation when ALL sampled fish are tested for PBT
#SE.phos = standard error (SE)
#CV.phos = Coefficient of variation
#BIAS.phos = relative bias estimate (NA if MONTE=FALSE)
#n1 = optimal number of visually marked spawners tested for PBT (when OPT=TRUE)

#top level function
phos.pbt.main<-function(phosi=.1*c(1/20,1/20,9/20,9/20),
Nsamp=1000,n=200,n1=50,
lambda=c(1,.95,.5,.5),
ppbt=c(.95,.95,.95,.95),
OPT=FALSE,
MONTE=FALSE,NSIM=1000){
check.inputs(phosi,Nsamp,n,n1,lambda,ppbt,OPT,MONTE,NSIM)
if(!OPT){
if(!MONTE){res<-phos.pbt.estimates(phosi=phosi,Nsamp=Nsamp,n=n,n1=n1,lambda=lambda,ppbt=ppbt,suppress=FALSE)}
if(MONTE){res<-phos.pbt.estimates2(NSIM=NSIM,phosi=phosi,Nsamp=Nsamp,n=n,n1=n1,lambda=lambda,ppbt=ppbt)}}
if(OPT){res<-optimize(phosi=phosi,Nsamp=Nsamp,n=n,lambda=lambda,ppbt=ppbt)}
final.res<-list(OPT=OPT,
MONTE=res$MONTE,
NSIM=res$NSIM,
phosi=res$phosi,
Nsamp=res$Nsamp,
n=res$n,
n1=res$n1,
lambda=res$lambda,
ppbt=res$ppbt,
phos=res$phos,
Ex1=res$Ex1,
Ex2=res$Ex2,
SE_MIN.phos=res$SE_MIN.phos,
CV_MIN.phos=res$CV_MIN.phos,
SE.phos=res$SE.phos,
CV.phos=res$CV.phos,
BIAS.phos=res$BIAS.phos)
return(final.res)
}

#check of feasibility of optimization
#avoid unusual case where hatcheries with zero expected tags recoveries
#and these hatcheries do not all use the same visible marking fraction
is.feas<-function(phosi,Nsamp,n,n1,lambda,ppbt){
n2<-n-n1
onelambda2<-TRUE
iii<-n1*lambda*ppbt+n2*(1-lambda)*ppbt==0
onelambda2<-sum(mean(lambda[iii])==lambda[iii])==sum(iii)
onelambda2<-onelambda2&(mean(lambda[iii])>0)
if((!onelambda2)&sum(iii))return(FALSE)
return(TRUE)
}

#Find the value of n1 that minimizes CV
optimize<-function(phosi,Nsamp,n,lambda,ppbt){
#loop over all feasible values of n1
#note that n1 might be a non-integer because the constraints are not necessarily integers
Ex1<-Nsamp*sum(lambda*phosi)
min.n1<-max(Ex1-(Nsamp-n),0)
max.n1<-min(Ex1,n)
min.int<-ceiling(min.n1)
max.int<-floor(max.n1)
if(min.int>=max.int){n1s<-min.n1}
if(min.int<max.int){n1s<-c(min.n1,min.int:max.int,max.n1)}
n1s<-unique(n1s)
nfeas<-length(n1s)
cv<-1.e10
icount<-0
for(ii in 1:nfeas){
if(is.feas(phosi=phosi,Nsamp=Nsamp,n=n,n1=n1s[ii],lambda=lambda,ppbt=ppbt)){
icount<-icount+1
res.new<-phos.pbt.estimates(phosi=phosi,Nsamp=Nsamp,n=n,n1=n1s[ii],lambda=lambda,ppbt=ppbt,suppress=TRUE)
if(res.new$CV.phos<cv){res<-res.new;n1<-res.new$n1;cv<-res.new$CV.phos}}
}
if(icount==0)stop("Error in optimize: no feasible values are available for n1")
return(res)
}

#make sure the inputs make sense
check.inputs<-function(phosi,Nsamp,n,n1,lambda,ppbt,OPT,MONTE,NSIM){
if(!is.logical(OPT))stop("OPT must be TRUE or FALSE")
if(!is.logical(MONTE))stop("MONTE must be TRUE or FALSE")
if(OPT){
if(MONTE){stop("For the optimization option, the program cannot be run in Monte Carlo mode")}}
if(MONTE){
if(floor(NSIM)!=NSIM){stop("NSIM must be a positive integer")}
if(NSIM<=0){stop("NSIM must be a positive integer")}}
if(floor(Nsamp)!=Nsamp){stop("Nsamp must be a positive integer")}
if(Nsamp<=0){stop("Nsamp must be a positive integer")}
#check dimension of inputs
k1<-length(phosi);k2<-length(lambda);k3<-length(ppbt)
mytest<-abs(k1-k2)+abs(k2-k3)
if(mytest>0) stop("dimensions of phosi, lambda, and ppbt must match")
#check constraints on ppbt, phosi, and lambda
if(sum(ppbt<0))stop("ppbts must all be greater than or equal to zero")
#check that each ppbt is less than or equal to one
if(sum(ppbt>1))stop("ppbts must all be less than or equal to 1.0")
#check that all lambdas are between zero and 1.0
if(sum(lambda<0))stop("lambdas must all be greater than or equal to zero")
if(sum(lambda>1))stop("lambdas must all be less than or equal to one")
#check that all phosi are between zero and 1.0
if(sum(phosi<=0))stop("phosis must all be greater than zero")
if(sum(phosi>1))stop("phosis must all be less than or equal to one")
#check that the subsample size is less than sample size
if(n>Nsamp)stop(paste("n must be less than or equal to Nsamp=",Nsamp))
#get expected values of observations
Ex1<-Nsamp*sum(lambda*phosi)
Ex2<-Nsamp-Ex1
#check subsample. Note in the case of OPT==TRUE, n1 might be changed
if(!OPT){
if(n1<0)stop("n1 must be nonnegative")
if(n<n1)stop("n must be greater than or equal to n1")
if(n1>Ex1)stop(paste("n1 must not exceed the expected number of VM spawners in sample=",Ex1))
if(n1<(n-Ex2))stop(paste("n1 must be >= n minus the expected number of ~VM spawners in sample=",n-Ex2))}
n2<-n-n1
if(!OPT){
iii<-n1*lambda*ppbt+n2*(1-lambda)*ppbt==0}
#note in the case of OPT that n1 might be changed to allow estimation
if(OPT){
iii<-(lambda*ppbt+(1-lambda)*ppbt==0)|(n==0)}
onelambda2<-sum(mean(lambda[iii])==lambda[iii])==sum(iii)
onelambda2<-onelambda2&(mean(lambda[iii])>0)
if((!onelambda2)&sum(iii)){
stop("Expected tag recoveries must not be zero when marking fractions differ or are zero")}

return(NULL)
}

#Theoretical results using maximum likelihood theory
phos.pbt.estimates<-function(phosi=.1*c(1/20,1/20,9/20,9/20),
Nsamp=1000,n=200,n1=50,lambda=c(1,.95,.5,.5),ppbt=c(.95,.95,.95,.95),
suppress){
n2<-n-n1
nhatch<-length(lambda)
#get expected values of observations
Ex1<-Nsamp*sum(lambda*phosi)
Ex2<-Nsamp-Ex1
Ey<-n1*lambda*ppbt*phosi/sum(lambda*phosi)
Ez<-n2*(1-lambda)*ppbt*phosi/(1-sum(lambda*phosi))
Ey2<-Ex1*lambda*ppbt*phosi/sum(lambda*phosi)
Ez2<-Ex2*(1-lambda)*ppbt*phosi/(1-sum(lambda*phosi))
if(sum(lambda)==0){
Ey<-0
Ez<-n2*ppbt*phosi
Ey2<-0
Ez2<-Ex2*ppbt*phosi
}
phos<-sum(phosi)
#get estimates using pbte routines
res<-phos.pbte.estimates(x1=Ex1,x2=Ex2,n1=n1,n2=n2,y=Ey,z=Ez,lambda,ppbt,suppress=suppress)
res2<-phos.pbte.estimates(x1=Ex1,x2=Ex2,n1=Ex1,n2=Ex2,y=Ey2,z=Ez2,lambda,ppbt,suppress=suppress)
myres<-list(MONTE=FALSE,
NSIM=NA,
phosi=phosi,
Nsamp=Nsamp,
n=n,
n1=n1,
lambda=lambda,
ppbt=ppbt,
phos=phos,
Ex1=Ex1,
Ex2=Ex2,
SE_MIN.phos=res2$SE.phos,
CV_MIN.phos=res2$CV.phos,
SE.phos=res$SE.phos,
CV.phos=res$CV.phos,
BIAS.phos=NA)
return(myres)
}

#Monte Carlo estimates of standard error and relative bias
phos.pbt.estimates2<-function(NSIM=1000,phosi=.1*c(1/20,1/20,9/20,9/20),
Nsamp=1000,n=n,n1=50,lambda=c(1,.95,.5,.5),ppbt=c(.95,.95,.95,.95)){
n2<-n-n1
nhatch<-length(phosi)
#get expected values of observations
Ex1<-Nsamp*sum(lambda*phosi)
Ex2<-Nsamp-Ex1
Ey<-n1*lambda*ppbt*phosi/sum(lambda*phosi)
Ez<-n2*(1-lambda)*ppbt*phosi/(1-sum(lambda*phosi))
Ey2<-Ex1*lambda*ppbt*phosi/sum(lambda*phosi)
Ez2<-Ex2*(1-lambda)*ppbt*phosi/(1-sum(lambda*phosi))
if(sum(lambda)==0){
Ey<-0
Ez<-n2*ppbt*phosi
Ey2<-0
Ez2<-Ex2*ppbt*phosi
}
phos<-sum(phosi)
res<-phos.pbte.estimates2(NBOOT=NSIM,x1=Ex1,x2=Ex2,n1=n1,n2=n2,y=Ey,z=Ez,lambda,ppbt)
res2<-phos.pbte.estimates2(NBOOT=NSIM,x1=Ex1,x2=Ex2,n1=Ex1,n2=Ex2,y=Ey2,z=Ez2,lambda,ppbt)

myres<-list(MONTE=TRUE,
NSIM=NSIM,
phosi=phosi,
Nsamp=Nsamp,
n=n,
n1=n1,
lambda=lambda,
ppbt=ppbt,
phos=phos,
Ex1=Ex1,
Ex2=Ex2,
SE_MIN.phos=res2$SE.phos,
CV_MIN.phos=res2$CV.phos,
SE.phos=res$SE.phos,
CV.phos=res$CV.phos,
BIAS.phos=res$BIAS.phos)
return(myres)
}

#Maximum Likelihood Theory results
phos.pbte.estimates<-function(x1,x2,n1,n2,y,z,lambda,ppbt,suppress){
nhatch<-length(lambda)
Nsamp<-x1+x2

#An important case for combining cells occurs when the expected w=x+y=0
#in this case, constant lambdas over these cells saves the estimation.
#note that this also takes care of the case with a single lambda for all hatcheries
iii<-n1*lambda*ppbt+n2*(1-lambda)*ppbt==0
onelambda2<-sum(lambda[iii]==mean(lambda[iii]))==sum(iii)
if((sum(iii)>1)&onelambda2){
if(!suppress)warning("collapsing cells with expected tag recoveries of zero into single cell since lambda is constant")
lambda1.new<-mean(lambda[iii])
ppbt1.new<-0.0
lambda.new<-c(lambda1.new,lambda[!iii])
ppbt.new<-c(ppbt1.new,ppbt[!iii])
y.new<-c(0,y[!iii])
z.new<-c(0,z[!iii])
nhatch.new<-length(lambda.new)
res<-phos.pbte.estimates(x1=x1,x2=x2,n1=n1,n2=n2,y=y.new,z=z.new,lambda=lambda.new,ppbt=ppbt.new,suppress=suppress)
phosi<-rep(NA,nhatch)
SE.phosi<-rep(NA,nhatch)
phosi[!iii]<-res$phosi[2:nhatch.new]
SE.phosi[!iii]<-res$SE.phosi[2:nhatch.new]

myres<-list(BOOT=FALSE,
NBOOT=NA,
Nsamp=Nsamp,
x1=x1,
x2=x2,
n=n1+n2,
n1=n1,
n2=n2,
y=y,
z=z,
lambda=lambda,
ppbt=ppbt,
phosi=phosi,
SE.phosi=SE.phosi,
phos=res$phos,
SE.phos=res$SE.phos,
CV.phos=res$CV.phos,
BIAS.phos=NA)
return(myres)
}


#get initial estimate of phosi
phosi.init<-init(x1,x2,n1,n2,y,z,lambda,ppbt,suppress=suppress)

myres<-get.estimates2(phosi.init,x1,x2,n1,n2,y,z,lambda,ppbt,suppress=suppress)
phos<-myres$phos
phos.var<-myres$phos.var
phosi<-myres$phosi
SE.phosi<-sqrt(myres$phosi.var)


SE.phos<-sqrt(phos.var)
CV.phos<-SE.phos/phos
SE.phos<-as.numeric(SE.phos)
CV.phos<-as.numeric(CV.phos)

myres<-list(BOOT=FALSE,
NBOOT=NA,
Nsamp=Nsamp,
x1=x1,
x2=x2,
n=n1+n2,
n1=n1,
n2=n2,
y=y,
z=z,
lambda=lambda,
ppbt=ppbt,
phosi=phosi,
SE.phosi=SE.phosi,
phos=phos,
SE.phos=SE.phos,
CV.phos=CV.phos,
BIAS.phos=NA)
return(myres)
}

#Fisher Information Matrix (general case)
getI<-function(Nsamp,n1,n2,lambda,ppbt,phosi){
Ex1<-Nsamp*sum(lambda*phosi)
Ex2<-Nsamp-Ex1
theta1<-n1/Ex1
theta2<-n2/Ex2

v<-lambda
I<- v%*%t(v)*Nsamp*(1-theta1)/sum(v*phosi)
I<-I +v%*%t(v)*Nsamp*(1-theta2)/(1-sum(v*phosi))

v<-(1-ppbt)*lambda
I<-I+v%*%t(v)*Nsamp*theta1/sum(v*phosi)

v<-lambda*(1-ppbt)+ppbt
I<-I+v%*%t(v)*Nsamp*theta2/(1-sum(v*phosi))

#fix diagonal
mydiag<-diag(I)+Nsamp*ppbt*(theta1*lambda+theta2*(1-lambda))/phosi
diag(I)<-mydiag
return(I)
}

#Fisher Information matrix used when all VM releases PBT
getI2<-function(Nsamp,n1,n2,lambda,ppbt,phosi){
Ex1<-Nsamp*sum(lambda*phosi)
Ex2<-Nsamp-Ex1
theta1<-n1/Ex1
theta2<-n2/Ex2

v<-lambda
I<- v%*%t(v)*Nsamp*(1-theta1)/sum(v*phosi)
I<-I +v%*%t(v)*Nsamp*(1-theta2)/(1-sum(v*phosi))

v<-ppbt
I<-I+v%*%t(v)*Nsamp*theta2/(1-sum(v*phosi))

#fix diagonal
mydiag<-diag(I)+Nsamp*ppbt*(theta1*lambda+theta2*(1-lambda))/phosi
diag(I)<-mydiag
return(I)
}

#Fisher Information matrix (used when all lambdas are zero)
getI3<-function(Nsamp,n1,n2,lambda,ppbt,phosi){
Ex2<-Nsamp
theta2<-n2/Ex2
v<-ppbt
I<-v%*%t(v)*Nsamp*theta2/(1-sum(v*phosi))
#fix diagonal
mydiag<-diag(I)+Nsamp*ppbt*theta2/phosi
diag(I)<-mydiag
return(I)
}

#Bootstrap estimates of standard error and bias
#consider cases where some of the phosi are missing
phos.pbte.estimates2<-function(NBOOT,x1,x2,n1,n2,y,z,lambda,ppbt){
nhatch<-length(lambda)
Nsamp<-x1+x2
#get MLE using theoretical results
res<-phos.pbte.estimates(x1,x2,n1,n2,y,z,lambda,ppbt,suppress=FALSE)
phosi<-res$phosi
phosi.orig<-res$phosi
phos<-res$phos
phos.sim<-rep(NA,NBOOT)
phosi.sim<-matrix(NA,nrow=NBOOT,ncol=nhatch)
if(!is.na(phos)){
for(ii in 1:NBOOT){
iii<-is.na(phosi)
phosi[iii]<-rep(phos-sum(phosi,na.rm=T),sum(iii))/sum(iii)
mysim<-pbtsim1(phosi,Nsamp,n1,n2,res$lambda,res$ppbt)
my.n1<-n1
my.n2<-n2
if(mysim$x1<n1)my.n1<-mysim$x1
if(mysim$x2<n2)my.n2<-mysim$x2
res<-phos.pbte.estimates(x1=mysim$x1,x2=mysim$x2,n1=my.n1,n2=my.n2,
y=mysim$y,z=mysim$z,lambda=lambda,ppbt=ppbt,suppress=TRUE)
phos.sim[ii]<-res$phos
phosi.sim[ii,]<-res$phosi

}}

SE.phosi<-apply(phosi.sim,c(2),var,na.rm=T)
SE.phosi<-sqrt(SE.phosi)
SE.phos<-sqrt(var(phos.sim,na.rm=T))
CV.phos<-SE.phos/phos
mymean<-mean(phos.sim,na.rm=T)
BIAS.phos<-(mymean-phos)/phos

myres<-list(BOOT=TRUE,
NBOOT=NBOOT,
Nsamp=Nsamp,
x1=x1,
x2=x2,
n=n1+n2,
n1=n1,
n2=n2,
y=y,
z=z,
lambda=lambda,
ppbt=ppbt,
phosi=phosi.orig,
SE.phosi=SE.phosi,
phos=phos,
SE.phos=SE.phos,
CV.phos=CV.phos,
BIAS.phos=BIAS.phos)
return(myres)
}

#simulate data when phosi available
pbtsim1<-function(phosi,Nsamp,n1,n2,lambda,ppbt){
m<-length(phosi)
#first get binomial sample of fish marked and unmarked
P<-sum(phosi*lambda)
x1<-rbinom(n=1,size=Nsamp,prob=P)
x2<-Nsamp-x1
#next use multinomial distribution to simulate
#how many fish are pbt and how many are not pbt
py<-ppbt*phosi*lambda/P
pz<-ppbt*phosi*(1-lambda)/(1-P)
if(P>0){y<-rmultinom(n=1,size=min(n1,x1),prob=c(py,1-sum(py)))}
if(P==0){y<-matrix(0,ncol=1,nrow=length(phosi))}
if(P<1){z<-rmultinom(n=1,size=min(n2,x2),prob=c(pz,1-sum(pz)))}
if(P==1){z<-matrix(0,ncol=1,nrow=length(phosi))}

return(list(Nsamp=Nsamp,x1=x1,x2=x2,n1=n1,n2=n2,y=y[1:m,1],z=z[1:m,1]))
}

#Use R.A. Fisher's scoring algorithm to estimate phosi
#phosi represents the intial guess on input
get.estimates2<-function(phosi.init,x1,x2,n1,n2,y,z,lambda,ppbt,suppress){
Nsamp<-x1+x2
NTRIAL<-100
tolx<-1.e-5
nhatch<-length(phosi.init)
nhatch.orig<-nhatch
w<-y+z
nas<-rep(NA,nhatch)
phosi<-nas
phosi.var<-nas
phos<-NA
phos.var<-NA

#in a rare case init can return zeroes even though the true estimate is not zero
#which would defeat Fisher's scoring method
jjj<-phosi.init==0
phosi.init[jjj]<-phosi.init[jjj]+.00001
#find zeroes that occur when lambda=0 and n1*lambda*ppbt+n2*(1-lambda)*ppbt>0
jjj<-(n1*lambda*ppbt+n2*(1-lambda)*ppbt>0)&(lambda==0)&(w==0)
if(sum(jjj)>=1){
y<-y[!jjj]
z<-z[!jjj]
ppbt<-ppbt[!jjj]
lambda<-lambda[!jjj]
phosi.init<-phosi.init[!jjj]
nhatch<-length(y)}

#check for special cases
#check to see if all VM releases PBT
pbttest<-sum((ppbt-1)*lambda)==0
#check to see if all lambdas are zero (special case)
lambdatest<-sum(lambda==0)==nhatch

#get the right Fisher Information function
if(lambdatest){
my.getI<-getI3
dlike<-dlike3
}
else{

if(pbttest){
my.getI<-getI2
dlike<-dlike2
}
else{
my.getI<-getI
dlike<-dlike1
}
}

#use R.A. Fisher's scoring algorithm to find where the partial derivatives of the
#log-likelihood are zero. Use Fisher Information matrix in Newton's method to approx -Hessian
phosi<-phosi.init
errf<-0.0
alpha<-0.9
for(ii in 1:NTRIAL){
I<-my.getI(Nsamp,n1,n2,lambda,ppbt,phosi)
df<-dlike(phosi,x1,x2,n1,n2,y,z,lambda,ppbt)
size<-prod(dim(I))
if(size==0){
if(!suppress)warning("dimension of I is 0 x 0")
return(list(phosi=nas,phosi.var=nas,phos=NA,phos.var=NA))}
if(is.na(rcond(I))){
if(!suppress)warning("condition number of I is NA")
return(list(phosi=nas,phosi.var=nas,phos=NA,phos.var=NA))}
if(rcond(I)<1.e-15){
if(!suppress)warning("computationally singular information matrix")
return(list(phosi=nas,phosi.var=nas,phos=NA,phos.var=NA))}
delx<-solve(I)%*%df
phosi<-phosi+delx*(1-alpha)
phosi<-abs(phosi)
errx<-sum(abs(delx))/sum(abs(phosi))
alpha<-alpha*alpha
if(errx<=tolx)break
}
if(ii==NTRIAL){
if(!suppress)warning("maximum number of iterations was reached")
return(list(phosi=nas,phosi.var=nas,phos=NA,phos.var=NA))}
phos<-sum(phosi)
e<-rep(1,nhatch)
myvar<-solve(I)
phos.var<-t(e)%*%myvar%*%e
phosi.var<-diag(myvar)
full.phosi.var<-rep(0,nhatch.orig)
full.phosi<-rep(0,nhatch.orig)
#recall that jjj represents hatcheries with estimates of phosi=0
#reduction occurs only when sum(jjj>=1)
if(sum(jjj)>=1){
full.phosi.var[!jjj]<-phosi.var
full.phosi[!jjj]<-phosi
}
if(sum(jjj)<1){
full.phosi.var<-phosi.var
full.phosi<-phosi

}
return(list(phosi=full.phosi,phosi.var=full.phosi.var,phos=phos,phos.var=phos.var))
}

#get gradient of the log likelihood function
#phosi is the current best estimate of phosi
#used in most general case
dlike1<-function(phosi,x1,x2,n1,n2,y,z,lambda,ppbt){
#estimate phosi
Nsamp<-x1+x2
sum1<-sum(lambda*phosi)
sum2<-sum((1-ppbt)*lambda*phosi)
sum3<-sum(phosi*ppbt)
res<-lambda*(x1-n1)/sum1-lambda*(Nsamp-x1-n2)/(1-sum1)+y/phosi
res<-res+(n1-sum(y))*(1-ppbt)*lambda/sum2+z/phosi
res<-res-(n2-sum(z))*(lambda*(1-ppbt)+ppbt)/(1-sum2-sum3)
return(res)
}

#get gradient of the log likelihood function
#phosi is the current best estimate of phosi
#used in the case where sum(lambda*(1-ppbt))===0
dlike2<-function(phosi,x1,x2,n1,n2,y,z,lambda,ppbt){
#estimate phosi
Nsamp<-x1+x2
sum1<-sum(lambda*phosi)
sum3<-sum(phosi*ppbt)
res<-lambda*(x1-n1)/sum1-lambda*(Nsamp-x1-n2)/(1-sum1)+y/phosi
res<-res+z/phosi
res<-res-(n2-sum(z))*ppbt/(1-sum3)
return(res)
}

#get gradient of the log likelihood function
#phosi is the current best estimate of phosi
#used when lambda=0 at all hatcheries
dlike3<-function(phosi,x1,x2,n1,n2,y,z,lambda,ppbt){
#estimate phosi
Nsamp<-x1+x2
sum3<-sum(phosi*ppbt)
res<-z/phosi
res<-res-(n2-sum(z))*ppbt/(1-sum3)
return(res)
}

#get initial estimates of phosi
#by equating x1,y,and z to their expectations
#this yields 2*nhatch +1 equations with nhatch unknowns
#which is, in general, overdetermined.
init<-function(x1,x2,n1,n2,y,z,lambda,ppbt,suppress){
Nsamp<-x1+x2
nhatch<-length(lambda)
A1<-lambda
B1<-x1/Nsamp
A2<--n1*diag(lambda*ppbt,ncol=nhatch,nrow=nhatch)
LAMBDAMAT<-matrix(lambda,ncol=nhatch,nrow=nhatch)
LAMBDAMAT<-t(LAMBDAMAT)
YDIAG<-diag(y,nrow=nhatch,ncol=nhatch)
ZDIAG<-diag(z,nrow=nhatch,ncol=nhatch)
A2<-A2+YDIAG%*%LAMBDAMAT
B2<-rep(0,nhatch)
A3<-n2*diag((1-lambda)*ppbt,nrow=nhatch,ncol=nhatch)
A3<-A3+ZDIAG%*%LAMBDAMAT
B3<-z
A<-rbind(A1,A2,A3)
B<-c(B1,B2,B3)
if(rcond(t(A)%*%A)<1.e-15){
if(!suppress)warning("matrix t(A)%A in init() is computationally singular")}
phosi<-solve(t(A)%*%A)%*%t(A)%*%B
return(abs(phosi))
}
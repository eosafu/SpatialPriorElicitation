library(Matrix) 
library(magic)
library(truncnorm)
library(sparseMVN)
library(mvtnorm)
library(tidyverse)
library(sp)
rm(list = ls())
source("utilityfunctions.R")
source("sampUtilFunctions.R")
######################
# Spatial windows
WInd <- list(
  global = c(0,1),
  window1= c(0.2,0.8),
  window2= c(0.4,0.9)
)
# Type of distance
LongLatA = T
LongLatField = T
######################
Admeas = NULL

set.seed(1234)
m=30
m0=50

# Discount and Inflation factor
discount= 2
infl = 0.2

# Sample location and obtain distance matrix
loc <-  generateLoc(m,WInd$window1) 
distmat <- spDists(loc,loc,longlat = LongLatField)

# sample the field
## Calculate the covariance matrix
### Set hyperparameters
tau   =10.0
nu    = 1 # smoothness
kappa = 5
Qtheta=cMatern(distmat,nu, kappa)

## Sample Field
### Set the true field values
###########################
# Bimodal spatial effect
peak1 = c(0.3,0.3) # Window 1 peak
peak2 =c(0.8,0.8) # Window 2 peak
thetaTrue=meanField(loc ,peak1,peak2)
thetaTrue = thetaTrue-mean(thetaTrue)
###########################

### Sample Field
Field=t(rmvnorm(1,thetaTrue,(1/tau)*Qtheta))

#####Fixed effect############
BetaTrue=c(0.5,3,-1,-2)     #
BetaTrue=matrix(c(BetaTrue))#
#############################

### Covariates
A=diag(rep(1,m))
x1=rep(1,m)
x2=rnorm(m,0,1)
x3=rnorm(m,0,1)
x4=rnorm(m,0,1)
I=diag(rep(1,m))

X=t(as.matrix(cbind(x1,x2,x3,x4)))

### Predictor
eta= t(X)%*%BetaTrue+A%*%Field

# See later 236:260
### Response
#Y=t(rmvnorm(1,eta,(1/tau)* I))

##### Historical Data ######
# Sample location and obtain distance matrix
#### Location of historical data
loc0   <-  generateLoc(m0,WInd$window2)
distmat0 <- spDists(loc0,loc0,longlat = LongLatField)
Qtheta0 =cMatern(distmat0,nu, kappa)

## Sample Field
### Set the true field values

###########################
thetaTrue0 = meanField(loc0, peak1, peak2)#
thetaTrue0 = thetaTrue0-mean(thetaTrue0)
###########################
# Inflate
Field0=t(rmvnorm(1,thetaTrue0,infl*(1/tau)*Qtheta0))

# discount
##########################
BetaTrue0=discount*c(0.5,3,-1,-2)     #
BetaTrue0=matrix(c(BetaTrue0))#
#############################

### Covariates

A0=diag(rep(1,m0))
x01=rep(1,m0)
x02=rnorm(m0,0,1)
x03=rnorm(m0,0,1)
x04=rnorm(m0,0,1)
I=diag(rep(1,m0))

X0=t(as.matrix(cbind(x01,x02,x03,x04)))

### Predictor
eta0= t(X0)%*%BetaTrue0+A0%*%Field0

### Response (See later 236:260)
#Y0=t(rmvnorm(1,eta0,infl*(1/tau)* I))

##############################################
######### Scalable Commensurate Prior ########
##############################################

####### Field location ##
#Number of fields 
nfl0 =6
nfr0 =6
nfl  =6
nfr  =5
fieldcoord0 =  expand.grid(seq((min(WInd$window2)),(max(WInd$window2)),length=nfl0),
                          seq((min(WInd$window2)),(max(WInd$window2)),length=nfr0))
fieldcoord0 = as.matrix(fieldcoord0)
fieldcoord =  expand.grid(seq((min(WInd$window1)),(max(WInd$window1)),length=nfl),
                           seq((min(WInd$window1)),(max(WInd$window1)),length=nfr))
fieldcoord = as.matrix(fieldcoord)


nu= 1    # Fixed smoothness
kappa= 5 # Fixed range
sigmatheta =10
q1=nrow(fieldcoord0)
p1=4
q=nrow(fieldcoord)
p=4
qstar1=p1+q1
qstar=p+q
a=a1=2 # Shape for tau and tau1
b=b1=1 # Scale for tau and tau1
n0=10
n=10
b3=2   # Proposal parameter related to lambda beta (Global environ.)
b4=2   # proposal parameter related to lambda theta (Global environ.)
a11=2  # Shape for lambda beta
b11=1  # Scale for lambda beta
a22=2  # Shape for lambda theta
b22=1  # Scale for lambda theta
##
tau1=0.5
lambdabeta=1
lambdatheta=1
####################
## Matrices A and A0
psi = 0.2
A0 = exp(-spDists(loc0,fieldcoord0,longlat = LongLatA)/psi)
A0 =A0/matrix(rowSums(A0),nrow=nrow(A0),ncol=ncol(A0))

A1 = exp(-spDists(loc,fieldcoord,longlat = LongLatA)/psi)
A1 = A1/matrix(rowSums(A1),nrow=nrow(A1),ncol=ncol(A1))


A0= as(A0,"sparseMatrix")
A =  as(A1,"sparseMatrix")

####################
#(1) Compute covariance matrices for the historical field
Dmat <- as.matrix(spDists(fieldcoord0,fieldcoord0,longlat = LongLatField))
Q=(sigmatheta)*cMatern(Dmat,nu=nu,kappa = kappa)

NN=NearestN(Dmat,5)

AD=computeAD(Q,NN)
H=AD$A
VV=AD$D
Qtild <- (sigmatheta)* solve(diag(rep(1,nrow(Q)))-H)%*%VV%*%t(solve(diag(rep(1,nrow(Q)))-H))
InvQtild <- (1/sigmatheta)*t(diag(rep(1,nrow(Q)))-H)%*%solve(VV)%*%(diag(rep(1,nrow(Q)))-H)

# Matrix for Beta
Qbeta <- diag(rep(100,p))

# Matrix for xi derived from the separate beta and theta cov matrices.
SigMak0    <- bdiag(Qbeta,Qtild)
invSigmak0 <- solve(SigMak0)

#(2)
# Considering the projection steps
# Deriving Sigma00, from the current data including betas

Dmat0 <- as.matrix(spDists(fieldcoord,fieldcoord,longlat = LongLatField))
Q0=(sigmatheta)*cMatern(Dmat0,nu=nu,kappa = kappa)

# Similar neighboring scheme is used
NN0=NearestN(Dmat0,5)
AD=computeAD(Q0,NN0)
H=AD$A
VV=AD$D
Qtild0 <- (sigmatheta)* solve(diag(rep(1,nrow(Q0)))-H)%*%VV%*%t(solve(diag(rep(1,nrow(Q0)))-H))

#(3) Now, derive Sigma-01 i.e current vs historical
Dmat01 <- as.matrix(spDists(fieldcoord,fieldcoord0,longlat = LongLatField))
Q01=(sigmatheta)*cMatern(Dmat01,nu=nu,kappa = kappa)
Q01 = as(Q01,"sparseMatrix")

res       <- list(matrix(0,nrow(fieldcoord),p1), Q01)
res       <- do.call(cbind,res)
res       <- list(matrix(0,p,p1+nrow(fieldcoord0)),res)
Qtild01   <- do.call(rbind,res)

# hyper-parameters
hparam     = list()
hparam[[1]]=invSigmak0 # structured Cov. matrix for historical data
hparam[[2]]=Qtild0 # To comupute cov. for the current data
hparam[[3]]=Qtild01 # To compute cross-cov. of current and historical data
hparam[[4]]=sigmatheta # scale parameters for the covariance matrices
hparam[[5]]=a # Shape for current precision (tau) 
hparam[[6]]=b # rate parameter for (tau)
hparam[[7]]=a1# shape for historical precision (tau1)
hparam[[8]]=b1# rate for historical precision   (tau1)

#for(n in c(1,5,10,20, 50)){ # For multiple replications
  
  numbIT=1         # For varying the number of repeat to calculate the Out-sample mse and cpo
  
  inCPO  <- vector("numeric",numbIT)
  outMSE <- vector("numeric",numbIT)
  TIME   <- vector("numeric",numbIT)
  LAMTHETA <-vector("numeric",numbIT)
  
  for(index in 1:numbIT){
    n0=n
    Y=list()
    V=list()
    I=diag(rep(1,m))
    for (i in 1:n) {
      Y[[i]] <- t(rmvnorm(1,eta,(1/tau)* I))
      X=as(t(as.matrix(cbind(x1,x2,x3,x4))),"sparseMatrix") 
      X <- list(t(X),A)
      V[[i]] <- do.call(cbind,X)
    }
    Y <- do.call(rbind,Y)
    V <- do.call(rbind,V)
    
    Y0=list()
    V0=list()
    I=diag(rep(1,m0))
    for (i in 1:n) {
      Y0[[i]] <- t(rmvnorm(1,eta0,(1/tau)* I))
      X0=as(t(as.matrix(cbind(x01,x02,x03,x04))),"sparseMatrix") 
      X0 <- list(t(X0),A0)
      V0[[i]] <- do.call(cbind,X0)
    }
    Y0 <- do.call(rbind,Y0)
    V0 <- do.call(rbind,V0)
    
    ########################
    ########################
    
    # Matrix for xi.
    param      <- list()
    #Parameters and htperparameters
    param[[1]]=as(matrix(rnorm(qstar,0,1)),"sparseMatrix") # current spatial process
    param[[2]]=tau # current precision
    param[[3]]=as(matrix(rnorm(qstar1,0,1)),"sparseMatrix") # historical spatial process
    param[[4]]=tau1 # historical Precision
    param[[5]]=lambdabeta # commensurate for beta
    param[[6]]=lambdatheta # commensurate for theta
    param[[7]]=a11   # Shape for commensurate lambdabeta
    param[[8]]=b11   # rate for commensurate lambdabet
    param[[9]]=a22 # shape for commensurate lambdatheta
    param[[10]]=b22# rate for commensurate lambdatheta
    
    
    ##########################################
    # Estimation through joint scalable spatial prior
    # Metropolis hasting within Gibbs steps.  
    ##########################################
    #Space to save
    Num    <- 10000
    burn   <- floor(Num/3)
    XI     <- as(matrix(0,qstar,Num),"sparseMatrix")
    Tau    <- vector("numeric",Num)
    XI.prime     <- as(matrix(0,qstar1,Num),"sparseMatrix")
    Tau1    <- vector("numeric",Num)
    Lambeta <- vector("numeric",Num)
    Lamtheta<- vector("numeric",Num)
    Aa11     <- vector("numeric",Num)
    Bb11     <- vector("numeric",Num)
    Aa22     <- vector("numeric",Num)
    Bb22     <- vector("numeric",Num)
    
    ##########################################
    start_time <- Sys.time()
    for (i in 1:Num) {
      
      ################
      Aux <-  tryCatch({
        rex<- ScalableComm(param=param,hparam=hparam)
      },
      warning=function(war){
        rex <- rex
        return(rex)
      },
      error=function(err){
        rex <- "error"
        return(rex)
      },
      finally = {
        rex <- rex
      }
      
      )
      if(class(Aux)[1]!="character"){
        
        XI[,i]       <-param[[1]] <-  Aux[[1]]
        Tau[i]       <-param[[2]] <-  Aux[[2]] 
        XI.prime[,i] <-param[[3]] <-  Aux[[3]]
        Tau1[i]      <-param[[4]] <-  Aux[[4]] 
        Lambeta[i]   <-param[[5]] <-  Aux[[5]]
        Lamtheta[i]  <-param[[6]] <-  Aux[[6]]
        
        Aa11[i]       <-param[[7]] <-  Aux[[7]]
        Bb11[i]       <-param[[8]] <-  Aux[[8]]
        Aa22[i]       <-param[[9]] <-  Aux[[9]]
        Bb22[i]       <-param[[10]]<-  Aux[[10]]
      }else{
        
        XI[,i]       <-NA
        Tau[i]       <-NA
        XI.prime[,i] <-NA
        Tau1[i]      <-NA
        Lambeta[i]   <-NA
        Lamtheta[i]  <-NA
        
        Aa11[i]       <-NA
        Bb11[i]       <-NA
        Aa22[i]       <-NA
        Bb22[i]       <-NA
        cat("Failed")
        ###############
      }
       if(i%%500==0){
        cat(i,"\n")
      }
    }
    
    end_time <- Sys.time()
    ##### In-sample ########
    
    inCPO[index]  <-    mean(log(CPO(n,m,100,Y,V,XI[,-c(1:burn)],Tau[-c(1:burn)])))
    TIME[index]   <-    end_time-start_time
    #### Out-sample ########
    
    OUTmse <-   function(mout = 100){
      set.seed(189123)
      
      mm = mout
      locnew <- as.matrix(expand.grid(seq(0,1,length=20),seq(0,1,length=20)))#
      distmat= spDists(locnew,locnew,longlat = LongLatField)
      Qtheta =cMatern(distmat,nu, kappa)
      thetaTrue=meanField(locnew, peak1, peak2)#
      thetaTrue=thetaTrue-mean(thetaTrue)
      Field= t(rmvnorm(1,thetaTrue,(1/tau)*Qtheta))
      
      #####
      BetaTrue=c(0.5,3,-1,-2)     
      BetaTrue=matrix(c(BetaTrue))
      
      ### Covariates
      mm=nrow(locnew)
      A=diag(rep(1,mm))
      x1=rep(1,mm)
      x2=rnorm(mm,0,1)
      x3=rnorm(mm,0,1)
      x4=rnorm(mm,0,1)
      I=diag(rep(1,mm))
      X=t(as.matrix(cbind(x1,x2,x3,x4)))
      
      ### Predictor
      etaout= t(X)%*%BetaTrue+A%*%Field
      ### True Response
      Yout=t(rmvnorm(1,etaout,(1/tau)* I))
      
      #### Estimated
      #psi=0.5
      A1 = exp(-spDists(locnew,fieldcoord,longlat = LongLatA)/psi)
      A1 = A1/matrix(rowSums(A1),nrow=nrow(A1),ncol=ncol(A1))
      A =  as(A1,"sparseMatrix")
      Vout <- list(t(X),A)
      Vout=do.call(cbind,Vout)
      # out-of-sample MSE
      AvXI  <- rowMeans(XI[,-c(1:burn)])
      outMSE=  sqrt(mean((Vout%*%AvXI - Yout)^2))
      return(outMSE)
    }
    outMSE[index] <- OUTmse(100)
    
  }
  inCPO =mean(inCPO)
  outMSE=mean(outMSE)
  TIME <- mean(TIME)
  admeas = data.frame(replication=n,
                      AVinCPO = inCPO,
                      AVoutMSE =outMSE,
                      discount = discount,
                      infl = infl,
                      AVTime=TIME
  )
  Admeas=rbind(Admeas,admeas)
#}  

#  }
#}

#####
set.seed(89123445)
mout = 50
# Sample location and obtain distance matrix
#locnew <-  generateLoc(m,WInd$global)
locnew <-  as.matrix(expand.grid(seq(0,1,length=50),seq(0,1,length=50)))#
distmat=   spDists(locnew,locnew,longlat = T)

# sample the field
## Calculate the covariance matrix
### Set hyperparameters

tau=1.00
nu= 1 # smoothness
kappa= 5
Qtheta =cMatern(distmat,nu, kappa)

## Sample Field
### Set the true field values

###########################
thetaTrue=meanField(locnew, peak1, peak2)#
thetaTrue=thetaTrue-mean(thetaTrue)
###########################

### Sample Field
Field= t(rmvnorm(1,thetaTrue,(1/tau)*Qtheta))

#############################
BetaTrue=c(0.5,3,-1,-2)     #
BetaTrue=matrix(c(BetaTrue))#
#############################

### Covariates
m= nrow(locnew)
A=diag(rep(1,m))
x1=rep(1,m)
x2=rnorm(m,0,1)
x3=rnorm(m,0,1)
x4=rnorm(m,0,1)
I=diag(rep(1,m))

X=t(as.matrix(cbind(x1,x2,x3,x4)))

### Predictor
eta= t(X)%*%BetaTrue+A%*%Field
### Response

Yout=eta +sqrt((1/tau))*rnorm(nrow(eta)) #t(rmvnorm(1,eta,(1/tau)* I))

##########################
## Projected/ predicted###
##########################
A1 = exp(-spDists(locnew,fieldcoord,longlat = LongLatA)/psi)
A1 = A1/matrix(rowSums(A1),nrow=nrow(A1),ncol=ncol(A1))
A =  as(A1,"sparseMatrix")
V <- list(t(X),A)
V=do.call(cbind,V)
AvXI = rowMeans(XI[,-(1:burn)])
YpredNoborow =(V%*%AvXI)


ln = data.frame(as.data.frame(locnew),pred=as.matrix(YpredNoborow),true=as.matrix(Yout))
p1 =ln %>%ggplot()+
  geom_tile(aes(x=Var1,y=Var2,fill=true))+
  scale_fill_viridis_c(option = "H", direction = 1,name="", guide=guide_colorbar(
    barheight = unit(40, units = "mm"),
    barwidth = unit(1, units = "mm")
  ))+labs(x="",y="",title ="")
p2 =ln %>%ggplot()+
  geom_tile(aes(x=Var1,y=Var2,fill=pred))+
  scale_fill_viridis_c(option = "H", direction = 1,name="", guide=guide_colorbar(
    barheight = unit(40, units = "mm"),
    barwidth = unit(1, units = "mm")
  ))+labs(x="",y="",title ="")

ggarrange(p1,p2,
          nrow = 1, ncol = 2)

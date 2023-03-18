library(Matrix) 
library(magic)
library(truncnorm)
library(sparseMVN)
library(mvtnorm)
library(sp)
library(tidyverse)
rm(list = ls())
source("utilityfunctions.R")
source("sampUtilFunctions.R")
######################
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

discount=2
infl=0.2

# Sample location and obtain distance matrix
loc <-  generateLoc(m,WInd$window1) 
distmat <- spDists(loc,loc,longlat =LongLatField) #

## Calculate the covariance matrix
### Set spatial hyperparameters

tau   = 10.00
nu    = 1 # smoothness
kappa = 5
Qtheta=cMatern(distmat,nu, kappa)

##  Sample Field
### Set the true field values

########################
# Bimodal spatial effect
peak1 = c(0.3,0.3)
peak2 =c(0.8,0.8)
thetaTrue = meanField(loc ,peak1,peak2)
thetaTrue = thetaTrue-mean(thetaTrue)
###########################

### Sample Field
Field=t(rmvnorm(1,thetaTrue,(1/tau)*Qtheta))

####### Fixed effect#########
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
### Response

#Y=t(rmvnorm(1,eta,(1/tau)* I))

##### Historical Data ######
# Sample location and obtain distance matrix

#### Location of historical data

loc0     <-  generateLoc(m0,WInd$window2)

#### Covaraince matrix of historical data
distmat0<- spDists(loc0,loc0,longlat = LongLatField)
Qtheta0 =cMatern(distmat0,nu, kappa)

## Sample Field
### Set the true field values

###### Fixed effect #######
thetaTrue0 = meanField(loc0, peak1, peak2)#
thetaTrue0 = thetaTrue0-mean(thetaTrue0)
###########################
Field0=t(rmvnorm(1,thetaTrue0,infl*(1/tau)*Qtheta0))

######## Discount ###########
BetaTrue0=discount*c(0.5,3,-1,-2)#
BetaTrue0=matrix(c(BetaTrue0))   #
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
### Response

#Y0=t(rmvnorm(1,eta0,infl*(1/tau)* I))

#########################################
####Commensurate Prior Borrowing ########
#########################################

####### Field location ########
#Number of fields 
nf  <- min(m,m0)
fieldcoord =  expand.grid(seq((min(WInd$global)+0.1),(max(WInd$global)-0.1),length=nf/5),
                          seq((min(WInd$global)+0.1),(max(WInd$global)-0.1),length=nf/6))
fieldcoord = as.matrix(fieldcoord)
q= nrow(fieldcoord)

nu         =1 # Fixed smoothness
kappa      =5 # Fixed range
sigmatheta =10
p          =4
qstar      =p+q
a          =2 # For shape for tau
b          =1 # scale for tau
a1=2  # Shape for lambda beta
b1=1  # Scale for lambda beta
a2=2  # Shape for lambda theta
b2=1  # Scale for lambda theta
n0         =10
n          =10

## Matrices A and A0
psi = 5
A0 = exp(-spDists(loc0,fieldcoord,longlat = LongLatA)/psi)
A0 =A0/matrix(rowSums(A0),nrow=nrow(A0),ncol=ncol(A0))
A0= as(A0,"sparseMatrix")

A1 = exp(-spDists(loc,fieldcoord,longlat = LongLatA)/psi)
A1 = A1/matrix(rowSums(A1),nrow=nrow(A1),ncol=ncol(A1))
A1 = as(A1,"sparseMatrix")

A00=as(A0,"sparseMatrix")
A1=as(A1,"sparseMatrix")

# Field covariance matrix

Dmat <- as.matrix(spDists(fieldcoord,fieldcoord,longlat = LongLatField))
Q=(sigmatheta)*cMatern(Dmat,nu=nu,kappa = kappa)
NN=NearestN(Dmat,5)
AD=computeAD(Q,NN)
H=AD$A
VV=AD$D
Qtild <- (sigmatheta)* solve(diag(rep(1,nrow(Q) ))-H)%*%VV%*%t(solve(diag(rep(1,nrow(Q)  ))-H))
InvQtild <- (1/sigmatheta)*t(diag(rep(1,nrow(Q) ))-H)%*%solve(VV)%*%(diag(rep(1,nrow(Q) ))-H)

# Matrix for Beta
Qbeta <- diag(rep(100,p))
# Matrix for xi.
SigMak     <-bdiag(Qbeta,Qtild)
invSigmak  <-solve(SigMak)

hparam      <- list()
hparam[[1]]=invSigmak
hparam[[2]]=sigmatheta
hparam[[3]]=a1
hparam[[4]]=b1
hparam[[5]]=a2
hparam[[6]]=b2

#for(n in c(1,5,10,20,50)){ # For multiple replications
  
  numbIT=1                  # For varying the number of repeat to calculate the Out-sample mse and cpo
  
  inCPO  <- vector("numeric",numbIT)
  outMSE <- vector("numeric",numbIT)
  TIME   <- vector("numeric",numbIT)
  
  for(index in 1:numbIT){
    n0=n
    Y=list()
    V=list()
    Y0=list()
    V0=list()
    
    IY0=diag(rep(1,m0))
    for (i in 1:n0) {
      Y0[[i]] <- t(rmvnorm(1,eta0,(1/tau)* IY0))
      X0      <- as(t(as.matrix(cbind(x01,x02,x03,x04))),"sparseMatrix")
      X0      <- list(t(X0),A00)
      V0[[i]] <- do.call(cbind,X0)
    }
    
    Y0=do.call(rbind,Y0)
    V0=do.call(rbind,V0)
    IY=diag(rep(1,m))
    for (i in 1:n) {
      Y[[i]] <- t(rmvnorm(1,eta,(1/tau)* IY))
      X=as(t(as.matrix(cbind(x1,x2,x3,x4))),"sparseMatrix")
      X <- list(t(X),A1)
      V[[i]] <- do.call(cbind,X)
    }
    Y=do.call(rbind,Y)
    V=do.call(rbind,V)
    ################
    # Matrix for xi.
    param      <- list()
    param[[1]]=as(matrix(rnorm(qstar,0,1)),"sparseMatrix")
    param[[2]]=0.5
    param[[3]]=1
    param[[4]]=1
    param[[5]]=a1
    param[[6]]=b1
    param[[7]]=a2
    param[[8]]=b2
    ##########################################
    # Estimation through joint Commensurate prior
    # Metropolis hasting within Gibbs steps.  
    ##########################################
    #Space to save
    Num    <- 10000
    burn   <- floor(Num/3)
    XI     <- as(matrix(0,qstar,Num),"sparseMatrix")
    Tau    <- vector("numeric",Num)
    Lambeta <- vector("numeric",Num)
    Lamtheta<- vector("numeric",Num)
    Aa1     <- vector("numeric",Num)
    Bb1     <- vector("numeric",Num)
    Aa2     <- vector("numeric",Num)
    Bb2     <- vector("numeric",Num)
    ##########################################
    start_time <- Sys.time()
    for (i in 1:Num) {
      ################
      Aux <-  tryCatch({
        rex<- commensurate(param=param,hparam=hparam)
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
        # Aux <- commensurate(param=param,hparam=hparam)
        XI[,i]       <-param[[1]] <-  Aux[[1]]
        Tau[i]       <-param[[2]] <-  Aux[[2]] 
        Lambeta[i]   <-param[[3]] <-  Aux[[3]]
        Lamtheta[i]  <-param[[4]] <-  Aux[[4]]
        
        Aa1[i]       <-param[[5]] <-  Aux[[5]]
        Bb1[i]       <-param[[6]] <-  Aux[[6]]
        Aa2[i]       <-param[[7]] <-  Aux[[7]]
        Bb2[i]       <-param[[8]] <-  Aux[[8]]
      }else{
        
        XI[,i]       <-NA
        Tau[i]       <-NA
        Lambeta[i]   <-NA
        Lamtheta[i]  <-NA
        
        Aa1[i]       <-NA
        Bb1[i]       <-NA
        Aa2[i]       <-NA
        Bb2[i]       <-NA
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
      tau=1.00
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
      outMSE=  sqrt(mean((Vout%*%AvXI - Yout)^2) ) #mean(log(CPO(1,mm,100,Yout,Vout,XI,Tau)))
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

set.seed(89123445)
mout = 50
m=mout
# Sample new location and obtain distance matrix
locnew <-  as.matrix(expand.grid(seq(0,1,length=100),seq(0,1,length=100)))#
distmat=   spDists(locnew,locnew,longlat = T)

## Calculate the covariance matrix
### Set hyperparameters

Qtheta =cMatern(distmat,nu, kappa)

## Sample Field
### Set the true field values

###########################
peak1 = c(0.3,0.3)
peak2 = c(0.8,0.8)
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

Yout=eta +sqrt((1/tau))*rnorm(nrow(eta)) 

psi=0.5
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

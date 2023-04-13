PostcondFullBorrowin <- function(param,hparam){
  # Function to sample from the posterior distribution under full borrowing
  # Input:
  # Parameter list (param)
  # hyper-parameter list (hparam)
  # Output: 
  # Updated parameter list
  # 
  
  xi    =param[[1]]
  tau   =param[[2]]
  
  invSigmak =hparam[[1]]
  sigmatheta=hparam[[2]]
  a         =hparam[[3]]
  b         =hparam[[4]]
  
  Sigmaxi=(1/tau)*solve(t(V)%*%V+invSigmak)
  CH     = Cholesky(symmpart(Sigmaxi)) 
  muxi   =(tau)*Sigmaxi%*%t(V)%*%Y
  xi     = rmvn.sparse(1,muxi,CH,prec=FALSE)
  xi     = t(xi)
  atau=(m*n+m0*n0+qstar+2*a)/2
  btau=0.5*(t(Y-V%*%xi)%*%(Y-V%*%xi)+ t(xi)%*%invSigmak%*%xi+2*b)
  tau=rgamma(1,shape = drop(atau),rate=drop(btau))
  
  param[[1]]=xi
  param[[2]]=tau
  return(param)
}


PostcondNoBorrowin <- function(param,hparam){
  # Function to sample from the posterior distribution under no borrowing
  # Input:
  # Parameter list (param)
  # hyper-parameter list (hparam)
  # Output: 
  # Updated parameter list
  # 
  
  xi    =param[[1]]
  tau   =param[[2]]
  
  invSigmak =hparam[[1]]
  sigmatheta=hparam[[2]]
  a         =hparam[[3]]
  b         =hparam[[4]]
  
  Sigmaxi=(1/tau)*solve(t(V)%*%V+invSigmak)
  CH     = Cholesky(Sigmaxi) 
  muxi   =(tau)*Sigmaxi%*%t(V)%*%Y
  xi     = rmvn.sparse(1,muxi,CH,prec=FALSE)
  xi     = t(xi)
  atau=(m*n+qstar+2*a)/2
  btau=0.5*(t(Y-V%*%xi)%*%(Y-V%*%xi)+ t(xi)%*%invSigmak%*%xi+2*b)
  tau=rgamma(1,shape = drop(atau),rate=drop(btau))
  
  param[[1]]=xi
  param[[2]]=tau
  return(param)
}


PostSamplepower <- function(param,hparam){
  # function to sample from the posterior distribution under power prior
  # Input:
  # Parameter list (param)
  # hyper-parameter list (hparam)
  # Output: 
  # Updated parameter list
  # 
  
  xi    =param[[1]]
  tau   =param[[2]]
  omega0=param[[3]]
  mub   =param[[4]]
  sigmab=param[[5]]
  
  invSigmak    =hparam[[1]]
  sigmatheta   =hparam[[2]]
  a            =hparam[[3]]
  b            =hparam[[4]]
  
  
  invD <- omega0*t(V0)%*%V0+invSigmak
  D    =solve(invD)
  B    <- omega0*t(V0)%*%Y0
  k    <- omega0*(t(Y0)%*%Y0)-t(B)%*%D%*%B+2*b
  v    <- 2*a+m0*n0*omega0
  
  # Xi
  Sigmaxi <- (1/tau)*solve( t(V)%*%V+ invD)
  CH      <- Cholesky(Sigmaxi) 
  muxi    <- tau*Sigmaxi%*%(t(V)%*%Y+B)
  xi      <- rmvn.sparse(1,muxi,CH,prec=FALSE)
  xi      <- t(xi)
  # tau
  atau <- (m*n+m0*n0*omega0+qstar+2*a)/2
  auxbtau <- t(xi-D%*%B)%*%invD%*%(xi-D%*%B)+k
  btau <- 0.5*(t(Y-V%*%xi)%*%(Y-V%*%xi)+auxbtau)
  tau  <- rgamma(1,shape=drop(atau),rate=drop(btau))
  
  # Metropolis Hasting
  # omega
  Logcondomega0<- function(omega0){
    auxnumwo   <- (0.5*v)*log(2*k)+(0.5*m0*n0*omega0)*log(tau)
    auxdenomwo <- lgamma(0.5*v)+0.5*log(det(D)) ### Use permutation mat here.
    shpw01     <- mub*(mub*(1-mub)/sigmab-1)
    shpw02     <- shpw01*(1-mub)/mub
    
    piw0 <- auxnumwo-auxdenomwo -0.5*tau*auxbtau+dbeta(omega0,10,10,log = TRUE)#dbeta(omega0,shpw01,shpw02,log = TRUE)
    return(piw0)
  }
  # draw wo from proposal dist.
  new.omega0 = rtruncnorm(1,0,1,omega0,0.3)
  
  aux.ratio <- (Logcondomega0(new.omega0)+ log(dtruncnorm(omega0,0,1,new.omega0,0.3)) )-
    (Logcondomega0(omega0)+ log(dtruncnorm(new.omega0,0,1,omega0,0.3)) )
  aux.ratio <- min(0,drop(aux.ratio))
  U <- runif(1)
  if(U<exp(aux.ratio)){
    omega0 =new.omega0
  }else{
    omega0 = omega0
  }
  
  # Update mub and sigmab
  
  ##############################
  ##############################
  
  condsigmab = function(sigmab){
    
    shpw01     <- mub*(mub*(1-mub)/sigmab-1)
    shpw02     <- shpw01*(1-mub)/mub
    
    dbeta(omega0,shpw01,shpw02,log = TRUE)#+log(1/(mub*(1-mub)))
  }
  
  mm=mub*(1-mub)
  new.sigmab = rtruncnorm(1,0,mm,0.5*mm,10) 
  # compute quantity for acceptance
  
  aux.ratio = (condsigmab(new.sigmab)+log(dtruncnorm(sigmab,0,mm,0.5*mm,10)))-
    
    (condsigmab(sigmab)+log(dtruncnorm(new.sigmab,0,mm,0.5*mm,10)))
  
  aux.ratio = min(0,aux.ratio)
  # test acceptance condition
  U = runif(1)
  
  # omega0 = ifelse(U < exp(aux.accept),new.omega0,omega0)
  if (U<exp(aux.ratio)){
    sigmab=new.sigmab
  }else{
    sigmab=sigmab
  }
  
  
  # posterior for muob
  
  condmub = function(mub){
    shpw01     <- mub*(mub*(1-mub)/sigmab-1)
    shpw02     <- shpw01*(1-mub)/mub
    
    dbeta(omega0,shpw01,shpw02,log = TRUE)#+log(1)
  }
  
  TEST=FALSE
  while (!TEST) {
    # proposal dist u(0.2,0.8)
    new.mub = runif(1,0.2,0.8)
    if(new.mub*(1-new.mub)>sigmab) TEST=TRUE
  }
  # compute quantity for acceptance
  aux.ratio = (condmub(new.mub)+dunif(mub,0.2,0.8,log = TRUE))-
    
    (condmub(mub)+dunif(new.mub,0.2,0.8,log = TRUE))
  aux.ratio = min(0,aux.ratio)
  # test acceptance condition
  U = runif(1)
  
  # omega0 = ifelse(U < exp(aux.accept),new.omega0,omega0)
  if (U<exp(aux.ratio)){
    mub=new.mub
  }else{
    mub=mub
  }
  ##################################
  
  param[[1]]= xi 
  param[[2]]=tau 
  param[[3]]= omega0
  param[[4]]= mub
  param[[5]]= sigmab
  return(param)
  
}

################
commensurate <- function(param,hparam){
  # Function to sample from the posterior distribution under commensurate borrowing
  # Input:
  # Parameter list (param)
  # hyper-parameter list (hparam)
  # Output: 
  # Updated parameter list
  # 
  
  xi         =param[[1]]
  tau        =param[[2]]
  lambdabeta =param[[3]]
  lambdatheta=param[[4]]
  a1         =param[[5]]
  b1         =param[[6]]
  a2         =param[[7]]
  b2         =param[[8]]
  
  invSigmak    =hparam[[1]]
  sigmatheta   =hparam[[2]]
  #a1           =hparam[[3]]
  #b1           =hparam[[4]]
  #a2           =hparam[[5]]
  #b2           =hparam[[6]]
  
  
  Itheta <- diag(rep(1,q))
  Ibeta <-  diag(rep(1,p))
  #Lambda <- sum(diag(t(V0)%*%V0))*bdiag(lambdabeta*Ibeta,lambdatheta*Itheta)
  Lambda <-  bdiag(lambdabeta*Ibeta,lambdatheta*Itheta)
  
  invD <- t(V0)%*%V0+Lambda
  
  D    <- solve(invD)
  
  invE=(invSigmak+Lambda-t(Lambda)%*%D%*%Lambda)
  E <- solve(invE)
  Faux <- t(Lambda)%*%D%*%t(V0)%*%Y0 
  
  Sigmaxi <- (1/tau)*solve(t(V)%*%V+invE)
  CH      <- Cholesky(Sigmaxi)
  muxi    <- tau*Sigmaxi %*% (t(V)%*%Y+Faux)
  xi     = rmvn.sparse(1,muxi,CH,prec=FALSE)
  xi     = t(xi)
  ## update tau
  atau <- (m*n+m0*n0+qstar+2*a)/2
  btau <- 0.5*(  t(Y-V%*%xi)%*%(Y-V%*%xi)+t(Y0)%*%Y0+
                   t(xi)%*%invE%*%xi-2*t(xi)%*%Faux-
                   t(Y0)%*%V0%*%D%*%t(V0)%*%Y0+2*b
  )
  tau=rgamma(1,shape = drop(atau),rate=drop(btau))
  
  ## update lambdas
  
  
  bivdenseLambda.beta <- function(lambdabeta){
     #Lambda <- diag(diag(t(V0)%*%V0)*0.01+0.001)%*%bdiag(lambdabeta*Ibeta,lambdatheta*Itheta)
    #Lambda <- sum(diag(t(V0)%*%V0))*bdiag(lambdabeta*Ibeta,lambdatheta*Itheta)
    Lambda <- bdiag(lambdabeta*Ibeta,lambdatheta*Itheta)
    invD <- t(V0)%*%V0+Lambda
    D    <- solve(invD)
    #invE=(invSigmak+Lambda-t(Lambda)%*%D%*%Lambda) # remove invSigmak to simplify
    invE=(Lambda-t(Lambda)%*%D%*%Lambda)
    #E <- solve(invE)
    Faux <- t(Lambda)%*%D%*%t(V0)%*%Y0 
    #flat prior for Lambdas
    aux= 0.5*log(det(D))-
      0.5*log(det(solve(Lambda)))-0.5*tau*(
        t(xi)%*%invE%*%xi-2*t(xi)%*%Faux-t(Y0)%*%V0%*%D%*%t(V0)%*%Y0
      )+dgamma(lambdabeta,shape=a1,rate=b1,log=TRUE)
    return(aux)
  }
  bivdenseLambda.theta <- function(lambdatheta){
    
    #Lambda <- (diag(diag(t(V0)%*%V0))+diag(rep(1,m+m0+p)))%*%bdiag(lambdabeta*Ibeta,lambdatheta*Itheta)
    Lambda <- bdiag(lambdabeta*Ibeta,lambdatheta*Itheta)
    invD   <- t(V0)%*%V0+Lambda
    D    <- solve(invD)
    #invE=(invSigmak+Lambda-t(Lambda)%*%D%*%Lambda) # remove invSigmak to simplify
    invE=(Lambda-t(Lambda)%*%D%*%Lambda)
    #E <- solve(invE)
    Faux <- t(Lambda)%*%D%*%t(V0)%*%Y0 
    #flat prior for Lambdas
    aux= 0.5*log(det(D))-
      0.5*log(det(solve(Lambda)))-0.5*tau*(
        t(xi)%*%invE%*%xi-2*t(xi)%*%Faux-t(Y0)%*%V0%*%D%*%t(V0)%*%Y0
      )+ dgamma(lambdatheta,shape=a2,rate=b2,log=TRUE)
    return(aux)
  }
  #fix the rate and swapped shape...to keep both borrowing close.
  new.lambdabeta <- rgamma(1,shape=5,rate=1)
  new.lambdatheta<- rgamma(1,shape=5,rate=1)
  
  aux.ratio = (bivdenseLambda.beta(new.lambdabeta)+
                 dgamma(lambdabeta,shape=5,rate=1,log = T))-
    #dgamma(lambdabeta,shape=b2*new.lambdatheta,rate=b2))-
    
    ( bivdenseLambda.beta(lambdabeta)+
        dgamma(new.lambdabeta,shape=5,rate=1,log = T))
  #  dgamma(new.lambdatheta,shape=b2*lambdabeta,rate=b2) )
  
  aux.ratio <- min(0,drop(aux.ratio))
  U <- runif(1)
  if(U<exp(aux.ratio)){
    lambdabeta =new.lambdabeta
    #lambdatheta=new.lambdatheta
  }else{
    lambdabeta =lambdabeta
    #lambdatheta=lambdatheta
  }
  
  #######
  aux.ratio = (bivdenseLambda.theta(new.lambdatheta)+
                 dgamma(lambdatheta,shape=5,rate=1,log=T))-
    #dgamma(lambdabeta,shape=b2*new.lambdatheta,rate=b2))-
    
    (bivdenseLambda.theta(lambdatheta)+
       dgamma(new.lambdatheta,shape=5,rate=1,log = T))
  #  dgamma(new.lambdatheta,shape=b2*lambdabeta,rate=b2) )
  
  aux.ratio <- min(0,drop(aux.ratio))
  U <- runif(1)
  if(U<exp(aux.ratio)){
    #lambdabeta =new.lambdabeta
    lambdatheta=new.lambdatheta
  }else{
    #lambdabeta =lambdabeta
    lambdatheta=lambdatheta
  }
  
  
  # Update hyper-param a1,b1,a2,b2
  
  new.a1 <- runif(1,0.5,10)
  aux.ratio <- (dgamma(lambdabeta,shape=new.a1,rate=b1 ,log=TRUE)+
               dunif(a1,0.5,10,log = TRUE))-
            (dgamma(lambdabeta,shape=a1,rate=b1 ,log=TRUE)+
             dunif(new.a1,0.5,10,log = TRUE))
  
  aux.ratio <- min(0,drop(aux.ratio))
  U <- runif(1)
  if(U<exp(aux.ratio)){
    a1=new.a1
  }else{
    a1=a1
  }
  
  new.b1 <- runif(1,0.5,10)
  aux.ratio <- (dgamma(lambdabeta,shape=a1,rate=new.b1 ,log=TRUE)+
                  dunif(b1,0.5,10,log = TRUE))-
    (dgamma(lambdabeta,shape=a1,rate=b1 ,log=TRUE)+
       dunif(new.b1,0.5,10,log = TRUE))
  
  aux.ratio <- min(0,drop(aux.ratio))
  U <- runif(1)
  if(U<exp(aux.ratio)){
    b1=new.b1
  }else{
    b1=b1
  }
  
  ######
  
  new.a2 <- runif(1,0.5,10)
  aux.ratio <- (dgamma(lambdatheta,shape=new.a2,rate=b2 ,log=TRUE)+
                  dunif(a2,0.5,10,log = TRUE))-
    (dgamma(lambdatheta,shape=a2,rate=b2 ,log=TRUE)+
       dunif(new.a2,0.5,10,log = TRUE))
  
  aux.ratio <- min(0,drop(aux.ratio))
  U <- runif(1)
  if(U<exp(aux.ratio)){
    a2=new.a2
  }else{
    a2=a2
  }
  
  new.b2 <- runif(1,0.5,10)
  aux.ratio <- (dgamma(lambdatheta,shape=a2,rate=new.b2 ,log=TRUE)+
                  dunif(b2,0.5,10,log = TRUE))-
    (dgamma(lambdatheta,shape=a2,rate=b2 ,log=TRUE)+
       dunif(new.b2,0.5,10,log = TRUE))
  
  aux.ratio <- min(0,drop(aux.ratio))
  U <- runif(1)
  if(U<exp(aux.ratio)){
    b2=new.b2
  }else{
    b2=b2
  }
  
  param[[1]]=xi
  param[[2]]=tau
  param[[3]]=lambdabeta
  param[[4]]=lambdatheta
  param[[5]]=a1
  param[[6]]=b1
  param[[7]]=a2
  param[[8]]=b2
  
  return(param)
  
}



ScalableComm <- function(param,hparam){
  # Function to sample from the posterior distribution under scalable commensurate borrowing
  # Input:
  # Parameter list (param)
  # hyper-parameter list (hparam)
  # Output: 
  # Updated parameter list
  # 
  #Parameters and htperparameters
  xi        =param[[1]] # current spatial process
  tau        =param[[2]] # current precision
  xi.prime   =param[[3]] # historical spatial process
  tau1       =param[[4]] # historical Precision
  lambdabeta =param[[5]] # commensurate for beta
  lambdatheta=param[[6]] # commensurate for theta
  a11        =param[[7]] # Shape for commensurate lambdabeta
  b11        =param[[8]] # rate for commensurate lambdabet
  a22        =param[[9]] # shape for commensurate lambdatheta
  b22        =param[[10]]# rate for commensurate lambdatheta
  
  # hyper-parameters
  invSigmak0 =hparam[[1]] # structured Cov. matrix for historical data
  Qtild0     =hparam[[2]] # To comupute cov. for the current data
  Qtild01    =hparam[[3]] # To compute cross-cov. of current and historical data
  sigmatheta =hparam[[4]] # scale parameters for the covariance matrices
  a          =hparam[[5]] # Shape for current precision (tau) 
  b          =hparam[[6]] # rate parameter for (tau)
  a1         =hparam[[7]] # shape for historical precision (tau1)
  b1         =hparam[[8]] # rate for historical precision   (tau1)
  
  
  
  Sigmaxi.prime <- (1/tau1)*solve(t(V0)%*%V0+invSigmak0)
  CH            <- Cholesky(Sigmaxi.prime)
  muxi.prime    <- tau1*Sigmaxi.prime%*%t(V0)%*%Y0
  xi.prime      <- rmvn.sparse(1,muxi.prime,CH,prec=FALSE)
  xi.prime      <- t(xi.prime)
  atau1 <- (m0*n0+qstar1+2*a1)/2
  btau1 <- 0.5*( t(Y0-V0%*%xi.prime)%*%(Y0-V0%*%xi.prime)+
                   t(xi.prime)%*%invSigmak0%*%xi.prime+2*b1)
  
  tau1 <- rgamma(1,shape=drop(atau1),rate = drop(btau1))
  # Now use tau1 to update the covariance matrix to have Sigma-prime
  Sigma.prime <- (1/tau)*SigMak0            # Sigma prime 
  
  ### Considering the projection steps
  # Deriving Sigma00, from the current data including betas
  Qbeta0 <- diag(rep(100,p))
  SigMa00    <- (1/tau1)*bdiag(Qbeta0,Qtild0) # (2) Sigma00
  
  ### Now, derive Sigma-01 i.e current vs historical
  SigMa01    <- (1/tau1)*Qtild01             # (3)
  
  ## Use these cov. matrices to find c and D
  invSigma.prime <- solve(Sigma.prime)
  c.t <- SigMa01%*%invSigma.prime
  c   <- t(c.t)
  D <- SigMa00- SigMa01%*%invSigma.prime%*%t(SigMa01)        # Sigma prime
  
  
  # Compute Lambda
  Itheta <- diag(rep(1,nrow(fieldcoord)))
  Ibeta <-  diag(rep(1,p))
  Lambda <- bdiag(lambdabeta*Ibeta,lambdatheta*Itheta)
  
  # Compute E and B
  invD <- solve(D)
  invE <- (Lambda+invD)
  E    <- solve(invE)
  
  # Compute K
  invK <- (c%*%invD%*%t(c)+invSigma.prime-c%*%invD%*%E%*%invD%*%t(c)) #+0.00001*diag(rep(1,m0+p1))
  K    <- solve(invK)         #ginv(as.matrix(invK), tol = sqrt(.Machine$double.eps))
  # Compute M
  invM <- (Lambda-Lambda%*%E%*%Lambda-Lambda%*%E%*%invD%*%t(c)%*%K%*%c%*%invD%*%E%*%Lambda)
  
  M=solve(invM)
  W <- Lambda %*%E%*%invD%*%t(c)%*%K%*%invSigma.prime%*%muxi.prime 
  
  ### Update current parameters
  SigMa00k    <- bdiag(Qbeta0,Qtild0) 
  invSigMa00k <- solve(SigMa00k)
  
  Sigmaxi <- solve(invM+tau*invSigMa00k+tau*t(V)%*%V)
  Sigmaxi=as(Sigmaxi,"sparseMatrix")
  CH      <- Cholesky(forceSymmetric(Sigmaxi))
  muxi    <- Sigmaxi%*%(tau*t(V)%*%Y+W)
  xi     = rmvn.sparse(1,muxi,CH,prec=FALSE)
  xi     = t(xi)
  
  # Update Gamma 
  
  atau <- (m*n+qstar+2*a)/2
  btau <- 0.5*( t(Y-V%*%xi)%*%(Y-V%*%xi)+t(xi)%*%invSigMa00k%*%xi+
                  2*b)
  tau=rgamma(1,shape = drop(atau),rate=drop(btau))
  
  # Update Lambdas
  
  
  bivdenseLambda<- function(lambdabeta,lambdatheta){
    Lambda <- bdiag(lambdabeta*Ibeta,lambdatheta*Itheta)
    invE <- (Lambda+invD)
    E    <- solve(invE)
    invK <- (c%*%invD%*%t(c)+invSigma.prime-c%*%invD%*%E%*%invD%*%t(c))
    K    <- solve(invK)       
    invM <- (Lambda-Lambda%*%E%*%Lambda-Lambda%*%E%*%invD%*%t(c)%*%K%*%c%*%invD%*%E%*%Lambda)
    W <- Lambda %*%E%*%invD%*%t(c)%*%K%*%invSigma.prime%*%muxi.prime 
    ##########
    #indp. gamma prior for Lambdas
    aux= 0.5*log(det(E))+0.5*log(det(K))-
      0.5*log(det(solve(Lambda)))-0.5*(
        t(xi)%*%invM%*%xi-2*t(xi)%*%W-t(muxi.prime)%*%invSigma.prime%*%(K%*%invSigma.prime-diag(rep(1,p1+nrow(fieldcoord0))))%*%muxi.prime 
      )+dgamma(lambdabeta,shape=a11,rate=b11,log=TRUE)+
      dgamma(lambdatheta,shape=a22,rate=b22,log=TRUE)
    return(aux)
  }
  #fix the rate and swapped shape...to keep both borrowing close.
  new.lambdabeta <- rgamma(1,shape=b3*lambdabeta,rate=b3)
  new.lambdatheta<- rgamma(1,shape=b4*lambdatheta,rate=b4)
  
  aux.ratio = (bivdenseLambda(new.lambdabeta,new.lambdatheta)+
                 dgamma(lambdabeta,shape=b3*new.lambdabeta,rate=b3,log = T)+
                 dgamma(lambdatheta,shape=b4*new.lambdatheta,rate=b4,log = T))-
    
    ( bivdenseLambda(lambdabeta,lambdatheta)+
        dgamma(new.lambdabeta,shape=b3*lambdabeta,rate=b3,log = T)+
        dgamma(new.lambdatheta,shape=b4*lambdatheta,rate=b4,log = T))
  
  
  aux.ratio <- min(0,drop(aux.ratio))
  U <- runif(1)
  if(U<exp(aux.ratio)){
    lambdabeta =new.lambdabeta
    lambdatheta=new.lambdatheta
  }else{
    lambdabeta =lambdabeta
    lambdatheta=lambdatheta
  }
  # Update hyper-param a11,b11,a22,b22
  
  new.a11 <- runif(1,0.5,10)
  aux.ratio <- (dgamma(lambdabeta,shape=new.a11,rate=b11 ,log=TRUE)+
                  dunif(a11,0.5,10,log = TRUE))-
    (dgamma(lambdabeta,shape=a11,rate=b11 ,log=TRUE)+
       dunif(new.a11,0.5,10,log = TRUE))
  
  aux.ratio <- min(0,drop(aux.ratio))
  U <- runif(1)
  if(U<exp(aux.ratio)){
    a11=new.a11
  }else{
    a11=a11
  }
  
  new.b11 <- runif(1,0.5,10)
  aux.ratio <- (dgamma(lambdabeta,shape=a11,rate=new.b11 ,log=TRUE)+
                  dunif(b11,0.5,10,log = TRUE))-
    (dgamma(lambdabeta,shape=a11,rate=b11 ,log=TRUE)+
       dunif(new.b11,0.5,10,log = TRUE))
  
  aux.ratio <- min(0,drop(aux.ratio))
  U <- runif(1)
  if(U<exp(aux.ratio)){
    b11=new.b11
  }else{
    b11=b11
  }
  
  ######
  
  new.a22 <- runif(1,0.5,10)
  aux.ratio <- (dgamma(lambdatheta,shape=new.a22,rate=b22 ,log=TRUE)+
                  dunif(a22,0.5,10,log = TRUE))-
    (dgamma(lambdatheta,shape=a22,rate=b22 ,log=TRUE)+
       dunif(new.a22,0.5,10,log = TRUE))
  
  aux.ratio <- min(0,drop(aux.ratio))
  U <- runif(1)
  if(U<exp(aux.ratio)){
    a22=new.a22
  }else{
    a22=a22
  }
  
  new.b22 <- runif(1,0.5,10)
  aux.ratio <- (dgamma(lambdatheta,shape=a22,rate=new.b22 ,log=TRUE)+
                  dunif(b22,0.5,10,log = TRUE))-
    (dgamma(lambdatheta,shape=a22,rate=b22 ,log=TRUE)+
       dunif(new.b22,0.5,10,log = TRUE))
  
  aux.ratio <- min(0,drop(aux.ratio))
  U <- runif(1)
  if(U<exp(aux.ratio)){
    b22=new.b22
  }else{
    b22=b22
  }
  
  ### Results
  param[[1]]=xi
  param[[2]]=tau 
  param[[3]]=xi.prime
  param[[4]]=tau1
  param[[5]]=lambdabeta
  param[[6]]=lambdatheta
  param[[7]]=a11
  param[[8]]=b11
  param[[9]]=a22
  param[[10]]=b22
  
  return(param)
  
}


SampleLoc <- function(n){
  # Input
  # n: number of points 
  # Output
  # distmat: distance matrix
  # loc: location/long and lat/ x,y
  
  pts <- cbind(s1 = sample(1:n / n - 0.5 / n)^2,
               s2 = sample(1:n / n - 0.5 / n)^2)
  
  dmat <- as.matrix(dist(pts))
  return(list(distmat=dmat,loc=pts))
}

cMatern <- function(h, nu, kappa) {
  ## Input
  # h: euclidean distance
  # nu: parameter nu 
  # kappa: parameter kappa
  # Vectorize == T
  ## Output
  # Correlation coefficient.
  
  ifelse(h > 0, besselK(h * kappa, nu) * (h * kappa)^nu / 
           (gamma(nu) * 2^(nu - 1)), 1)
}

NearestN <- function(dmat,neigb){
  
  ## Input
  #  dmat: distance matrix
  #  neigb: number of neigbors
  ## Output: List of nearest neigbhors
  
  N=list()
  n=nrow(dmat)
  diag(dmat) <- NA
  
  for(i in 1 : n){
    roworder <-  (dmat[i,] %>% order) # order the row according to dist.
    roworderid <- roworder<i           # determin j < i
    roworder <-  roworder[roworderid]  # Select index of j <i
    
    if(i<=neigb){
      
      N[[i]] <- roworder
    }else {
      N[[i]] <- roworder[1:neigb]
    }
    
  } 
  return(N)
}


computeAD <- function(C,N){
  
  ## Input
  # C : dense covariance matrix
  # N: a list of Neigbhors
  ## Output
  ## Sparse matrix of A and D.
  
  n=nrow(C)
  A=D=matrix(0,n,n)
  
  for(i in 1:(n-1)) {
    A[i+1,N[[i+1]]] = solve(C[N[[i+1]],N[[i+1]]], C[N[[i+1]],i+1])
    D[i+1,i+1] = C[i+1,i+1] - sum(C[i+1,N[[i+1]]]*A[i+1,N[[i+1]]])
  }
  D[1,1]=C[1,1]
  A[1,]=0
  
  return(list(A=as(A,"sparseMatrix"),
              D=as(D,"sparseMatrix")))
}


# optimize loess function for interpolation
Optimizeloess <- function(ObData){
  
  ## Input
  # obData : data frame containing Y,x1,x2
  # Y is the model prediction and x1,x2 is a pair of coordinate
  ## Output
  ## Optimum parameter for span in loess(...,span)
  
  sse=99999
  calcSSE <- function(x){
    loessMod <- try(loess(Y ~ x1+x2, data=ObData, span=x), silent=T)
    res <- try(loessMod$residuals, silent=T)
    if(class(res)!="try-error"){
      if((sum(res, na.rm=T) > 0)){
        sse <- sum(res^2)  
      }
    }else{
      sse <- 99999
    }
    return(sse)
  }
  
  # Run optim to find span that gives min SSE, starting at 0.5
  res <- optim(par=c(0.5), calcSSE, method="SANN")
  
  return(res)
}

FullauxDIC <- function(XI,XIbar){
  XI=XI
FulauxDIC <- function(XI){
  # Auxiliary function to compute DIC
  # Input 
  # XI: Parameter vector for FULL model
  # Output
  # dic: D(u)=-2log(L(u))
  
  f <- function(j){
    
    D1  <-  V[1:(m*n),1:(m+p)]%*%XI[1:(m+p),j]
    L <- dmvnorm(t(matrix(Y[1:(m*n),1])),
                 matrix(D1),
                 Tau[j]^(-1)*diag(rep(1,m*n))
                 ,log = TRUE)
    return(-2*L)
  }
  rexDIC <- unlist(lapply(1:ncol(XI), f))
  return(rexDIC)
}
  
DIC <- 2*mean(FulauxDIC(XI))-FulauxDIC(matrix(XIbar))
 
return(DIC)
}


FullWAIC <- function(XI=XI){
  # Function to compute WAIC for full model
  # Input:
  # XI: Parameter
  
  auxestLPPD <- function(k){
    # Function to compute the Lppd for the kth replication
    # k : kth replication
    
    f<- function(j){
      D1  <-          V[(1:m+(k-1)*m),1:(m+p)]%*%XI[1:(m+p),j]
      D2  <- t(matrix(Y[(1:m+(k-1)*m),1]))
      L <- dmvnorm(D2,
                   matrix(D1),
                   Tau[j]^(-1)*diag(rep(1,m))
                   ,log = FALSE)
      return(L)
    }
    rexlppd <- unlist(lapply(1:ncol(XI), f))
    pwaic   <- var(log(rexlppd))
    rexlppd <-log(mean(rexlppd))
    
    return(list(rexlppd,
                pwaic))
  }
  
  EstLPPD=NULL
  Pwaic=NULL
  #EstLPPD <- unlist(lapply(1:n,auxestLPPD))
  for (k in 1:n) {
    rex = auxestLPPD(k)
    EstLPPD = c(EstLPPD,rex[[1]])
    Pwaic   = c(Pwaic,rex[[2]])
  }
  EstLPPD <- sum(EstLPPD)
  Pwaic   <- sum(Pwaic)
  return(-2*(EstLPPD-Pwaic))
}

generateLoc <- function(n,x){
  # Input
  # n: number of points 
  # x: spatial window
  # Output
  # loc: location/long and lat/ x,y
  
  pts <- cbind(s1 = (sample(seq(min(x),max(x), length=1000),n,replace = F)),
               s2 = (sample(seq(min(x),max(x), length=1000),n,replace = F)))
  
  return(pts)
}

# Simulate a bimodal mean field given two peaks
meanField = function(l,peak1,peak2) {
  # INPUT 
  # l: matrix of coordinate
  # peaks: vector of peak location
  
  peak1 =t(as.matrix(peak1))
  peak2 =t(as.matrix(peak2))
  
  aux1 = spDists(l, peak1,longlat = FALSE)
  aux2 = spDists(l, peak2,longlat = FALSE)
  v= unlist(lapply(1:nrow(l),function(k) min(aux1[k],aux2[k])))
  v= 10*exp(-v/0.5)
  v=v-mean(v)
  return(v)
}


CPO <- function(n,m,M,Y,V,XI,Tau){
  # CPO for each replication.
  # INPUT:
  # n Number of replication
  # m number of locations
  # M number of Monte Carlos draws
  # Y rbind response for all replications
  # V Rbind for all desing matrix including projection A
  # XI parameter vector from MCMC
  # Tau MCMC draws of Tau of the likelihood
  #OUTPUT
  #####
  CPO = vector("numeric",length = n)
  I   = diag(x=1,nrow = m) 
  for (k in 1:n) {
    aux = sample(1:ncol(XI),M)
    rex = 0
    for (l in 1:M) {
      
      eta= (V[(1:m)+(m*(k-1)),]%*%XI[,aux[l]])
      rex = rex+exp(-dmvnorm(t(matrix(Y[(1:m)+(m*(k-1)),])), mean = as.matrix(eta), sigma = (1/Tau[aux[l]])* I,log = T))
      
    }
    rex= rex/M
    rex= 1/rex
    CPO[k]=rex
  }
  return(CPO)
}

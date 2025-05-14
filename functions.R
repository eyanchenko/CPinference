#' @title Borgatti and Everett core-periphery metric
#' @description This function evaluates the BE CP metric
#' @param A adjacency matrix
#' @param C CP labels (C[i] = 1 if node i in core, and 0 otherwise)
#' @return value of BE metric of A at C
#' @export
obj.fun <- function(A,C){
  
  if(nrow(A)!=ncol(A)){
    stop("A must be a square matrix.")
  }
  
  if(nrow(A)!=length(C)){
    stop("Number of rows/columns in A must be the same as number of elements in C.")
  }
  
  if( sum(sort(unique(C))==c(0,1))!=2 ){
    stop("C must only contain 0's and 1's.")
  }
  
  n = length(C)
  CP <- C #idealized CP matrix to compare
  
  CP[C==1] = 2
  CP[C==0] = 1
  
  CP = CP%*%t(CP)%%2
  CP =  1 - CP
  diag(CP) <- 0
  
  #don't allow partitions of all core or all periphery
  if( sum(CP)==0 || sum(CP)==n*(n-1) ){
    rho=0
  }else{
    rho = cor(A[upper.tri(A)], CP[upper.tri(CP)])
  }
  
  return(rho)
}

#' @title Model-agnostic core-periphery population parameter
#' @description This function evaluates the CP population parameter
#' @param P data-generating probability matrix
#' @param C CP labels (C[i] = 1 if node i in core, and 0 otherwise)
#' @return value of CP populations parameter of P at C
#' @export
rho.fun <- function(P, C){
  if(nrow(P)!=ncol(P)){
    stop("P must be a square matrix.")
  }
  
  if(nrow(P)!=length(C)){
    stop("Number of rows/columns in P must be the same as number of elements in C.")
  }
  
  if( sum(sort(unique(C))==c(0,1))!=2 ){
    stop("C must only contain 0's and 1's.")
  }
  
  if(sum(C)==length(C) || sum(C)==0){
    stop("C must contain at least one 1 and one zero.")
  }
  
  n = length(C)
  k = sum(C)
  
  CP <- C #idealized CP matrix to compare
  
  CP[C==1] = 2
  CP[C==0] = 1
  
  CP = CP%*%t(CP)%%2
  CP =  1 - CP
  diag(CP) <- 0
  
  p     = mean(P[upper.tri(P)])
  d =  (0.5*k*(k-1)+(n-k)*k) / choose(n,2)
  
  
  rho = (0.5*sum(P*CP) - 0.5*n*(n-1)*p*d) / (0.5*n*(n-1)*sqrt(p*(1-p)*d*(1-d)))
  
  return(rho)
}

#' @title Intersection tests
#' @description This function carries out the intersection tests under either 
#' Erdos-Renyi (ER) or Chung-Lu (CL) null models.
#' @param A adjacency matrix
#' @param C CP labels (C[i] = 1 if node i in core, and 0 otherwise)
#' @param null null model. Must be either "ER" for Erdos-Renyi or "CL" for Chung-Lu
#' @return Test statistics, rejection thresholds and rejection decision
#' 0 means "fail to reject" while 1 means "reject"
#' @export
anal_rejc <- function(A, C, null=c("ER", "CL")){
  
  if(nrow(A)!=ncol(A)){
    stop("A must be a square matrix.")
  }
  
  if(nrow(A)!=length(C)){
    stop("Number of rows/columns in A must be the same as number of elements in C.")
  }
  
  if( sum(sort(unique(C))==c(0,1))!=2 ){
    stop("C must only contain 0's and 1's.")
  }
  
  if(sum(C)==length(C) || sum(C)==0){
    stop("C must contain at least one 1 and one zero.")
  }
  
  if(!(null %in% c("ER", "CL"))){
    stop("Null hypothesis must either be Erdos-Renyi (ER) or Chung-Lu (CL) model.")
  }
  
  # Find value of test stat
  n = dim(A)[1] # number of nodes
  k <- sum(C)
  alpha = k/n
  p  <- sum(A)/(n*(n-1))
  
  if(k*p < 1){
    stop("Network is too sparse and/or core is too small for asymptotic results to apply.")
  }
  
  if(null=="ER"){
    # Compute test statistics
    T1 <- obj.fun(A, C)
    C <- as.logical(C)
    p11h <- sum(A[C,C])/(k*(k-1))
    p12h <- sum(A[C,!C])/(k*(n-k))
    T2 <- p11h - p12h
    
    # Compute cutoffs
    C1 <- sqrt(log(k*p)/n) 
    C2 <- 2*sqrt(2)*p*log(n)/sqrt(k)
    
  }else if(null=="CL"){
    # Compute test statistics
    T1 <- obj.fun(A, C)
    
    theta  <- colSums(A)/sqrt(sum(A))
    Phat <- theta%*%t(theta)
    diag(Phat) <- 0
    C1 <- rho.fun(Phat, C)
    
    C <- as.logical(C)
    p11h <- sum(A[C,C])/(k*(k-1))
    p12h <- sum(A[C,!C])/(k*(n-k))
    T2 <- p11h - p12h
    
    ## Find cutoffs
    
    C1 = C1 + sqrt(alpha)*log(k) / (n^(1.5)*sqrt(p))
    C1 = C1 + sqrt(log(k*p))/n
    
    C2 = 2*sqrt(2)*p*log(n)/sqrt(k)
  }
  
  return(list(T1=T1, C1=C1, T2=T2, C2=C2, reject= as.numeric(T1 > C1 ) * as.numeric(T2 > C2 )) )
  
}


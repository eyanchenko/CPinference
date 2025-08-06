cppFunction('IntegerVector borgattiCpp(IntegerMatrix A){

  double n = A.nrow();
  double m = n*(n-1)/2;
  int k = 0;
  
  
  // Initialize vector 
  
  IntegerVector C(n);
  
  
  for(int i = 0; i < n; i++){
    
    if(rand()%10 > 5){
      C(i) = 1;
      k++;
    }else{
      C(i) = 0;
    }
     
  }
  
  int kk = k;
  
  
  double sum_A = 0; 
  double sum_CP = 0.5*k*(k-1)+(n-k)*k; 
  double sum_ACP = 0;
  
  for (int i = 1; i < n; i++){
      for(int j = 0; j < i; j++){
        sum_A   += A(i,j);
        sum_ACP += A(i,j) * (C(i) + C(j) - C(i)*C(j)); 
      }
    }
  
  double obj_last = (sum_ACP - sum_A / m * sum_CP) / 
    (sqrt(sum_A - sum_A/m*sum_A) * sqrt(sum_CP - sum_CP/m*sum_CP) );
  
                   
  int ind = 1;
  int n_iter = 0;
  int max_iter = 100;
  
  IntegerVector Ctest = clone(C);
  
  IntegerVector v(n);
  for(int i=0; i<n;i++){v(i) = i;}
  
  double obj_new = 0;
  
  while(ind > 0 & n_iter < max_iter){
  
  ind = 0;
  /* Randomly shuffle node order */
  int N = n;
  for(int a=0;a<n;a++){
    int b = round(rand()%N);
    int c = v(a); v(a)=v(b); v(b)=c;
  }
  
  
  
   for(int i = 0; i < n; i++){
    
     Ctest = clone(C);
     
     Ctest(v(i)) = 1 - Ctest(v(i));
     
    /* Update size of core for test vector, kk */
    
     if(Ctest(v(i))==0){
       kk = k - 1;
     }else{
       kk = k + 1;
     }

     double sum_CP_test = 0.5*kk*(kk-1) + (n-kk)*kk; 
     double sum_ACP_test = sum_ACP;
    
    /* Update value of objective function */
    
     for(int j = 0; j < n; j++){
        sum_ACP_test += A(v(i),j) * (Ctest(v(i)) + Ctest(j) - Ctest(v(i))*Ctest(j)); 
        sum_ACP_test -= A(v(i),j) * (C(v(i)) + C(j) - C(v(i))*C(j)); 
     }

    
     obj_new = (sum_ACP_test - sum_A / m * sum_CP_test) / 
      (sqrt(sum_A - sum_A/m*sum_A) * sqrt(sum_CP_test - sum_CP_test/m*sum_CP_test) );

    
     if(obj_new > obj_last){
       C(v(i)) = Ctest(v(i));
       sum_ACP = sum_ACP_test;
       ind++;
       obj_last = obj_new;
       k = kk;
      
       
     }
    
    
    }
    
    n_iter++;
    
  }
  
  return( C );
}')

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
  
  if(sum(C)==length(C) || sum(C)==0){
    stop("C must contain at least one 1 and one zero.")
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

generateP <- function(n, p11, p12, p22, prop=0.50){
  
  P = matrix(p22, ncol=n, nrow=n)
  k = round(n*prop)
  
  P[1:k, 1:k] <- p11
  P[1:k, (k+1):n] <- P[(k+1):n, 1:k] <- p12
  diag(P) <- 0
  
  return(P)
}
generatePDCBM <- function(theta, p11, p12, p22, prop=0.1){
  n = length(theta)
  P = matrix(p22, ncol=n, nrow=n)
  k = round(n*prop)
  
  P[1:k, 1:k] <- p11
  P[1:k, (k+1):n] <- P[(k+1):n, 1:k] <- p12
  
  P <- theta%*%t(theta) * P
  
  diag(P) <- 0
  
  P[P > 1] <- 1
  
  return(P)
}
generateA <- function(n, p11, p12, p22, prop=0.50){
  m = as.integer(round(n^2*prop^2, digits=0))
  A11 = matrix(rbinom(m, 1, p11), ncol=as.integer(round(n*prop)))
  A11[lower.tri(A11)] <- 0
  diag(A11) <- 0
  A11 = A11 + t(A11)
  
  m = as.integer(round(n^2*(1-prop)^2))
  A22 = matrix(rbinom(m, 1, p22), ncol=as.integer(round(n*(1-prop))), nrow = as.integer(round(n*(1-prop))))
  A22[lower.tri(A22)] <- 0
  diag(A22) <- 0
  A22 = A22 + t(A22)
  
  m = as.integer(round(n^2*prop*(1-prop)))
  A12 <- matrix(rbinom(m, 1, p12), ncol=as.integer(round(n*(1-prop))), nrow = as.integer(round(n*prop)))
  A21 <- t(A12)
  
  A = cbind(rbind(A11, A21), rbind(A12, A22))
  return(A)
}
generateDCBM <- function(theta, p11, p12, p22, q=0.1){
  #theta: different degrees
  #pij: block-block probabilities
  #q: proportion of nodes in block 1
  
  n <- length(theta)
  
  Omega11 <- matrix(p11, ncol = n*q, nrow = n*q)
  Omega22 <- matrix(p22, ncol = n*(1-q), nrow = n*(1-q))
  Omega12 <- matrix(p12, ncol = n*(1-q), nrow = n*q)
  
  Omega <- rbind(cbind(Omega11, Omega12), cbind(t(Omega12), Omega22) )
  
  P = as.vector(theta %*% t(theta) * Omega)
  
  P[P > 1] <- 1
  
  A <- matrix(rbinom(n^2, 1, P), ncol=n, nrow=n)
  
  A[lower.tri(A, diag=TRUE)] <- 0
  A = A + t(A)
  
  return(A) 
  
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

# From R Programs for Computing Truncated Distributions
qtrunc <- function(p, spec, a = -Inf, b = Inf, ...) {
  tt <- p
  G <- get(paste("p", spec, sep = ""), mode = "function")
  Gin <- get(paste("q", spec, sep = ""), mode = "function")
  tt <- Gin(G(a, ...) + p*(G(b, ...) - G(a, ...)), ...)
  return(tt)
}
rtrunc <- function(n, spec, a = -Inf, b = Inf, ...) {
  x <- u <- runif(n, min = 0, max = 1)
  x <- qtrunc(u, spec, a = a, b = b,...)
  return(x)
}

# Code up a Gibbs sampler with uniform I(p11 > p12 > p12) prior with Poisson approx

# Compute log-likelihood
llike <- function(A, p, C){
  n = dim(A)[1]
  k = sum(C)
  
  C = as.logical(C)
  
  n11 = 0.5*k*(k-1)
  n12 = k*(n-k)
  n22 = 0.5*(n-k)*(n-k-1)
  
  O11 = sum(A[C, C])/2
  O12 = sum(A[C, !C])
  O22 = sum(A[!C, !C])/2
  
  return(list(llike = -n11*p[1] + O11*log(p[1]) - n12*p[2] + O12*log(p[2]) - n22*p[3] + O22*log(p[3]),
              O = c(O11, O12, O22)))
}

llike_pot <- function(A, p, O, C1, C2){
  n = dim(A)[1]
  C1 = as.logical(C1)
  m = sum(O)

  k = sum(C2)
  
  C2 = as.logical(C2)
  
  n11 = 0.5*k*(k-1)
  n12 = k*(n-k)
  n22 = 0.5*(n-k)*(n-k-1)
  
  # Node which was changed
  idx = (1:n)[!(C1==C2)]
  
  if(C1[idx]){
    
    O[1] = O[1] - sum(A[C1, idx])
    O[3] = O[3] + sum(A[!C1, idx])
    O[2] = m - O[1] - O[3]
  }else{
    O[1] = O[1] + sum(A[C1, idx])
    O[3] = O[3] - sum(A[!C1, idx])
    O[2] = m - O[1] - O[3]    
  }
  
  return(list(llike = -n11*p[1] + O[1]*log(p[1]) - n12*p[2] + O[2]*log(p[2]) - n22*p[3] + O[3]*log(p[3]),
              O = c(O[1], O[2], O[3])))
  
} 

cppFunction('IntegerVector loopCpp(IntegerMatrix A, IntegerVector C, NumericVector p, IntegerVector OO, double llike){

  double n = A.nrow();
  double m = OO(0) + OO(1) + OO(2);
  IntegerVector C2 = clone(C);
  IntegerVector v(n);
  for(int i=0; i<n;i++){v(i) = i;}
  
  int k = 0;
      
  for(int i=0; i< n; i++){
    k += C(i);
  }
  
  double llike_cur = llike;

  
  /* Randomly shuffle node order */
  int N = n;
  for(int a=0;a<n;a++){
    int b = round(rand()%N);
    int c = v(a); v(a)=v(b); v(b)=c;
  }
  
   for(int x = 0; x < n; x++){
    
     C2 = clone(C);
     int idx = v(x);
     C2(idx) = 1 - C2(idx);
  
     // Compute new log-likelihood
      int k_new = k;
      
      if(C2(idx)==1){
        k_new++;
      }else{
        k_new--;
      }

  
      int n11 = k_new*(k_new-1)/2;
      int n12 = k_new*(n-k_new);
      int n22 = (n-k_new)*(n-k_new-1)/2;

      
      IntegerVector OO_prop(3);
      OO_prop(0) = OO(0);
      OO_prop(1) = OO(1);
      OO_prop(2) = OO(2);
          
      
      if(C(idx)==1){
       for(int i=0; i<n; i++){
        OO_prop(0) -= C(i)*A(i, idx);
        OO_prop(2) += (1-C(i))*A(i, idx);
       }
      }else{
       for(int i=0; i<n; i++){
        OO_prop(0) += C(i)*A(i, idx);
        OO_prop(2) -= (1-C(i))*A(i, idx);
       }
      }
      
      OO_prop(1) = m - OO_prop(0) - OO_prop(2);
     
      double llike_prop =  -n11*p(0) + OO_prop(0) *log(p(0)) - n12*p(1) + OO_prop(1) *log(p(1))- n22*p(2) + OO_prop(2) *log(p(2)); 
     
      double a = llike_prop - llike_cur;
      double r = ((double) rand() / (RAND_MAX));
      
      if(log(r) < a){ 
        C(idx) = 1 - C(idx);
        llike_cur = llike_prop;
        OO(0) = OO_prop(0);
        OO(1) = OO_prop(1);
        OO(2) = OO_prop(2);
        k = k_new;
      }
      
   }
  
  C.push_back(OO(0));
  C.push_back(OO(1));
  C.push_back(OO(2));
  
  return(C);
}')


runMCMC <- function(A, iters = 1000, burnin = 100){
  
  # Uniform prior on C
  # Uniform(0,1) prior on p11 > p12 > p22
  
  n = dim(A)[1]
  
  keepers <- matrix(0, ncol = 3 + n, nrow = iters)
  
  # Initialize
  C <- rbinom(n, 1, 0.1)
  p <- c(1.1*mean(A), mean(A), 0.9*mean(A))
  
  out <- llike(A, p, C)
  llike_cur = out$llike
  OO = out$O
  
  idx <- 1
  for(iter in 1:(burnin+iters)){
    
    ## update C
    # randomize node order
    # n.seq <- sample(1:n)
    # for(i in n.seq){
    #   
    #   C.prop <- C
    #   C.prop[i] <- 1 - C.prop[i] # switch label
    #   
    #   out_prop <- llike_pot(A, p, OO, C, C.prop)
    #   llike_prop <- out_prop$llike
    #   OO_prop    <- out_prop$O
    #   
    #   a <- min(0, llike_prop - llike_cur)
    #   
    #   if(log(runif(1)) < a){ # accept the move
    #     C[i] <- 1 - C[i]
    #     llike_cur = llike_prop
    #     OO = OO_prop
    #   }
    # }
    
    out <- loopCpp(A, C, p, OO, llike_cur)
    
    C   <- out[1:n]
    OO  <- out[(n+1):(n+3)]
    
    k   <- sum(C)
    n11 = 0.5*k*(k-1)
    n12 = k*(n-k)
    n22 = 0.5*(n-k)*(n-k-1)
    llike_cur <- -n11*p[1] + OO[1] *log(p[1]) - n12*p[2] + OO[2] *log(p[2])- n22*p[3] + OO[3] *log(p[3])
    
    ## update p
    
    p[1] <- rtrunc(1, "gamma", p[2], 1,    shape = OO[1] + 1, scale = 1/n11)
    p[2] <- rtrunc(1, "gamma", p[3], p[1], shape = OO[2] + 1, scale = 1/n12)
    p[3] <- rtrunc(1, "gamma", 0, p[2],    shape = OO[3] + 1, scale = 1/n22)
    
    C <- as.numeric(C)
    
    if(iter > burnin){
      keepers[idx,1:n] <- C
      keepers[idx, (n+1):(n+3)] <- p
      idx <- idx + 1
      
      if(idx%%1000==0){
        print(idx)
      }
      
    }
    
  }
  
  
  return(keepers)
  
}

# Classification accuracy
class_acc <- function(Cpred, Ctrue){
  sum(Cpred == Ctrue) / length(Cpred)
}


#source_python('~/Documents/Research/Srijan/Cp/cp_be.py')


C_rombach <- function(A, k){
  n = dim(A)[1]
  
  out <- unlist(cp_rombach(A))
  
  # ## Assign to core over clear delienation
  # M     <- matrix(0, nrow=n, ncol=3)
  # M[,1] <- 1:n   # node index
  # M[,2] <- out # proportion of samples in core
  # 
  # M <- M[order(M[,2],decreasing=TRUE),] # order by proportions
  # 
  # d   <- -diff(M[,2])  
  # idx <- which.max(d)  # find largest difference
  # M[1:idx,3] <- 1     # all nodes above this are assigned to core
  # 
  # M <- M[order(M[,1],decreasing=FALSE),]
  # 
  # C <- M[,3]
  
  ## Assign to core over clear delienation
  M     <- matrix(0, nrow=n, ncol=3)
  M[,1] <- 1:n   # node index
  M[,2] <- out # proportion of samples in core
  
  M <- M[order(M[,2],decreasing=TRUE),] # order by proportions
  
  # Assign top k nodes to core
  M[1:k, 3] <- 1
  M <- M[order(M[,1],decreasing=FALSE),]
  C <- M[,3]
  
  return(C)
}



C_cur <- function(A){
  unlist(cp_cur(A))
}

C_SBM <- function(A, iters, burn){
  n = dim(A)[1]
  out <- runMCMC(A, iters, burn) 
  
  out <- colSums(out[,1:n])/iters
  
  C <- numeric(n)
  C[out > 0.5] <- 1
  
  return(C)
  
}














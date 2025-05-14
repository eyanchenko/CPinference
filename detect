# Rcpp code to detect CP structure using greedy algorithm (Algorithm 1 from the manuscript).

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

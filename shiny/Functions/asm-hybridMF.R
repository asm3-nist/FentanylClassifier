asm_HybridMF <- function(Q,L,dM,a){
  
  M = dM;
  
  if(M>0){  # Query MW larger than Library MW
    A = c(Q,numeric(M)); 
    B = c(L,numeric(M)); 
    H = B;
  } else {  # Library MW larger than Query MW
    A = c(numeric(abs(M)),Q); 
    B = c(numeric(abs(M)),L); 
    H = B;  
  }
  maxmz = length(H);
  
  H0 = H;
  hMF0 = asm_SimpleMF(A,H);
  
  # plot(A,type="h",lwd=3,col=col1,main = paste("hMF = ",hMF0,sep=""),
  #      xlab="m/z",ylab="abundance",ylim=c(-32,32));
  # matplot(-H,type="h",lwd=3,col=col2,add=TRUE)
  # 
  
  if(M==0){
    result = c(hMF0,H^(1/a));
    return(result)
  }
  
  
  j = 1;
  k = 50;
  iter = 0;
  
  
  while(j <= k && iter < 100){
  
    alpha = abs(H-A);          # a vector of intensity discrepencies
    k = sum(alpha>10);
    if(k<=0){
      hMF = hMF0
      H = H0;
      break
    }
    
    
    omz = order(-alpha)[j];    # mz value of largest intensity discrepency
    omz_M_minus1 = omz - 1L;    # mz value of M-1 isotope of largest intensity discrepency
    omz_M_plus1 = omz + 1L;     # mz value of M+1 isotope of largest intensity discrepency
    omz_M_plus2 = omz + 2L;     # mz value of M+2 isotope of largest intensity discrepency

    temp = asm_hMF_IntensitySplitter(H,A,omz,M)
    A = temp[1,];
    H = temp[2,];
    
    hMF = asm_SimpleMF(A,H);
    epsilon = hMF-hMF0;
    if (epsilon < 1e-16){
      hMF = hMF0
      H = H0;
    }
    
    temp = asm_hMF_IntensitySplitter(H,A,omz_M_minus1,M)
    A = temp[1,];
    H = temp[2,];
    
    hMF = asm_SimpleMF(A,H);
    epsilon = hMF-hMF0;
    if (epsilon < 1e-16){
      hMF = hMF0
      H = H0;
    }
    
    if((maxmz-omz_M_plus1)>0){
    temp = asm_hMF_IntensitySplitter(H,A,omz_M_plus1,M)
    A = temp[1,];
    H = temp[2,];
    hMF = asm_SimpleMF(A,H);
    epsilon = hMF-hMF0;
    if (epsilon < 1e-16){
      hMF = hMF0
      H = H0;
    }
    }
    
    if((maxmz-omz_M_plus2)>1){
    temp = asm_hMF_IntensitySplitter(H,A,omz_M_plus2,M)
    A = temp[1,];
    H = temp[2,];
    hMF = asm_SimpleMF(A,H);
    epsilon = hMF-hMF0;
    if (epsilon < 1e-16){
      hMF = hMF0
      H = H0;
    }
    
    
    if (epsilon < 1e-16){
      hMF = hMF0
      H = H0;
      j = j+1;
    } else {
      H0 = H;
      hMF0 = hMF;
      j = 1;
    }
    }
  
  iter = iter+1;  
  }
  
  
  
  if(M>0){
    H = H;
  } else {
    H = H[(abs(M)+1):length(H)];  
  }
    
  result = c(hMF,H^(1/a))
  return(result);
  
}
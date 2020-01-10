asm_hMF_IntensitySplitter <- function(H,A,x,M){
  
  Hprime = H[x] - A[x];
  
  if(Hprime>0){
    Hdonor = H[x] 
    Hreceiver = H[(x+M)]
    Hp = H[x] - A[x];
    
    a = A[(x+M)] - (Hreceiver + Hp);
    
    if (a > 0){
      H[x+M] = Hreceiver + Hp;
      H[x]  = Hdonor - Hp;  
    } else {
      H[x+M] = Hreceiver + (Hp + a);
      H[x] = Hdonor - (Hp + a);
    }
    
  } 

  result = rbind(A,H);
  
  return(result)
}
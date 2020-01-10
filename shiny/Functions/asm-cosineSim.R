asm_CosineSim <- function(Q,L){
  
  num = sum(Q*L);
  
  den1 = sqrt(sum(Q^2))
  
  den2 = sqrt(sum(L^2));
  
  MF = num/(den1*den2)
  
  return(MF)
}
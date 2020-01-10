asm_SimpleMF <- function(Q,L){
  
  num = (sum(Q*L))^2;
  
  den1 = sum(Q^2)
  
  den2 = sum(L^2);
  
  MF = 999*num/(den1*den2)
  
  MF = round(MF)
  
  return(MF)
}
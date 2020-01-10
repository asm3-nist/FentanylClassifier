asm_muIndex <- function(HitListF,HitListB,FentMap){
  
  topRow = c(999,HitListF[,hMF]);
  
  leftCol = cbind(HitListB[,hMF],FentMap);
  
  HSimMap = rbind(topRow,leftCol);
  
  HSimMap = 0.5*(HSimMap+t(HSimMap));
  
  HDistMap = 1 - HSimMap/999;
  
  D = HDistMap;
  N = dim(D)[1]
  M_index = numeric(N)
  for(i in 1:length(M_index)){
    M_index[i] = asm_CosineSim(D[,1],D[,i])
  }
  
  result = M_index[-1];
  result = round(result,2)
  return(result)
  
}
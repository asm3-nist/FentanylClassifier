asmFentMassEstimator <- function(q){
  
  Q_relevant = asm_TopPeakSelector(q,800);
  
  bp = max(Q_relevant[,mz]);
  
  bp = bp + 91;
  
  return(bp)
}

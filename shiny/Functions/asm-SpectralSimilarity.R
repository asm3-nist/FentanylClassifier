asm_SpectralSimilarityMeasures <- function(q,mwq,l,mwl,threshold_ab,ab_scaling,wplot){ 
  
  mwq = as.numeric(mwq);
  mwl = as.numeric(mwl);

  DeltaMass = mwq - mwl;
  
  Q_relevant = asm_TopPeakSelector(q,threshold_ab);
  L_relevant = asm_TopPeakSelector(l,threshold_ab);

  max_mz = max(Q_relevant[,mz],L_relevant[,mz]);
  
  Q_padded = asm_SpectraPadder(Q_relevant,max_mz);
  L_padded = asm_SpectraPadder(L_relevant,max_mz);
  
  Q_padded_hat = Q_padded^ab_scaling
  L_padded_hat = L_padded^ab_scaling
  
  sMF = asm_SimpleMF(Q_padded_hat,L_padded_hat);
  
  HybridResults = asm_HybridMF(Q_padded_hat,L_padded_hat,DeltaMass,ab_scaling)
  hMF = HybridResults[1];
  
  H = (HybridResults[2:length(HybridResults)])
  
  if(wplot == 1){ asm_HeadTailSpectraPlotter(Q_padded,L_padded,H) }
  
 Result = c(hMF,sMF,DeltaMass);
 return(Result);
}

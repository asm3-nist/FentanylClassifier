asm_TopPeakSelector <- function(x,MF1){ 

  a = strsplit(x[[1]][1],";")[[1]]
  b = length(a)
  mz = numeric(b);
  ab = numeric(b)
  
  for(i in 1:b){
    c = strsplit(a[i]," ")[[1]];
    d = length(c);
    if(d>1){
    mz[i] = as.numeric(strsplit(a[i]," ")[[1]][d-1]);
    ab[i] = as.numeric(strsplit(a[i]," ")[[1]][d]);
    }
  }
  
  intermediateTable = as.data.table(cbind(mz,ab));
  intermediateTable = intermediateTable[order(-ab)];
  return(intermediateTable[ab>=MF1])
}




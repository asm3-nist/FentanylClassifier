simpleAlg_mw = function(x){ 
  
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
  
  paddedPeaks = asm_SpectraPadder(cbind(mz,ab),max(mz));
  
  peaks = cbind(seq(1:max(mz)),paddedPeaks);
  
  
  lowest_ab = 3;
  
  numPeaks = dim(peaks)[1];
  
  if (length(numPeaks)== 0) { return(peaks[1]); }
  
  for(i in numPeaks:2){

    m0 = peaks[i,1];
    a0 = peaks[i,2];
    
    m1 = peaks[(i-1),1];
    md1 = m0 - m1;
    a1 = peaks[(i-1),2];
    
    if((i-2) < 1){
      md2 = 0;
    } else{
      m2 = peaks[(i-2),1];
      md2 = m0 - m2;  
      a2 = peaks[(i-2),2];
    }
    
    if (a0 < lowest_ab) next
    
    
    if (md1 > 2){ return(m0)} # no peaks within 2 Da
    
    
    if (md1 == 2){ # possible Cl Br peaks
      
      if (a1 < as.integer(a0/2)){ return(m0)}
      
    } else if (md1 == 1){ # check for c13 isotope peaks
      
      isocalc = (a1*0.011*m1) / 14;
      
      if(as.integer(3*isocalc) > a0 | (a0 - as.integer(isocalc)) < lowest_ab){
        if(i==2){ 
          return(m1);}
        next;
      }
      
      if(md2 > 2){
        return(m0);
      } else {
        if (i-2 >= 1){
          if (a2 < as.integer(a0/2)){
            return(m0);
          }
        }
      }
    }
  }
  
  return(peaks[numPeaks,1]);
}
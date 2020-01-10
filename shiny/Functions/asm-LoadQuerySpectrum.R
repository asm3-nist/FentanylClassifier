asm_load_QuerySpectrum <- function(InpDirName,values){
  
  data = readLines(InpDirName)
  
  iName = grep("Name: ",data);
  if(length(iName[[1]])==0){
    qName = "Unknown";
  } else {
    qName = strsplit(data[iName],": ")[[1]][2];  
  }
  
  sdata = data[-grep(":",data)]
  
  a = length(sdata);
  c = sdata[a];
  d = strsplit(c," ");
  
  while(length(d[[1]])==0){
    a = a - 1;
    c = sdata[a];
    d = strsplit(c," ");
  }
  
  sdata = sdata[1:a];
  if(length(grep(";",sdata[1]))==0)
  {
    tdata = cat(sdata[1],";",sep="");
    for(i in 2:a){
      tdata = paste(tdata," ",sdata[i],";",sep="")
    }
    sdata = tdata
  }
  q = asm_PeakListCreator(sdata)
  
  iMW = grep("MW: ",data);
  if(length(iMW)==0){
    #mwq = simpleAlg_mw(q)+91; # estimate
    mwq = asmFentMassEstimator(q); # estimate
    values$mwisest <- 1; 
  } else {
    mwq = as.numeric(strsplit(data[iMW],": ")[[1]][2])
    values$mwisest <- 0;
  }
  
  a = c(mwq,q);

  return(a);
  
}
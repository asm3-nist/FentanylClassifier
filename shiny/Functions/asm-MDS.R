asm_MDS <- function(HitListF,HitListB,FentMap){
  
topRow = c(999,HitListF[,hMF]);

leftCol = cbind(HitListB[,hMF],FentMap);

HSimMap = rbind(topRow,leftCol);

HSimMap = 0.5*(HSimMap+t(HSimMap));

HDistMap = 1 + 1e-8 - HSimMap/999;

fit = isoMDS(HDistMap,k=2);
x = fit$points[,1];
y = fit$points[,2];

result = cbind(x,y);
return(result)

}
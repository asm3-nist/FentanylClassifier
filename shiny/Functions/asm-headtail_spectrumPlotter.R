asm_HeadTailSpectraPlotter <- function(Q_padded,L_padded,H){
  
  source("Routines/asm-colours.R")
  
  a = 10*floor(which(Q_padded>0)[1]/10);
  b = 10*floor(which(H>0)[1]/10);
  h = 10*floor(which(L_padded>0)[1]/10)
  c = min(a,b,h);
  
  e = 10*ceiling(which(Q_padded>0)[length(which(Q_padded>0))]/10);
  f = 10*ceiling(which(H>0)[length(which(H>0))]/10);
  g = 10*ceiling(which(L_padded>0)[length(which(L_padded>0))]/10)
  
  d = max(e,f,g)
  
  plot(Q_padded,type="h",lwd=3,col=colb,
       xlab="m/z",ylab="abundance",ylim=c(-1000,1000),xlim=c(c,d),axes=FALSE);
  axis(1,label=seq(c,d,20),at=seq(c,d,20),cex=0.7)
  axis(2,label=c(100,50,0,50,100),at=c(-1000,-500,0,500,1000),cex=0.7)
  abline(h=0,lwd=1,col="black")
  
  # textx = which(Q_padded > 100);
  # textx1 = numeric(length(textx));
  # textx1[1] = 1
  # 
  # for(i in 2:length(textx1)){
  #   if((textx[i] - textx[(i-1)]) > 1 ){
  #       textx1[i] = 1
  #   }
  # }
  # 
  # for(i in (length(textx1)):2){
  #   if(Q_padded[textx[i]] > Q_padded[textx[i-1]]){
  #       textx1[i] = 1;
  #       textx1[i-1] = 0;
  #   }
  # }
  # 
  # textx = textx*textx1
  # textx = textx[which(textx>0)]
  # texty = Q_padded[textx]
  # text(textx,texty+50,textx,cex=0.5)
  
  
  
  matplot(-L_padded,type="h",lwd=3,col=col2,add=TRUE)
  matplot(-H,type="h",lwd=3,col=col1,add=TRUE)
  
  # textxl = which(L_padded > 100);
  # textxh = which(H > 100)
  # textxL = setdiff(textxl,textxh)
  
 #  textx=textxL;
 #  textx1 = numeric(length(textx));
 #  textx1[1] = 1
 #  
 #  for(i in 2:length(textx1)){
 #    if((textx[i] - textx[(i-1)]) > 1 ){
 #      textx1[i] = 1
 #    }
 #  }
 #  
 #  for(i in (length(textx1)):2){
 #    if(L_padded[textx[i]] > L_padded[textx[i-1]]){
 #      textx1[i] = 1;
 #      textx1[i-1] = 0;
 #    }
 #  }
 #  
 #  textx = textx*textx1
 #  textx = textx[which(textx>0)]
 #  textxL = textx;
 #  
 #  
 #  textx = textxh;
 #  textx1 = numeric(length(textx));
 #  textx1[1] = 1
 #  
 #  for(i in 2:length(textx1)){
 #    if((textx[i] - textx[(i-1)]) > 1 ){
 #      textx1[i] = 1
 #    }
 #  }
 #  
 #  for(i in (length(textx1)):2){
 #    if(H[textx[i]] > H[textx[i-1]]){
 #      textx1[i] = 1;
 #      textx1[i-1] = 0;
 #    }
 #  }
 #  
 #  textx = textx*textx1
 #  textx = textx[which(textx>0)]
 #  textxh = textx;
 # 
 # if(length(textxL!=0)){
 #       textyL = -L_padded[textxL]
 #       text(textxL,(textyL-50),textxL,cex=0.5,col=col2)
 #  }
 #  
 #  textyh = -H[textx]
 #  text(textxh,(textyh-50),textxh,cex=0.5,col=colb)

  
}
asm_strucPlotter_2 <-function(B,A,C,D,E){
  
  #par(mai=c(1.02,0.82,0.82,0.42))
  par(mai=c(0.1,0.1,0.1,0.1))
  
  # COMMON SCAFFOLD
  x = c(0,1,1,2,3,3,2,-1,-1,-2)
  y = c(0,-1,-3,0,-1,-3,-4,-1,-3,0)
  l = c("N","C","C","C","C","N","C","C","O","C");
  
  ox = x;
  oy = y;
  ol = l;
  
  plot(x,y,xlim=c(-6,10),ylim=c(-8,9.5),
       pch=15,cex=3,col="white",
       axes=FALSE,xlab="",ylab="")
  
  epsilon = 0.05
  segments(x[1],y[1],x[2],y[2])
  segments(x[2],y[2],x[3],y[3])
  segments(x[2],y[2],x[4],y[4])
  segments(x[4],y[4],x[5],y[5])
  segments(x[5],y[5],x[6],y[6])
  segments(x[6],y[6],x[7],y[7])
  segments(x[3],y[3],x[7],y[7])
  segments(x[1],y[1],x[8],y[8])
  segments(x[8]-epsilon,y[8],x[9]-epsilon,y[9])
  segments(x[8]+epsilon,y[8],x[9]+epsilon,y[9])
  segments(x[8],y[8],x[10],y[10])
  
  
  
  ## SITE B STRUCTURE
  
  bx = as.numeric(B$a_blockInfo[[1]][,1]);
  by = as.numeric(B$a_blockInfo[[1]][,2]);
  bl = B$a_blockInfo[[1]][,3];
  
  bi = as.numeric(B$b_blockInfo[[1]][,1]);
  bj = as.numeric(B$b_blockInfo[[1]][,2]);
  bt = as.numeric(B$b_blockInfo[[1]][,3]);
  
  txb = x[6]+bx;
  tyb = y[6]+by;
  
  ox = c(ox,txb);
  oy = c(oy,tyb);
  ol = c(ol,bl)
  
  segments(x[6],y[6],txb[1],tyb[1])
  for(i in 1:length(bi)){
    if(bt[i]==2){
      if(bx[bi[i]]==bx[bj[i]]){ # vertical bond
        x1 = txa[bi[i]]+epsilon;
        x2 = txa[bi[i]]-epsilon;
        segments(x1,tya[bi[i]],x1,tya[bj[i]],col="blue")
        segments(x2,tya[bi[i]],x2,tya[bj[i]],col="blue")
      } else if(by[bi[i]]==by[bj[i]]){ # horizontal bond
        y1 = tya[bi[i]]+epsilon;
        y2 = tya[bi[i]]-epsilon;
        segments(txa[bi[i]],y1,txa[bj[i]],col="blue")
        segments(txa[bi[i]],y2,txa[bj[i]],col="blue")
      } else {
        y1 = tya[bi[i]]+1.5*epsilon;
        y2 = tya[bj[i]]+1.5*epsilon;
        segments(txa[bi[i]],y1,txa[bj[i]],y2,col="blue")
        y1 = tya[bi[i]]-1.5*epsilon;
        y2 = tya[bj[i]]-1.5*epsilon;
        segments(txa[bi[i]],y1,txa[bj[i]],y2,col="blue")
      }
    } else {
      segments(txb[bi[i]],tyb[bi[i]],txb[bj[i]],tyb[bj[i]]) 
    }
  }

  
  
  ## SITE A STRUCTURE
  
  bx = as.numeric(A$a_blockInfo[[1]][,1]);
  by = as.numeric(A$a_blockInfo[[1]][,2]);
  bl = A$a_blockInfo[[1]][,3];
  
  bi = as.numeric(A$b_blockInfo[[1]][,1]);
  bj = as.numeric(A$b_blockInfo[[1]][,2]);
  bt = as.numeric(A$b_blockInfo[[1]][,3]);
  
  i = length(txb);
  if(txb[i]%%2==1){
    txa = txb[i]+bx;
    tya = tyb[i]+by;  
  } else {
    txa = txb[i]+bx;
    tya = tyb[i]-by;
  }
  
  
  ox = c(ox,txa);
  oy = c(oy,tya);
  ol = c(ol,bl)
  
  segments(txb[i],tyb[i],txa[1],tya[1])
  for(i in 1:length(bi)){
    if(bt[i]==2){
      if(bx[bi[i]]==bx[bj[i]]){ # vertical bond
        x1 = txa[bi[i]]+epsilon;
        x2 = txa[bi[i]]-epsilon;
        segments(x1,tya[bi[i]],x1,tya[bj[i]],col="blue")
        segments(x2,tya[bi[i]],x2,tya[bj[i]],col="blue")
      } else if(by[bi[i]]==by[bj[i]]){ # horizontal bond
        y1 = tya[bi[i]]+epsilon;
        y2 = tya[bi[i]]-epsilon;
        segments(txa[bi[i]],y1,txa[bj[i]],col="blue")
        segments(txa[bi[i]],y2,txa[bj[i]],col="blue")
      } else {
        y1 = tya[bi[i]]+1.5*epsilon;
        y2 = tya[bj[i]]+1.5*epsilon;
        segments(txa[bi[i]],y1,txa[bj[i]],y2,col="blue")
        y1 = tya[bi[i]]-1.5*epsilon;
        y2 = tya[bj[i]]-1.5*epsilon;
        segments(txa[bi[i]],y1,txa[bj[i]],y2,col="blue")
       }
    } else {
      segments(txa[bi[i]],tya[bi[i]],txa[bj[i]],tya[bj[i]]) 
    }
  }
    
  
  
  # SITE C STRUCTURE
    
  bx = as.numeric(C$a_blockInfo[[1]][,1]);
  by = as.numeric(C$a_blockInfo[[1]][,2]);
  bl = C$a_blockInfo[[1]][,3];
  
  bi = as.numeric(C$b_blockInfo[[1]][,1]);
  bj = as.numeric(C$b_blockInfo[[1]][,2]);
  bt = as.numeric(C$b_blockInfo[[1]][,3]);
  bLoc = as.numeric(C$b_blockInfo[[1]][,4]);  #UNIQUE TO C STRUCTURES
  
  txc = x[bLoc]+bx;
  tyc = y[bLoc]+by;
  
  ox = c(ox,txc);
  oy = c(oy,tyc);
  ol = c(ol,bl)   
    
  segments(x[bLoc],y[bLoc],txc[1],tyc[1])  
  if(is.na(bi[1])==FALSE){
    for(i in 1:length(bi)){
      if(bt[i]==2){
        if(bx[bi[i]]==bx[bj[i]]){ # vertical bond
          x1 = txc[bi[i]]+epsilon;
          x2 = txc[bi[i]]-epsilon;
          segments(x1,tyc[bi[i]],x1,tyc[bj[i]],col="blue")
          segments(x2,tyc[bi[i]],x2,tyc[bj[i]],col="blue")
        } else if(by[bi[i]]==by[bj[i]]){ # horizontal bond
          y1 = tyc[bi[i]]+epsilon;
          y2 = tyc[bi[i]]-epsilon;
          segments(txc[bi[i]],y1,txc[bj[i]],col="blue")
          segments(txc[bi[i]],y2,txc[bj[i]],col="blue")
        } else {
          y1 = tyc[bi[i]]+1.5*epsilon;
          y2 = tyc[bj[i]]+1.5*epsilon;
          segments(txc[bi[i]],y1,txc[bj[i]],y2,col="blue")
          y1 = tyc[bi[i]]-1.5*epsilon;
          y2 = tyc[bj[i]]-1.5*epsilon;
          segments(txc[bi[i]],y1,txc[bj[i]],y2,col="blue")
        }
      } else {
        segments(txc[bi[i]],tyc[bi[i]],txc[bj[i]],tyc[bj[i]]) 
      }
    }
  }
  
  
  
  # SITE D STRUCTURE
  
  bx = as.numeric(D$a_blockInfo[[1]][,1]);
  by = as.numeric(D$a_blockInfo[[1]][,2]);
  bl = D$a_blockInfo[[1]][,3];
  
  bi = as.numeric(D$b_blockInfo[[1]][,1]);
  bj = as.numeric(D$b_blockInfo[[1]][,2]);
  bt = as.numeric(D$b_blockInfo[[1]][,3]);

  
  txd = x[1]+bx;
  tyd = y[1]+by;
  
  ox = c(ox,txd);
  oy = c(oy,tyd);
  ol = c(ol,bl)   
  
  segments(x[1],y[1],txd[1],tyd[1])  
  if(is.na(bi[1])==FALSE){
    for(i in 1:length(bi)){
      if(bt[i]==2){
        if(bx[bi[i]]==bx[bj[i]]){ # vertical bond
          x1 = txd[bi[i]]+epsilon;
          x2 = txd[bi[i]]-epsilon;
          segments(x1,tyd[bi[i]],x1,tyd[bj[i]],col="blue")
          segments(x2,tyd[bi[i]],x2,tyd[bj[i]],col="blue")
        } else if(by[bi[i]]==by[bj[i]]){ # horizontal bond
          y1 = tyd[bi[i]]+epsilon;
          y2 = tyd[bi[i]]-epsilon;
          segments(txd[bi[i]],y1,txd[bj[i]],col="blue")
          segments(txd[bi[i]],y2,txd[bj[i]],col="blue")
        } else {
          y1 = tyd[bi[i]]+1.5*epsilon;
          y2 = tyd[bj[i]]+1.5*epsilon;
          segments(txd[bi[i]],y1,txd[bj[i]],y2,col="blue")
          y1 = tyd[bi[i]]-1.5*epsilon;
          y2 = tyd[bj[i]]-1.5*epsilon;
          segments(txd[bi[i]],y1,txd[bj[i]],y2,col="blue")
        }
      } else {
        segments(txd[bi[i]],tyd[bi[i]],txd[bj[i]],tyd[bj[i]]) 
      }
    }
  }
  
  
  
  # SITE E STRUCTURE
  
  bx = as.numeric(E$a_blockInfo[[1]][,1]);
  by = as.numeric(E$a_blockInfo[[1]][,2]);
  bl = E$a_blockInfo[[1]][,3];
  
  bi = as.numeric(E$b_blockInfo[[1]][,1]);
  bj = as.numeric(E$b_blockInfo[[1]][,2]);
  bt = as.numeric(E$b_blockInfo[[1]][,3]);
  
  
  txe = x[10]+bx;
  tye = y[10]+by;
  
  ox = c(ox,txe);
  oy = c(oy,tye);
  ol = c(ol,bl)   
  
  segments(x[10],y[10],txe[1],tye[1])  
  if(is.na(bi[1])==FALSE){
    for(i in 1:length(bi)){
      if(bt[i]==2){
        if(bx[bi[i]]==bx[bj[i]]){ # vertical bond
          x1 = txe[bi[i]]+epsilon;
          x2 = txe[bi[i]]-epsilon;
          segments(x1,tye[bi[i]],x1,tye[bj[i]],col="blue")
          segments(x2,tye[bi[i]],x2,tye[bj[i]],col="blue")
        } else if(by[bi[i]]==by[bj[i]]){ # horizontal bond
          y1 = tye[bi[i]]+epsilon;
          y2 = tye[bi[i]]-epsilon;
          segments(txe[bi[i]],y1,txe[bj[i]],col="blue")
          segments(txe[bi[i]],y2,txe[bj[i]],col="blue")
        } else {
          y1 = tye[bi[i]]+1.5*epsilon;
          y2 = tye[bj[i]]+1.5*epsilon;
          segments(txe[bi[i]],y1,txe[bj[i]],y2,col="blue")
          y1 = tye[bi[i]]-1.5*epsilon;
          y2 = tye[bj[i]]-1.5*epsilon;
          segments(txe[bi[i]],y1,txe[bj[i]],y2,col="blue")
        }
      } else {
        segments(txe[bi[i]],tye[bi[i]],txe[bj[i]],tye[bj[i]]) 
      }
    }
  }
  
  
  
  # ADD LABELS
  NoC = which(ol!="C");
  matpoints(ox[NoC],oy[NoC],pch=15,cex=3,col="white")
  text(ox[NoC],oy[NoC],ol[NoC],col="red")
  
  
}
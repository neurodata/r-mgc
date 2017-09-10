FindLargestRectangles <-function(I,crit,minSize){
  # originally in Matlab by Jaroslaw Tuszynski 
  # https://www.mathworks.com/matlabcentral/fileexchange/28155-inscribed-rectangle/content/FindLargestRectangles.m
  # finds largest rectangle regions within all points set to 1.
  # input: I       - B/W boolean matrix or output of FindLargestSquares
  #        minSize - [height width] - minimum width and height on regions of 
  #                  interest (used to restrict final choise)
  #        crit    - Optimazation Criteria parameters to optimize:
  #                   crit(1)*height + crit(2)*width + crit(3)*height*width
  # output: 
  #         C    - value of the optimization criteria "crit" calculated for 
  #                each pixel 
  #         W, H - for each pixel I(r,c) return height and width of the largest 
  #                all-white rectangle with its upper-left corner at I(r,c)
  #         M    - Mask the largest all-white rectangle of the image
  # converted to R by Cencheng Shen
  if (missing(crit)){
    crit=c(0,0,1); 
  }
  if (missing(minSize)){
    minSize=c(2,2); 
  }
  p=crit;
  nR=nrow(I);
  nC=ncol(I);
  if (minSize[1]<1){
    minSize[1]= floor(minSize[1]*nR);
  }
  if (minSize[2]<1){
    minSize[2]= floor(minSize[2]*nC);
  }
  
  if (max(I) - min(I)==1){
    S = FindLargestSquares(I);
  } else {
    S = I;
  }
  n=max(S);
  W = S; 
  H = S;
  C = ((p[1]+p[2]) + p[3]*S) * S; 
  minH = max(minSize[1],1);
  minW = max(minSize[2],1);
  
  hight2width = matrix(0,n+1,1);  
  for (r in (1 : nR)){
    hight2width = matrix(0,n+1,1);  
    for (c in (nC: 1)){
      s=S[r,c];
      if(s>0){
        MaxCrit=C[r,c];
        for (hight in (s:1)){
          width=hight2width[hight];
          width = max(width+1,s);
          hight2width[hight] = width;
          Crit = p[1]*hight + p[2]*width + p[3]*width*hight;
          if (Crit>MaxCrit){
            MaxCrit = Crit;    
            W[r,c]  = width;
            H[r,c]  = hight;
          }
        }
        C[r,c]  = MaxCrit;
      }
      hight2width[(s+1):(n+1)] = 0;
    }       
  } 
  
  width2height = matrix(0,n+1,1);  
  for (c in (1 : nC)){
    width2height = matrix(0,n+1,1);  
    for (r in (nR: 1)){
      s=S[r,c];
      if(s>0){
        MaxCrit=C[r,c];
        for (width in (s:1)){
          hight=width2height[width];
          hight = max(hight+1,s);
          width2height[width] = hight;
          Crit = p[1]*hight + p[2]*width + p[3]*width*hight;
          if (Crit>MaxCrit){
            MaxCrit = Crit;    
            W[r,c]  = width;
            H[r,c]  = hight;
          }
        }
        C[r,c]  = MaxCrit;
      }
      width2height[(s+1):(n+1)] = 0;
    }       
  }
  
  C[H<minH | W<minW ] = 0;
  M = matrix(FALSE,nR,nC);
  if (max(C)>0){
    pos=which(C==max(C));
    for (i in (1:length(pos))){
      r = ((pos[i]-1) %% nR) + 1
      c = floor((pos[i]-1) / nR) + 1
      M[ r:(r+H[r,c]-1), c:(c+W[r,c]-1) ] = TRUE;
    }
  }
  result=list(C=C,H=H,W=W,M=M);
  return(result);
}

FindLargestSquares <-function(I){
  nr=nrow(I);
  nc=ncol(I);
  S=(I>0);
  for (r in ((nr-1):1)){
    for (c in ((nc-1):1)){
      if (S[r,c]){
        a = S[r  ,c+1];
        b = S[r+1,c  ];
        d = S[r+1,c+1];
        S[r,c] = min(c(a,b,d)) + 1;
      }
    }
  }
  return(S)
}
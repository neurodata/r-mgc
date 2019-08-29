#' This function generates the 20 simulations.
#'
#' @param type is a number between 1 to 20, 1-5 are monotone / close to linear relationships, 6-19 are non-monotone and strongly nonlinear, 20 is independent
#' @param n is the sample size,
#' @param d is the dimensionality that is an integer no smaller than 1.
#' @param dependent specifies dependent or independent relationship by 1 or 0
#' @param noise specifies the amount of noise
#' @return result: list containing x and y, for which x is n*d and y is n*d or n*1 depending on the set-up.
#' @export
#'

library(MASS)
#source('GenerateSimulations.R')
GenerateSimulations = function(type,n,d,dependent, noise){

  eps= rnorm(n, 0, 1);
  A=matrix(0,d,1);
  for (i in (1:d)){
    A[i]=1/i;
  }

  x=matrix(0,n,d);
  for (i in (1:d)){
    x[,i]=runif(n, -1, 1);
  }
  xA=x%*%A;
  if (dependent==0){
    x=matrix(0,n,d);
    for (i in (1:d)){
      x[,i]=runif(n, -1, 1);
    }
  }

  ## 1. Linear
  switch(type,
         { #Linear
           y=xA+noise*eps;
         },
         { #Exponential
           x=matrix(0,n,d);
           for (i in (1:d)){
             x[,i]=runif(n, 0, 3);
           }
           xA=x%*%A;
           y=exp(xA)+10*noise*eps;
           if (dependent==0){
             x=matrix(0,n,d);
             for (i in (1:d)){
               x[,i]=runif(n, 0, 3);
             }
           }
         },
         { #Cubic
           y=128*(xA-1/3)^3+48*(xA-1/3)^2-12*(xA-1/3)+80*noise*eps;
         },
         { #Joint Normal
           Sigma=diag(d*2);
           rho=1/d/2;
           Sigma[1:d,(d+1):(2*d)]=rho;
           Sigma[(d+1):(2*d),1:d]=rho;
           x=mvrnorm(n, matrix(0,1,2*d), Sigma);
           y=x[,(d+1):(2*d)]+0.5*noise*matrix(eps,n,d);
           if (dependent==0){
             x=mvrnorm(n, matrix(0,1,2*d), Sigma);
           }
           x=as.matrix(x[,1:d]);
         },
         { #Step Function
           if (d>1){
             noise=1;
           }
           y=(xA>0)+1*noise*eps;
         },
         { #Quadratic
           y=(xA)^2+0.5*noise*eps;
         },
         { #W Shape
           z=matrix(0,n,d);
           for (i in (1:d)){
             z[,i]=runif(n, 0, 1);
           }
           y=4*( ( xA^2 - 1/2 )^2 + z%*%A/500 )+0.5*noise*eps;
         },
         { #Spiral
           if (d>1){
             noise=1;
           }
           cc=0.4;
             rx=runif(n, 0, 5);
             ry=rx;
             rx=matrix(rx,n,d);
             z=rx;
           x[,1]=cos(z[,1]*pi);
           if (d>1){
           for (i in (1:(d-1))){
             x[,i+1]=x[,i]*cos(z[,i+1]*pi);
             x[,i]=x[,i]*sin(z[,i+1]*pi);
           }}
           x=rx*x;
           y=ry*sin(z[,1]*pi);
             y=y+cc*(d)*noise*mvrnorm(n, 0, 1);
           if (dependent==0){
               rx=runif(n, 0, 5);
               rx=matrix(rx,n,d);
               z=rx;
             x[,1]=cos(z[,1]*pi);
             if (d>1){
             for (i in (1:(d-1))){
               x[,i+1]=x[,i]*cos(z[,i+1]*pi);
               x[,i]=x[,i]*sin(z[,i+1]*pi);
             }}
             x=rx*x;
           }
         },
         { #Uncorrelated Binomial
           if (d>1){
             noise=1;
           }
           x=matrix(0,n,d);
           for (i in (1:d)){
             x[,i]=rbinom(n, 1, 0.5);
           }
           x=x+0.5*noise*mvrnorm(n, matrix(0,1,d), diag(d));
           y=(rbinom(n, 1, 0.5)*2-1);
           y=x%*%A*y+0.5*noise*eps;
           if (dependent==0){
             x=matrix(0,n,d);
             for (i in (1:d)){
               x[,i]=rbinom(n, 1, 0.5);
             }
             x=x+0.5*noise*mvrnorm(n, matrix(0,1,d), diag(d));
           }
         },
         { #Log(X^2)
           x=mvrnorm(n, matrix(0,1,d), diag(d));
           y=log(x^2)+3*noise*matrix(eps,n,d);
           if (dependent==0){
             x=mvrnorm(n, matrix(0,1,d), diag(d));
           }
         },
         { #Fourth root
           y=abs(xA)^(0.25)+noise/4*eps;
         },
         { #Sine 1/2
           x=matrix(runif(n, -1, 1),n,d);
           if (noise>0 || d>1){
             x=x+0.02*(d)*mvrnorm(n, matrix(0,1,d), diag(d));
           }
           theta=4;cc=1;
           y=sin(theta*pi*x)+cc*noise*matrix(eps,n,d);
           if (dependent==0){
             x=matrix(runif(n, -1, 1),n,d);
             if (noise>0 || d>1){
               x=x+0.02*(d)*mvrnorm(n, matrix(0,1,d), diag(d));
             }
           }
         },
         { #Sine 1/8
           x=matrix(runif(n, -1, 1),n,d);
           if (noise>0 || d>1){
             x=x+0.02*(d)*mvrnorm(n, matrix(0,1,d), diag(d));
           }
           theta=16;cc=0.5;
           y=sin(theta*pi*x)+cc*noise*matrix(eps,n,d);
           if (dependent==0){
             x=matrix(runif(n, -1, 1),n,d);
             if (noise>0 || d>1){
               x=x+0.02*(d)*mvrnorm(n, matrix(0,1,d), diag(d));
             }
           }
         },
         { #Square
           u=matrix(runif(n, -1, 1),n,d);
           v=matrix(runif(n, -1, 1),n,d);
           theta=-pi/8;
           eps=0.05*(d)*mvrnorm(n, matrix(0,1,d), diag(d));
           x=u*cos(theta)+v*sin(theta)+eps;
           y=-u*sin(theta)+v*cos(theta);
           if (dependent==0){
             u=matrix(runif(n, -1, 1),n,d);
             v=matrix(runif(n, -1, 1),n,d);
             eps=0.05*(d)*mvrnorm(n, matrix(0,1,d), diag(d));
             x=u*cos(theta)+v*sin(theta)+eps;
           }
         },
         { #Two Parabolas
           y=( xA^2  + 2*noise*runif(n, 0, 1))*(rbinom(n,1,0.5)-0.5);
         },
         { #Circle
           if (d>1){
             noise=1;
           }
           cc=0.4;
           rx=matrix(1,n,d);
           z=matrix(0,n,d);
           for (i in (1:d)){
             z[,i]=runif(n, -1, 1);
           }
           ry=matrix(1,n,1);
           x[,1]=cos(z[,1]*pi);
           if (d>1){
           for (i in (1:(d-1))){
             x[,i+1]=x[,i]*cos(z[,i+1]*pi);
             x[,i]=x[,i]*sin(z[,i+1]*pi);
           }}
           x=rx*x;
           y=ry*sin(z[,1]*pi);
           x=x+cc*noise*rx*mvrnorm(n, matrix(0,1,d), diag(d));
           if (dependent==0){
             z=matrix(0,n,d);
             for (i in (1:d)){
               z[,i]=runif(n, -1, 1);
             }
             x[,1]=cos(z[,1]*pi);
             if (d>1){
             for (i in (1:(d-1))){
               x[,i+1]=x[,i]*cos(z[,i+1]*pi);
               x[,i]=x[,i]*sin(z[,i+1]*pi);
             }}
             x=rx*x;
             x=x+cc*noise*rx*mvrnorm(n, matrix(0,1,d), diag(d));
           }
         },
         { #Ecllipse
           if (d>1){
             noise=1;
           }
           cc=0.4;
           rx=matrix(5,n,d);
           z=matrix(0,n,d);
           for (i in (1:d)){
             z[,i]=runif(n, -1, 1);
           }
           ry=matrix(1,n,1);
           x[,1]=cos(z[,1]*pi);
           if (d>1){
             for (i in (1:(d-1))){
               x[,i+1]=x[,i]*cos(z[,i+1]*pi);
               x[,i]=x[,i]*sin(z[,i+1]*pi);
             }}
           x=rx*x;
           y=ry*sin(z[,1]*pi);
           x=x+cc*noise*rx*mvrnorm(n, matrix(0,1,d), diag(d));
           if (dependent==0){
             z=matrix(0,n,d);
             for (i in (1:d)){
               z[,i]=runif(n, -1, 1);
             }
             x[,1]=cos(z[,1]*pi);
             if (d>1){
               for (i in (1:(d-1))){
                 x[,i+1]=x[,i]*cos(z[,i+1]*pi);
                 x[,i]=x[,i]*sin(z[,i+1]*pi);
               }}
             x=rx*x;
             x=x+cc*noise*rx*mvrnorm(n, matrix(0,1,d), diag(d));
           }
         },
         { #Diamond
           u=matrix(runif(n, -1, 1),n,d);
           v=matrix(runif(n, -1, 1),n,d);
           theta=-pi/4;
           eps=0.05*(d)*mvrnorm(n, matrix(0,1,d), diag(d));
           x=u*cos(theta)+v*sin(theta)+eps;
           y=-u*sin(theta)+v*cos(theta);
           if (dependent==0){
             u=matrix(runif(n, -1, 1),n,d);
             v=matrix(runif(n, -1, 1),n,d);
             eps=0.05*(d)*mvrnorm(n, matrix(0,1,d), diag(d));
             x=u*cos(theta)+v*sin(theta)+eps;
           }
         },
         { #Multiplicative Noise
           x=mvrnorm(n, matrix(0,1,d), diag(d));
           y=mvrnorm(n, matrix(0,1,d), diag(d));
           y=x*y;
           if (dependent==0){
             x=mvrnorm(n, matrix(0,1,d), diag(d));
           }
         },
         { #Independent clouds
           z=matrix(0,n,d);
           for (i in (1:d)){
             z[,i]=rbinom(n, 1, 0.5);
           }
           x=mvrnorm(n, matrix(0,1,d), diag(d))/3+(z-0.5)*2;
           z=matrix(0,n,d);
           for (i in (1:d)){
             z[,i]=rbinom(n, 1, 0.5);
           }
           y=mvrnorm(n, matrix(0,1,d), diag(d))/3+(z-0.5)*2;
         }
  )

  result=list(x=x,y=y);
  return(result)
}

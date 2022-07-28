library(ncvreg)
library(GAGA)
library(mvtnorm)
rm(list = ls())

set.seed(1234)


Nlambda=100
p_size_list = c(300,600)
sample_size=1000
expnum = 10

Mean=0
Sd=1

rr = 0.5

TC_GAGA = NULL
TC_GAGA_QR = NULL
TC_ALASSO_CV = NULL
TC_SCAD_CV = NULL
TC_MCP_CV = NULL

Mean_TC_GAGA=NULL
Mean_TC_GAGA_QR=NULL
Mean_TC_ALASSO_CV=NULL
Mean_TC_SCAD_CV=NULL
Mean_TC_MCP_CV=NULL

for(ii in 1:length(p_size_list)){
  
  cat("iteration = ", ii, "\n")
  
  
  p_size=p_size_list[ii];
  cov_mat=matrix(1:p_size*p_size,p_size,p_size) 
  for(i in 1:p_size){for(j in 1:p_size) {if(i!=j){cov_mat[i,j]=0.5}else{cov_mat[i,j]=1}}}
  
  for(iter in 1:expnum){
    #cat("iter is:",iter,"\n");
  
    zeroNum = round(rr*p_size)
    ind1 = sample(1:p_size,p_size)
    ind2 = ind1[1:zeroNum]
    beta_true = runif(p_size,0,5)
    
    beta_true[ind2] = 0
    
    ind3 = ind1[(zeroNum+1):p_size]
    
    pos_true=ind3;
    pos_false=ind2;
   
    X=rmvnorm(n=sample_size, mean=rep(0, nrow(cov_mat)), sigma=cov_mat)
    
    raodong=rnorm(sample_size,mean=Mean,sd=Sd)
    
    y=X%*%beta_true+raodong
    
    ##GaGa
    mratio = 2
    timeb = proc.time();
    EW = GAGA(X,y,alpha = mratio, family = "gaussian")
    timee = proc.time();
    TC_GAGA[iter]=timee[3]-timeb[3];
    
    
    timeb = proc.time();
    EW2 = GAGA(X,y,alpha = mratio, family = "gaussian",QR_flag = TRUE)
    timee = proc.time();
    TC_GAGA_QR[iter]=timee[3]-timeb[3];
    
    
    
    ##Adaptive LASSO
    timeb = proc.time();
    LSE=lm(y~X-1)## Linear Regression to create the Adaptive Weights Vector
    weight=abs(LSE$coefficients)^1# Using gamma = 1
    XW=X%*%diag(weight)##
    fit_ALASSO=ncvreg(XW, y,family="gaussian", penalty="lasso",nlambda =Nlambda)
    timee = proc.time();
    TC_ALASSO_CV[iter]=timee[3]-timeb[3];
    
    
    ## SCAD 
    timeb = proc.time();
    cvfit_SCAD= cv.ncvreg(X, y,family="gaussian", penalty="SCAD",lambda.min=0.00001,nlambda =Nlambda, nfolds=10)
    timee = proc.time();
    TC_SCAD_CV[iter]=timee[3]-timeb[3];
    
    
    ## MCP
    timeb = proc.time();
    cvfit_MCP= cv.ncvreg(X, y,family="gaussian", penalty="MCP",lambda.min=0.00001,nlambda =Nlambda, nfolds=10)
    timee = proc.time();
    TC_MCP_CV[iter]=timee[3]-timeb[3];
    
    
    
  }#for(iter in 1:expnum)
  Mean_TC_GAGA[ii]=mean(TC_GAGA)
  Mean_TC_GAGA_QR[ii]=mean(TC_GAGA_QR)
  Mean_TC_ALASSO_CV[ii]=mean(TC_ALASSO_CV)
  Mean_TC_SCAD_CV[ii]=mean(TC_SCAD_CV)
  Mean_TC_MCP_CV[ii]=mean(TC_MCP_CV)
  
  
}

TimeCost=data.frame(Mean_TC_GAGA,Mean_TC_GAGA_QR,Mean_TC_ALASSO_CV,Mean_TC_SCAD_CV,Mean_TC_MCP_CV)
write.table(TimeCost,"timecost.txt",row.names = FALSE)


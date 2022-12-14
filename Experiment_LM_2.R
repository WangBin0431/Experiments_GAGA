library(ncvreg)
library(mvtnorm)
library(GAGA)
set.seed(1234)

Nlambda=100
p_size = 100
sample_size=500
expnum = 30

Mean=0
Sd=2
rr = 0.5 #Proportion of value zero in beta
R1=1
R2=5

cov_mat=matrix(1:p_size*p_size,p_size,p_size) 
for(i in 1:p_size){for(j in 1:p_size) {if(i!=j){cov_mat[i,j]=0.5}else{cov_mat[i,j]=1}}}

ERR_ALASSO=matrix(1:expnum*Nlambda,expnum,Nlambda)
ERR_SCAD=matrix(1:expnum*Nlambda,expnum,Nlambda)
ERR_MCP=matrix(1:expnum*Nlambda,expnum,Nlambda)

ACC_ALASSO=matrix(1:expnum*Nlambda,expnum,Nlambda)
ACC_SCAD=matrix(1:expnum*Nlambda,expnum,Nlambda)
ACC_MCP=matrix(1:expnum*Nlambda,expnum,Nlambda)

ERR_GAGA = NULL
ACC_GAGA = NULL
ERR_GAGA_QR = NULL
ACC_GAGA_QR = NULL
ERR_ALASSO_CV = NULL
ERR_SCAD_CV = NULL
ERR_MCP_CV = NULL
ACC_ALASSO_CV = NULL
ACC_SCAD_CV = NULL
ACC_MCP_CV = NULL

for(iter in 1:expnum){#
  cat("iter is:",iter,"\n");
  
  #Set true beta
  zeroNum = round(rr*p_size)
  ind = sample(1:p_size,zeroNum)
  beta_true = runif(p_size,0.0*R2,R2)
  beta_true[ind] = 0
  
  X=R1*rmvnorm(n=sample_size, mean=rep(0, nrow(cov_mat)), sigma=cov_mat)
  noise=rnorm(sample_size,mean=Mean,sd=Sd)
  
  y=X%*%beta_true+noise
  
  
  ##Adaptive LASSO
  LSE=lm(y~X-1)## Linear Regression to create the Adaptive Weights Vector
  weight=abs(LSE$coefficients)^1# Using gamma = 1
  XW=X%*%diag(weight)##
  # fit_ALASSO <- ncvreg(XW, y,family="gaussian", penalty="lasso",nlambda =Nlambd)
  if(iter==1){
    cvfit_ALASSO<- cv.ncvreg(XW, y,family="gaussian", penalty="lasso",lambda.min=0.00001,nlambda =Nlambda, nfolds=10)
    labmda_ALASSO=cvfit_ALASSO$lambda
  }
  cvfit_ALASSO<- cv.ncvreg(XW, y,family="gaussian", penalty="lasso",lambda.min=0.00001,lambda =labmda_ALASSO, nfolds=10)
  fit_ALASSO = cvfit_ALASSO$fit
  
  for(j in 1:Nlambda){
    tmp = fit_ALASSO$beta[,j][-1]
    tmp = weight*tmp
    ERR_ALASSO[iter,j] = norm(beta_true-tmp,type = '2')
    ACC_ALASSO[iter,j] = cal.w.acc(as.character(tmp!=0),as.character(beta_true!=0))
  }
  tmp = fit_ALASSO$beta[,cvfit_ALASSO$min][-1]
  tmp = weight*tmp
  ERR_ALASSO_CV[iter]=norm(beta_true-tmp,type = '2')
  ACC_ALASSO_CV[iter]= cal.w.acc(as.character(tmp!=0),as.character(beta_true!=0))
  
  
  ## SCAD 
  # fit_SCAD <- ncvreg(X, y,family="gaussian", penalty="SCAD",nlambda =Nlambda)
  if(iter==1){
    cvfit_SCAD<- cv.ncvreg(X, y,family="gaussian", penalty="SCAD",lambda.min=0.00001,nlambda =Nlambda, nfolds=10)
    labmda_SCAD=cvfit_SCAD$lambda
  }
  cvfit_SCAD <- cv.ncvreg(X, y,family="gaussian", penalty="SCAD",lambda.min=0.00001,lambda =labmda_SCAD, nfolds=10)
  fit_SCAD = cvfit_SCAD$fit
  for(j in 1:Nlambda){ERR_SCAD[iter,j]=norm(beta_true-fit_SCAD$beta[,j][-1],type = '2')}
  for(j in 1:Nlambda){
    ACC_SCAD[iter,j]= cal.w.acc(as.character(fit_SCAD$beta[,j][-1]!=0),as.character(beta_true!=0))
  }
  tmp = fit_SCAD$beta[,cvfit_SCAD$min][-1]
  ERR_SCAD_CV[iter]=norm(beta_true-tmp,type = '2')
  ACC_SCAD_CV[iter]=cal.w.acc(as.character(tmp!=0),as.character(beta_true!=0))
  
  
  ## MCP
  # fit_MCP <- ncvreg(X, y,family="gaussian", penalty="MCP",nlambda =Nlambda)
  if(iter==1){
    cvfit_MCP<- cv.ncvreg(X, y,family="gaussian", penalty="MCP",lambda.min=0.00001,nlambda =Nlambda, nfolds=10)
    labmda_MCP=cvfit_MCP$lambda
  }
  cvfit_MCP <- cv.ncvreg(X, y,family="gaussian", penalty="MCP",lambda.min=0.00001,lambda =labmda_MCP, nfolds=10)
  fit_MCP = cvfit_MCP$fit
  for(j in 1:Nlambda){ERR_MCP[iter,j]=norm(beta_true-fit_MCP$beta[,j][-1],type = '2')}#
  for(j in 1:Nlambda){
    ACC_MCP[iter,j]=cal.w.acc(as.character(fit_MCP$beta[,j][-1]!=0),as.character(beta_true!=0))
  }
  tmp = fit_MCP$beta[,cvfit_MCP$min][-1]
  ERR_MCP_CV[iter]=norm(beta_true-tmp,type = '2')
  pos1=which(tmp!=0);
  pos2=which(tmp==0);
  ACC_MCP_CV[iter]=cal.w.acc(as.character(tmp!=0),as.character(beta_true!=0))
  
  ##GAGA
  mratio = 2
  fit_gaga1 = GAGA(X,y,alpha = mratio,itrNum = 100, family = "gaussian")
  EW = fit_gaga1$beta
  ERR_GAGA[iter]=norm(beta_true-EW,type = '2')
  ACC_GAGA[iter]=cal.w.acc(as.character(EW!=0),as.character(beta_true!=0))
  
  fit_gaga2 = GAGA(X,y,alpha = mratio,itrNum = 100, family = "gaussian",QR_flag = TRUE)
  EW2 = fit_gaga2$beta
  ERR_GAGA_QR[iter]=norm(beta_true-EW2,type = '2')
  ACC_GAGA_QR[iter]=cal.w.acc(as.character(EW2!=0),as.character(beta_true!=0))
  
}#for(iter in 1:expnum)

Mean_ERR_ALASSO=colMeans(ERR_ALASSO)
Mean_ERR_SCAD=colMeans(ERR_SCAD)
Mean_ERR_MCP=colMeans(ERR_MCP)

Mean_ACC_ALASSO=colMeans(ACC_ALASSO)
Mean_ACC_SCAD=colMeans(ACC_SCAD)
Mean_ACC_MCP=colMeans(ACC_MCP)

Mean_ERR_GAGA = mean(ERR_GAGA)
Mean_ERR_GAGA_QR = mean(ERR_GAGA_QR)
Mean_ACC_GAGA = mean(ACC_GAGA)
Mean_ACC_GAGA_QR = mean(ACC_GAGA_QR)

Mean_ERR_ALASSO_CV = mean(ERR_ALASSO_CV)
Mean_ERR_SCAD_CV = mean(ERR_SCAD_CV)
Mean_ERR_MCP_CV = mean(ERR_MCP_CV)
Mean_ACC_ALASSO_CV = mean(ACC_ALASSO_CV)
Mean_ACC_SCAD_CV = mean(ACC_SCAD_CV)
Mean_ACC_MCP_CV = mean(ACC_MCP_CV)


Std_ACC_MCP=sqrt(apply(ACC_MCP, 2, var))
Std_ACC_SCAD=sqrt(apply(ACC_SCAD, 2, var))
Std_ACC_ALASSO=sqrt(apply(ACC_ALASSO, 2, var))
Std_ERR_MCP=sqrt(apply(ERR_MCP, 2, var))
Std_ERR_SCAD=sqrt(apply(ERR_SCAD, 2, var))
Std_ERR_ALASSO=sqrt(apply(ERR_ALASSO, 2, var))




ERR=c(Mean_ERR_ALASSO,Mean_ERR_SCAD,Mean_ERR_MCP,Mean_ERR_GAGA,Mean_ERR_GAGA_QR,Mean_ERR_ALASSO_CV,Mean_ERR_SCAD_CV,Mean_ERR_MCP_CV)
ACC=c(Mean_ACC_ALASSO,Mean_ACC_SCAD,Mean_ACC_MCP,Mean_ACC_GAGA,Mean_ACC_GAGA_QR,Mean_ACC_ALASSO_CV,Mean_ACC_SCAD_CV,Mean_ACC_MCP_CV)
Algorithms=factor(c(rep('ALASSO',Nlambda),rep('SCAD',Nlambda),rep('MCP',Nlambda),'GAGA','GAGA_QR','ALASSO_CV','SCAD_CV','MCP_CV'),
                  levels=c('GAGA','GAGA_QR','ALASSO','SCAD','MCP','ALASSO_CV','SCAD_CV','MCP_CV'))
ERR_ACC=data.frame(ERR,ACC,Algorithms)#



library(ggplot2)
g1=ggplot(ERR_ACC,aes(ERR,ACC,shape=Algorithms,color=Algorithms))+ geom_point(size=3)+scale_shape_manual(values=seq(0,15))#shape>6,
g1

ERR2=c(ERR_ALASSO_CV,ERR_SCAD_CV,ERR_MCP_CV,ERR_GAGA,ERR_GAGA_QR)
ACC2=c(ACC_ALASSO_CV,ACC_SCAD_CV,ACC_MCP_CV,ACC_GAGA,ACC_GAGA_QR)
Algorithms=factor(c(rep('ALASSO_CV',expnum),rep('SCAD_CV',expnum),rep('MCP_CV',expnum),
                    rep('GAGA',expnum),rep('GAGA_QR',expnum)),
                  levels=c('GAGA','GAGA_QR','ALASSO_CV','SCAD_CV','MCP_CV'))
ERR_ACC2=data.frame(ERR2,ACC2,Algorithms)#

g2=ggplot(ERR_ACC2, aes(x=Algorithms, y=ERR2,fill=Algorithms))+ylab("ERR") + geom_boxplot()#
g2

g3=ggplot(ERR_ACC2, aes(x=Algorithms, y=ACC2,fill=Algorithms))+xlab("Algorithms")+ylab("ACC") + geom_boxplot()#
g3

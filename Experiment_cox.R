# Demo of Gaga dealing with Cox model

setwd("I:/GAGA/R")
library(glmnet)
library(mvtnorm)

set.seed(1234)
source("calculate.R")
source("cox_GAGA.R")


Nlambda=100
p_size = 50
sample_size = 500
test_size = 1000

R1 = 3
R2 = 1
rate = 0.5 #Proportion of value zero in beta
censoringRate = 0.25 #Proportion of censoring data in observation data
Num = 10 # Total number of experiments

ERR_glmnet = NULL
ACC_glmnet = NULL
ERR_GAGA = NULL
ACC_GAGA = NULL
Cidx_glmnet = NULL
Cidx_GAGA = NULL

for(ii in 1:Num ){
  cat("iter is:",ii,"\n");
  #Set true beta
  zeroNum = round(rate*p_size)
  ind1 = sample(1:p_size,p_size)
  ind2 = ind1[1:zeroNum]
  beta_true = runif(p_size,-R2,R2)
  beta_true[ind2] = 0

  #Generate training samples
  cov_mat=matrix(1:p_size*p_size,p_size,p_size) ##covariance matrix
  for(i in 1:p_size){for(j in 1:p_size) {if(i!=j){cov_mat[i,j]=0.0}else{cov_mat[i,j]=1}}}
  X = R1*rmvnorm(n=sample_size, mean=rep(0, nrow(cov_mat)), sigma=cov_mat)
  z = X%*%beta_true
  u = runif(sample_size,0,1)
  t = ((-log(1-u)/(3*exp(z)))*100)^(0.1)
  cs = rep(0,sample_size)
  csNum = round(censoringRate*sample_size)
  ind1 = sample(1:sample_size,sample_size)
  ind2 = ind1[1:csNum]
  cs[ind2] = 1
  t[ind2] = runif(csNum,0,0.8)*t[ind2]
  y = cbind(t,1 - cs)
  colnames(y) = c("time", "status")
  
  #Estimation
  cvfit = cv.glmnet(X,y,family = "cox",type.measure = "C",nfolds = 10,nlambda=100)
  Lmd = cvfit$lambda.min
  fit = glmnet(X,y,family = "cox",lambda = Lmd)
  Eb1 = fit$beta
  
  Eb2 = cox_GAGA(X,y,ratio=2,itrNum=20,fdiag=TRUE)
 
  err1 = norm(Eb1-beta_true,type="2")
  err2 = norm(Eb2-beta_true,type="2")

  ERR_glmnet[ii] = err1
  ERR_GAGA[ii] = err2

  ACC_glmnet[ii] = calculate.w.accuracy(as.character(Eb1!=0),as.character(beta_true!=0))
  ACC_GAGA[ii] = calculate.w.accuracy(as.character(Eb2!=0),as.character(beta_true!=0))
  
  #Prediction#######################################################################################################
  #Generate test samples
  X = R1*rmvnorm(n=test_size, mean=rep(0, nrow(cov_mat)), sigma=cov_mat)
  z = X%*%beta_true
  u = runif(test_size,0,1)
  t = ((-log(1-u)/(3*exp(z)))*100)^(0.1)
  cs = rep(0,test_size)
  csNum = round(censoringRate*test_size)
  ind1 = sample(1:test_size,test_size)
  ind2 = ind1[1:csNum]
  cs[ind2] = 1
  t[ind2] = runif(csNum,0,0.8)*t[ind2]
  y = cbind(t,1 - cs)
  colnames(y) = c("time", "status")
  
  pred1 = predict(fit, newx = X)
  Cidx_glmnet[ii] = apply(pred1, 2, Cindex, y=y)
  tmpfit = fit
  tmpfit$beta = Matrix(Eb2,sparse = TRUE)
  pred2 = predict(tmpfit, newx = X)
  Cidx_GAGA[ii] = apply(pred2, 2, Cindex, y=y)
  
  

}

mean_ERR_glmnet = mean(ERR_glmnet)
mean_ERR_GAGA = mean(ERR_GAGA)
mean_ACC_glmnet = mean(ACC_glmnet)
mean_ACC_GAGA = mean(ACC_GAGA)
mean_Cidx_glmnet = mean(Cidx_glmnet)
mean_Cidx_GAGA = mean(Cidx_GAGA)

#Plotting
library(ggplot2)

ERR=c(ERR_glmnet,ERR_GAGA)
ACC=c(ACC_glmnet,ACC_GAGA)
CIDX = c(Cidx_glmnet,Cidx_GAGA)
Algorithms=factor(c(rep('glmnet_COX',Num),rep('GAGA_COX',Num)),
                  levels=c('glmnet_COX','GAGA_COX'))
ERR_ACC_CIDX=data.frame(ERR,ACC,CIDX,Algorithms)

g1=ggplot(ERR_ACC_CIDX, aes(x=Algorithms, y=ERR,fill=Algorithms))+ylab("ERR") + geom_boxplot()#ERR box
g1

g2=ggplot(ERR_ACC_CIDX, aes(x=Algorithms, y=ACC,fill=Algorithms))+xlab("Algorithms")+ylab("ACC") + geom_boxplot()#ACC box
g2

g3=ggplot(ERR_ACC_CIDX, aes(x=Algorithms, y=CIDX,fill=Algorithms))+ylab("C-index") + geom_boxplot()#C-index box
g3





# library("R.matlab")
# #filename ="F:/文档/GAGA/拟牛顿法求解logistic回归20190724/Example100.mat"
# writeMat(filename, X = X, y = y, Eb1 = Eb1, Ebb=Eb2, beta_true = beta_true)
# 



# library(survival)
# data(CoxExample)
# y = CoxExample$y
# X = CoxExample$x
# y = y[1:400,]
# X = X[1:400,]
# cvfit = cv.glmnet(X,y,family = "cox",type.measure = "C", nfolds = 10,nlambda=100)
# Lmd = cvfit$lambda.min
# Eb1 = coef(cvfit,s = Lmd)
# fit = glmnet(X,y,family = "cox",lambda = Lmd)
# Eb5 = cox_GAGA(X,y,ratio=1,itrNum=20,lamda_0=5)
# y = CoxExample$y
# X = CoxExample$x
# y = y[401:1000,]
# X = X[401:1000,]
# pred = predict(fit, newx = X)
# cidx = apply(pred, 2, Cindex, y=y)
# tmpfit = fit
# tmpfit$beta = Matrix(Eb5,sparse = TRUE)
# pred2 = predict(tmpfit, newx = X)
# cidx2 = apply(pred2, 2, Cindex, y=y)


# library("R.matlab")
# filename ="I:/王晓飞/GAGA/生存模型/tmp2.mat"
# mMAT = readMat(filename)
# 
# t = mMAT$t
# cs = 1 - mMAT$cs
# X = mMAT$X
# y = cbind(t,cs)
# colnames(y) = c("time", "status")
# beta_true = mMAT$beta
# Eb2 = mMAT$Eb2
# 
# Eb5 = cox_GAGA(X,y,ratio=2,itrNum=20)
# err_GAGA2 = norm(Eb5-beta_true,type="2")
# ACC_GAGA2 = calculate.w.accuracy(as.character(Eb5!=0),as.character(beta_true!=0))
# 
# 
# 
# cvfit = cv.glmnet(X,y,family = "cox",type.measure = "C",nfolds = 10,nlambda=100)
# Lmd = cvfit$lambda.min
# fitx = glmnet(X,y,family = "cox",lambda = Lmd)
# Eb4 = fitx$beta
# 
# err_GAGA = norm(Eb2-beta_true,type="2")
# err_glmnet = norm(Eb4-beta_true,type="2")
# 
# ACC_glmnet = calculate.w.accuracy(as.character(Eb4!=0),as.character(beta_true!=0))
# ACC_GAGA = calculate.w.accuracy(as.character(Eb2!=0),as.character(beta_true!=0))

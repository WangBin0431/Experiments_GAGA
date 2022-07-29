# Demo of Gaga dealing with Cox model
library(glmnet)
library(ncvreg)
library(mvtnorm)
library(GAGA)
library(survival)
set.seed(1234)

Nlambda=100
p_size = 50
sample_size = 500
test_size = 1000

R1 = 3
R2 = 1
rate = 0.5 #Proportion of value zero in beta
censoringRate = 0.2 #Proportion of censoring data in observation data
Num = 10 # Total number of experiments

ERR_glmnet = NULL
ACC_glmnet = NULL
ERR_MCP = NULL
ACC_MCP = NULL
ERR_SCAD = NULL
ACC_SCAD = NULL
ERR_GAGA = NULL
ACC_GAGA = NULL
Cidx_glmnet = NULL
Cidx_MCP = NULL
Cidx_SCAD = NULL
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
  Eb0 = fit$beta
  
  y_surv=Surv(y[,"time"],y[,"status"])
  
  cvfit_MCP <- cv.ncvsurv(X, y_surv, penalty = "MCP",nfolds=10,nlambda = 100)
  Eb1 = cvfit_MCP$fit$beta[,cvfit_MCP$min]

  cvfit_SCAD <- cv.ncvsurv(X, y_surv, penalty = "SCAD",nfolds=10,nlambda = 100)
  Eb2 = cvfit_SCAD$fit$beta[,cvfit_SCAD$min]
 
  fit_gaga = GAGA(X,y,alpha=3,family = "cox", itrNum=20,fdiag=TRUE)
  Eb3 = fit_gaga$beta

  ERR_glmnet[ii] = norm(Eb0-beta_true,type="2")
  ERR_MCP[ii] = norm(Eb1-beta_true,type="2")
  ERR_SCAD[ii] = norm(Eb2-beta_true,type="2")
  ERR_GAGA[ii] = norm(Eb3-beta_true,type="2")

  ACC_glmnet[ii] = cal.w.acc(as.character(Eb0!=0),as.character(beta_true!=0))
  ACC_MCP[ii] = cal.w.acc(as.character(Eb1!=0),as.character(beta_true!=0))
  ACC_SCAD[ii] = cal.w.acc(as.character(Eb2!=0),as.character(beta_true!=0))
  ACC_GAGA[ii] = cal.w.acc(as.character(Eb3!=0),as.character(beta_true!=0))
  
  #Prediction
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
  
  pred0 = predict(fit, newx = X)
  Cidx_glmnet[ii] = apply(pred0, 2, Cindex, y=y)
  
  pred1 = predict(cvfit_MCP, X,lambda=cvfit_MCP$lambda.min)
  Cidx_MCP[ii] = apply(as.matrix(pred1), 2, Cindex, y=y)
  
  pred2 = predict(cvfit_SCAD, X,lambda=cvfit_MCP$lambda.min)
  Cidx_SCAD[ii] = apply(as.matrix(pred2), 2, Cindex, y=y)

  pred3 = predict.GAGA(fit_gaga, X)
  Cidx_GAGA[ii] = cal.cindex(pred3,y)
}

mean_ERR_glmnet = mean(ERR_glmnet)
mean_ERR_MCP = mean(ERR_MCP)
mean_ERR_SCAD = mean(ERR_SCAD)
mean_ERR_GAGA = mean(ERR_GAGA)
mean_ACC_glmnet = mean(ACC_glmnet)
mean_ACC_MCP = mean(ACC_MCP)
mean_ACC_SCAD = mean(ACC_SCAD)
mean_ACC_GAGA = mean(ACC_GAGA)
mean_Cidx_glmnet = mean(Cidx_glmnet)
mean_Cidx_MCP = mean(Cidx_MCP)
mean_Cidx_SCAD = mean(Cidx_SCAD)
mean_Cidx_GAGA = mean(Cidx_GAGA)

#Plotting
library(ggplot2)

ERR=c(ERR_glmnet,ERR_MCP,ERR_SCAD,ERR_GAGA)
ACC=c(ACC_glmnet,ACC_MCP,ACC_SCAD,ACC_GAGA)
CIDX = c(Cidx_glmnet,Cidx_GAGA)
Algorithms=factor(c(rep('glmnet_COX',Num),rep('MCP_COX',Num),rep('SCAD_COX',Num),rep('GAGA_COX',Num)),
                  levels=c('glmnet_COX','MCP_COX','SCAD_COX','GAGA_COX'))
ERR_ACC_CIDX=data.frame(ERR,ACC,CIDX,Algorithms)

g1=ggplot(ERR_ACC_CIDX, aes(x=Algorithms, y=ERR,fill=Algorithms))+ylab("ERR") + geom_boxplot()#ERR box
g1

g2=ggplot(ERR_ACC_CIDX, aes(x=Algorithms, y=ACC,fill=Algorithms))+xlab("Algorithms")+ylab("ACC") + geom_boxplot()#ACC box
g2

g3=ggplot(ERR_ACC_CIDX, aes(x=Algorithms, y=CIDX,fill=Algorithms))+ylab("C-index") + geom_boxplot()#C-index box
g3

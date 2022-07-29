#Demo of Gaga dealing with logistic regression
library(glmnet)
library(ncvreg)
library(mvtnorm)
library(GAGA)

set.seed(2022)
Nlambda=100
p_size = 100
sample_size = 500
test_size = 1000

rate = 0.5 #Proportion of value zero in beta
Num = 10 # Total number of experiments

R1 = 5
R2 = 3

ERR_glmnet = NULL
ACC_glmnet = NULL
PACC_glmnet = NULL
ERR_MCP = NULL
ACC_MCP = NULL
PACC_MCP = NULL
ERR_SCAD = NULL
ACC_SCAD = NULL
PACC_SCAD = NULL
ERR_GAGA = NULL
ACC_GAGA = NULL
PACC_GAGA = NULL

for(ii in 1:Num ){
  cat("iter is:",ii,"\n");
  #Set true beta
  zeroNum = round(rate*p_size)
  ind1 = sample(1:p_size,p_size)
  ind2 = ind1[1:zeroNum]
  beta_true = runif(p_size,0.2*R2,R2)
  
  beta_true[ind2] = 0
  ind3 = ind1[(zeroNum+1):p_size]
  
  beta_pos_true=ind3;
  beta_pos_false=ind2;
  
  cov_mat=matrix(1:p_size*p_size,p_size,p_size) ##covariance matrix
  for(i in 1:p_size){for(j in 1:p_size) {if(i!=j){cov_mat[i,j]=0.0}else{cov_mat[i,j]=1}}}
  
  #Generate training samples
  X = R1*rmvnorm(n=sample_size, mean=rep(0, nrow(cov_mat)), sigma=cov_mat)
  # for(i in 1:p_size){
  #   X[1:sample_size,i] = runif(sample_size,-R1,R1)
  # }
  X[1:sample_size,1]=1
  t = 1/(1+exp(-X%*%beta_true))
  tmp = runif(sample_size,0,1)
  y = rep(0,sample_size)
  y[t>tmp] = 1
  
  #Estimation
  cvfit_glm = cv.glmnet(X[,-1], y, family = "binomial", type.measure = "class", nfolds = 10, nlambda =Nlambda)
  Lmd = cvfit_glm$lambda.min
  fit_glm = glmnet(X[,-1], y, family = "binomial", lambda = Lmd)
  Eb0 = rbind(fit_glm$a0, fit_glm$beta)
  
  cvfit_MCP <-
    cv.ncvreg(
      X[,-1],
      y,
      family = "binomial",
      penalty = "MCP",
      nlambda = Nlambda,
      nfolds = 10
    )
  Eb1 = cvfit_MCP$fit$beta[,cvfit_MCP$min]
  
  cvfit_SCAD <-
    cv.ncvreg(
      X[,-1],
      y,
      family = "binomial",
      penalty = "SCAD",
      nlambda = Nlambda,
      nfolds = 10
    )
  Eb2 = cvfit_SCAD$fit$beta[,cvfit_SCAD$min]
  
  fit_gaga = GAGA(X, y, family = "binomial", alpha=1)
  Eb3 = fit_gaga$beta

  #Prediction#######################################################################################################
  #Generate test samples
  X_t = R1*rmvnorm(n=test_size, mean=rep(0, nrow(cov_mat)), sigma=cov_mat)
  X_t[1:test_size,1]=1
  t = 1/(1+exp(-X_t%*%beta_true))
  tmp = runif(test_size,0,1)
  y_t = rep(0,test_size)
  y_t[t>tmp] = 1

  Ey0 = predict(fit_glm, newx=X_t[,-1], type="class")
  Ey1 = predict(cvfit_MCP, X_t[,-1], type="class", lambda=cvfit_MCP$lambda.min)
  Ey2 = predict(cvfit_SCAD, X_t[,-1], type="class", lambda=cvfit_SCAD$lambda.min)
  Ey3 = predict.GAGA(fit_gaga,newx=X_t)
  

  ERR_glmnet[ii] = norm(Eb0-beta_true,type="2")
  ERR_MCP[ii] = norm(Eb1-beta_true,type="2")
  ERR_SCAD[ii] = norm(Eb2-beta_true,type="2")
  ERR_GAGA[ii] = norm(Eb3-beta_true,type="2")
  
  ACC_glmnet[ii] = cal.w.acc(as.character(Eb0!=0),as.character(beta_true!=0))
  ACC_MCP[ii] = cal.w.acc(as.character(Eb1!=0),as.character(beta_true!=0))
  ACC_SCAD[ii] = cal.w.acc(as.character(Eb2!=0),as.character(beta_true!=0))
  ACC_GAGA[ii] = cal.w.acc(as.character(Eb3!=0),as.character(beta_true!=0))
  
  PACC_glmnet[ii] = cal.w.acc(as.character(Ey0),as.character(y_t))
  PACC_MCP[ii] = cal.w.acc(as.character(Ey1),as.character(y_t))
  PACC_SCAD[ii] = cal.w.acc(as.character(Ey2),as.character(y_t))
  PACC_GAGA[ii] = cal.w.acc(as.character(Ey3),as.character(y_t))

}

mean_ERR_glmnet = mean(ERR_glmnet)
mean_ERR_MCP = mean(ERR_MCP)
mean_ERR_SCAD = mean(ERR_SCAD)
mean_ERR_GAGA = mean(ERR_GAGA)

mean_ACC_glmnet = mean(ACC_glmnet)
mean_ACC_MCP = mean(ACC_MCP)
mean_ACC_SCAD = mean(ACC_SCAD)
mean_ACC_GAGA = mean(ACC_GAGA)

mean_PACC_glmnet = mean(PACC_glmnet)
mean_PACC_MCP = mean(PACC_MCP)
mean_PACC_SCAD = mean(PACC_SCAD)
mean_PACC_GAGA = mean(PACC_GAGA)

#Plotting
library(ggplot2)

ERR=c(ERR_glmnet,ERR_MCP,ERR_SCAD,ERR_GAGA)
ACC=c(ACC_glmnet,ACC_MCP,ACC_SCAD,ACC_GAGA)
PACC = c(PACC_glmnet,PACC_MCP,PACC_SCAD,PACC_GAGA)

Algorithms=factor(c(rep('glmnet_logistic',Num),rep('MCP_logistic',Num),rep('SCAD_logistic',Num),rep('GAGA_logistic',Num)),
                  levels=c('glmnet_logistic','MCP_logistic','SCAD_logistic','GAGA_logistic'))
ERR_ACC_PACC=data.frame(ERR,ACC,PACC,Algorithms)

g1=ggplot(ERR_ACC_PACC, aes(x=Algorithms, y=ERR,fill=Algorithms))+ylab("ERR") + geom_boxplot()#ERR box
g1

g2=ggplot(ERR_ACC_PACC, aes(x=Algorithms, y=ACC,fill=Algorithms))+xlab("Algorithms")+ylab("ACC") + geom_boxplot()#ACC box
g2

g3=ggplot(ERR_ACC_PACC, aes(x=Algorithms, y=PACC,fill=Algorithms))+xlab("Algorithms")+ylab("PACC") + geom_boxplot()#ACC box
g3




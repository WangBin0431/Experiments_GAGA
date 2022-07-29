#Demo of Gaga dealing with Poisson regression
library(glmnet)
library(GAGA)
library(mvtnorm)
set.seed(1234)

Nlambda=100
p_size = 100
sample_size = 200
rate = 0.5 #Proportion of value zero in beta
Num = 10 # Total number of experiments


R1 = 1
R2 = 5

ERR_glmnet = NULL
ACC_glmnet = NULL
ERR_GAGA = NULL
ACC_GAGA = NULL

for(ii in 1:Num ){
  cat("iter is:",ii,"\n");
  #Set true beta
  zeroNum = round(rate*p_size)
  ind1 = sample(1:p_size,p_size)
  ind2 = ind1[1:zeroNum]
  beta_true = runif(p_size,0.2*R2,R2)
  beta_true[ind2] = 0
  
  #Generate training samples
  X = R1*rnorm(sample_size * p_size)
  
  X[sample.int(sample_size * p_size, size = 0.5 * sample_size * p_size)] = 0
  X = matrix(X, nrow = sample_size)
  X[1:sample_size,1]=1
  y = X%*%beta_true + 2*rnorm(sample_size)
  y = as.vector(y)
  
  
  #Estimation
  cvfit = cv.glmnet(X[,-1], y, family = "gaussian", type.measure = "mse",nfolds = 10, lambda.min=0.00001,nlambda =Nlambda)
  Lmd = cvfit$lambda.min
  fit1 = glmnet(X[,-1], y, family = "gaussian", lambda = Lmd)
  Eb1 = rbind(fit1$a0,fit1$beta)
  
  fit2 = GAGA(X, y,alpha=2, family = "gaussian")
  Eb2 = fit2$beta
  
  err1 = norm(Eb1-beta_true,type="2")
  err2 = norm(Eb2-beta_true,type="2")
  
  ERR_glmnet[ii] = err1
  ERR_GAGA[ii] = err2
  
  ACC_glmnet[ii] = cal.w.acc(as.character(Eb1!=0),as.character(beta_true!=0))
  ACC_GAGA[ii] = cal.w.acc(as.character(Eb2!=0),as.character(beta_true!=0))
  
}
mean_ERR_glmnet = mean(ERR_glmnet)
mean_ERR_GAGA = mean(ERR_GAGA)
mean_ACC_glmnet = mean(ACC_glmnet)
mean_ACC_GAGA = mean(ACC_GAGA)

library(ggplot2)

ERR=c(ERR_glmnet,ERR_GAGA)
ACC=c(ACC_glmnet,ACC_GAGA)
Algorithms=factor(c(rep('glmnet_LM',Num),rep('GAGA_LM',Num)),
                  levels=c('glmnet_LM','GAGA_LM'))
ERR_ACC=data.frame(ERR,ACC,Algorithms)

g1=ggplot(ERR_ACC, aes(x=Algorithms, y=ERR,fill=Algorithms))+ylab("ERR") + geom_boxplot()#ERR box
g1

g2=ggplot(ERR_ACC, aes(x=Algorithms, y=ACC,fill=Algorithms))+xlab("Algorithms")+ylab("ACC") + geom_boxplot()#ACC box
g2
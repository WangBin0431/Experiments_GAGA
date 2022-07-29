#Demo of Gaga dealing with multinomial
library(glmnet)
library(GAGA)
library(mvtnorm)
set.seed(1234)

Nlambda=100
p_size = 20
C = 3
sample_size = 500
test_size = 1000

rate = 0.5 #Proportion of value zero in beta
Num = 10 # Total number of experiments

R1 = 1
R2 = 5
lv = 0.2

ERR_glmnet = NULL
ACC_glmnet = NULL
PACC_glmnet = NULL

ERR_GAGA = NULL
ACC_GAGA = NULL
PACC_GAGA = NULL

for(ii in 1:Num ){
  cat("iter is:",ii,"\n");
  #Set true beta
  beta_true = matrix(rep(0,p_size*C),c(p_size,C))
  zeroNum = round(rate*p_size)
  for(jj in 1:C){
    ind1 = sample(1:p_size,p_size)
    ind2 = ind1[1:zeroNum]
    tmp = runif(p_size,lv*R2,R2)
    tmp[ind2] = 0
    beta_true[,jj] = tmp
  }
  
  cov_mat=matrix(1:p_size*p_size,p_size,p_size) 
  for(i in 1:p_size){for(j in 1:p_size) {if(i!=j){cov_mat[i,j]=0.0}else{cov_mat[i,j]=1}}}
  
  #Generate training samples
  X = R1*rmvnorm(n=sample_size, mean=rep(0, nrow(cov_mat)), sigma=cov_mat)
  
  # for(i in 1:p_size){
  #   X[1:sample_size,i] = runif(sample_size,-R1,R1)
  # }

  X[1:sample_size,1]=1
  z = X%*%beta_true
  t = exp(z)/(1+rowSums(exp(z)))
  t = cbind(t,1-rowSums(t))
  tt = t(apply(t,1,cumsum))
  tt = cbind(rep(0,sample_size),tt)
  
  # y = matrix(rep(0,sample_size*(C+1)),c(sample_size,C+1))
  y = rep(0,sample_size)
  for(jj in 1:sample_size){
    tmp = runif(1,0,1)
    for(kk in 1:(C+1)){
      if((tmp>tt[jj,kk])&&(tmp<=tt[jj,kk+1])){
        y[jj] = kk
        break
      }
    }
    # y[jj] = which.max(t[jj,])
  } 
 
  
  #Estimation
  cvfit = cv.glmnet(X[,-1], y, family = "multinomial", type.measure = "class",nfolds = 10, lambda.min=0.00001,nlambda =Nlambda)

  Lmd = cvfit$lambda.min

  fitx = glmnet(X[,-1], y, family = "multinomial",lambda = Lmd)
  tmp = fitx$beta
  Eb1 = matrix(rep(1,C*p_size-C),c(p_size-1,C))
  for(kk in 1:C){
    Eb1[,kk] = as.matrix(tmp[[kk]])
  }
  Eb1 = rbind(as.numeric(fitx$a0[1:C]),Eb1)
 
  fit_gaga = GAGA(X, y,alpha=2,family = "multinomial")
  Eb2 = fit_gaga$beta
  
  #Prediction
  #Generate test samples
  X_t = R1*rmvnorm(n=test_size, mean=rep(0, nrow(cov_mat)), sigma=cov_mat)
  X_t[1:test_size,1]=1
  z = X_t%*%beta_true
  t = exp(z)/(1+rowSums(exp(z)))
  t = cbind(t,1-rowSums(t))
  tt = t(apply(t,1,cumsum))
  tt = cbind(rep(0,test_size),tt)
  
  y_t = rep(0,test_size)
  #y_tt = rep(0,test_size)
  for(jj in 1:test_size){
    tmp = runif(1,0,1)
    for(kk in 1:(C+1)){
      if((tmp>tt[jj,kk])&&(tmp<=tt[jj,kk+1])){
        y_t[jj] = kk
        break
      }
    }
    # y_t[jj] = which.max(t[jj,])
  } 
  
  Ey1 = predict(fitx, newx = X_t[,-1], type = "class")

  Ey2 = predict.GAGA(fit_gaga, newx = X_t)
  
  ERR_glmnet[ii] = norm(Eb1-beta_true,type="2")
  ERR_GAGA[ii] = norm(Eb2-beta_true,type="2")
  
  ACC_glmnet[ii] = cal.w.acc(as.character(Eb1!=0),as.character(beta_true!=0))
  ACC_GAGA[ii] = cal.w.acc(as.character(Eb2!=0),as.character(beta_true!=0))
  PACC_glmnet[ii] = cal.w.acc(as.character(Ey1),as.character(y_t))
  PACC_GAGA[ii] = cal.w.acc(as.character(Ey2),as.character(y_t))

}

mean_ERR_glmnet = mean(ERR_glmnet)
mean_ERR_GAGA = mean(ERR_GAGA)
mean_ACC_glmnet = mean(ACC_glmnet)
mean_ACC_GAGA = mean(ACC_GAGA)
 
mean_PACC_glmnet = mean(PACC_glmnet)
mean_PACC_GAGA = mean(PACC_GAGA)



#Plotting
library(ggplot2)

ERR=c(ERR_glmnet,ERR_GAGA)
ACC=c(ACC_glmnet,ACC_GAGA)
PACC = c(PACC_glmnet,PACC_GAGA)

Algorithms=factor(c(rep('glmnet_multinomial',Num),rep('GAGA_multinomial',Num)),
                  levels=c('glmnet_multinomial','GAGA_multinomial'))
ERR_ACC_PACC=data.frame(ERR,ACC,PACC,Algorithms)

g1=ggplot(ERR_ACC_PACC, aes(x=Algorithms, y=ERR,fill=Algorithms))+ylab("ERR") + geom_boxplot()#ERR box
g1

g2=ggplot(ERR_ACC_PACC, aes(x=Algorithms, y=ACC,fill=Algorithms))+xlab("Algorithms")+ylab("ACC") + geom_boxplot()#ACC box
g2

g3=ggplot(ERR_ACC_PACC, aes(x=Algorithms, y=PACC,fill=Algorithms))+xlab("Algorithms")+ylab("PACC") + geom_boxplot()#ACC box
g3



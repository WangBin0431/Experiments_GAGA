library(glmnet)
library(GAGA)
library(mvtnorm)
set.seed(1234)

Nalpha=10
nfolds = 10

p_size = 30
sample_size = 300
rate = 0.8 #Proportion of value zero in beta

R1 = 3
R2 = 1
rate = 0.5 #Proportion of value zero in beta
censoringRate = 0.2 #Proportion of censoring data in observation data

#Set true beta
zeroNum = round(rate*p_size)
ind = sample(1:p_size,zeroNum)
beta_true = runif(p_size,0.0*R2,R2)
beta_true[ind] = 0

cov_mat=matrix(1:p_size*p_size,p_size,p_size) ##covariance matrix
for(i in 1:p_size){for(j in 1:p_size) {if(i!=j){cov_mat[i,j]=0.0}else{cov_mat[i,j]=1}}}
#Generate training samples
X = R1*rmvnorm(n=sample_size, mean=rep(0, nrow(cov_mat)), sigma=cov_mat)
X[,1] = 1 
t = 1/(1+exp(-X%*%beta_true))
tmp = runif(sample_size,0,1)
y = rep(0,sample_size)
y[t>tmp] = 1

alpha_list = seq(from=0.5,to=2,length=Nalpha)
foldid  = sample(rep(seq(nfolds), length = sample_size))

mean_pacc_list = NULL
tmp_pacc_list = rep(0,nfolds)

mean_acc_list = NULL
tmp_acc_list = rep(0,nfolds)

mean_err_list = NULL
tmp_err_list = rep(0,nfolds)

tmpid = 1
for(alpha in alpha_list){
  cat("alpha is:",alpha,"\n");
  for(index in seq(nfolds)){
    X_train = X[foldid!=index,]
    y_train = y[foldid!=index]
    X_test = X[foldid==index,]
    y_test = y[foldid==index]
    fit_gaga = GAGA(X_train,y_train,family = "binomial",alpha = alpha)
    Ey = predict.GAGA(fit_gaga,X_test)
    tmp_pacc_list[index] = cal.w.acc(as.character(Ey),as.character(y_test))
    tmp_acc_list[index] = cal.w.acc(as.character(fit_gaga$beta!=0),as.character(beta_true!=0))
    tmp_err_list[index] = norm(fit_gaga$beta - beta_true, type = "2")/sqrt(length(beta_true))
  }
  mean_pacc_list[tmpid] = mean(tmp_pacc_list)
  mean_acc_list[tmpid] = mean(tmp_acc_list)
  mean_err_list[tmpid] = mean(tmp_err_list)
  tmpid = tmpid + 1
}


library(ggplot2)
g1=ggplot(data.frame(mean_pacc_list,alpha_list),aes(alpha_list,mean_pacc_list))+ geom_point(size=3)+xlab("alpha") + geom_line()
g1

g2=ggplot(data.frame(mean_acc_list,alpha_list),aes(alpha_list,mean_acc_list))+ geom_point(size=3)+xlab("alpha") + geom_line()
g2

g3=ggplot(data.frame(mean_err_list,alpha_list),aes(alpha_list,mean_err_list))+ geom_point(size=3)+xlab("alpha") + geom_line()
g3




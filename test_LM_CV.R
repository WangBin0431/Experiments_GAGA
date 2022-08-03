library(glmnet)
library(GAGA)
library(mvtnorm)
set.seed(1234)

Nalpha=20
nfolds = 10

p_size = 100
sample_size = 500
rate = 0.5 #Proportion of value zero in beta

R1 = 1
R2 = 3

#Set true beta
zeroNum = round(rate*p_size)
ind = sample(1:p_size,zeroNum)
beta_true = runif(p_size,0.0*R2,R2)
beta_true[ind] = 0

cov_mat=matrix(1:p_size*p_size,p_size,p_size) ##covariance matrix
for(i in 1:p_size){for(j in 1:p_size) {if(i!=j){cov_mat[i,j]=0.0}else{cov_mat[i,j]=1}}}
#Generate training samples
X = R1*rmvnorm(n=sample_size, mean=rep(0, nrow(cov_mat)), sigma=cov_mat)
y = X%*%beta_true + rnorm(sample_size,0,sd=2)


alpha_list = seq(from=0.5,to=10,length=Nalpha)
foldid  = sample(rep(seq(nfolds), length = sample_size))

mean_perr_list = NULL
tmp_perr_list = rep(0,nfolds)

mean_acc_list = NULL
tmp_acc_list = rep(0,nfolds)

mean_err_list = NULL
tmp_err_list = rep(0,nfolds)

tmpid = 1
for(alpha in alpha_list){
  for(index in seq(nfolds)){
    X_train = X[foldid!=index,]
    y_train = y[foldid!=index,]
    X_test = X[foldid==index,]
    y_test = y[foldid==index,]
    fit_gaga = GAGA(X_train,y_train,family = "gaussian",alpha = alpha)
    Ey = predict.GAGA(fit_gaga,X_test)
    tmp_perr_list[index] = norm(Ey - y_test, type = "2")/sqrt(length(y_test))
    tmp_acc_list[index] = cal.w.acc(as.character(fit_gaga$beta!=0),as.character(beta_true!=0))
    tmp_err_list[index] = norm(fit_gaga$beta - beta_true, type = "2")/sqrt(length(beta_true))
  }
  mean_perr_list[tmpid] = mean(tmp_perr_list)
  mean_acc_list[tmpid] = mean(tmp_acc_list)
  mean_err_list[tmpid] = mean(tmp_err_list)
  tmpid = tmpid + 1
}


library(ggplot2)
g1=ggplot(data.frame(mean_perr_list,alpha_list),aes(alpha_list,mean_perr_list))+ geom_point(size=3)+xlab("alpha") + geom_line()
g1

g2=ggplot(data.frame(mean_acc_list,alpha_list),aes(alpha_list,mean_acc_list))+ geom_point(size=3)+xlab("alpha") + geom_line()
g2

g3=ggplot(data.frame(mean_err_list,alpha_list),aes(alpha_list,mean_err_list))+ geom_point(size=3)+xlab("alpha") + geom_line()
g3




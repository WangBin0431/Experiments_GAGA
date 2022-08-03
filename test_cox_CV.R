library(glmnet)
library(GAGA)
library(mvtnorm)
set.seed(1234)

Nalpha=10
nfolds = 10

p_size = 30
sample_size = 300
rate = 0.5 #Proportion of value zero in beta

R1 = 1
R2 = 1

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


alpha_list = seq(from=0.5,to=10,length=Nalpha)
foldid  = sample(rep(seq(nfolds), length = sample_size))

mean_cidx_list = NULL
tmp_cidx_list = rep(0,nfolds)

mean_acc_list = NULL
tmp_acc_list = rep(0,nfolds)

mean_err_list = NULL
tmp_err_list = rep(0,nfolds)

tmpid = 1
for(alpha in alpha_list){
  cat("alpha is:",alpha,"\n");
  for(index in seq(nfolds)){
    X_train = X[foldid!=index,]
    y_train = y[foldid!=index,]
    X_test = X[foldid==index,]
    y_test = y[foldid==index,]
    fit_gaga = GAGA(X_train,y_train,family = "cox",alpha = alpha)
    Ey = predict.GAGA(fit_gaga,X_test)
    tmp_cidx_list[index] = cal.cindex(Ey,y_test)
    tmp_acc_list[index] = cal.w.acc(as.character(fit_gaga$beta!=0),as.character(beta_true!=0))
    tmp_err_list[index] = norm(fit_gaga$beta - beta_true, type = "2")/sqrt(length(beta_true))
  }
  mean_cidx_list[tmpid] = mean(tmp_cidx_list)
  mean_acc_list[tmpid] = mean(tmp_acc_list)
  mean_err_list[tmpid] = mean(tmp_err_list)
  tmpid = tmpid + 1
}


library(ggplot2)
g1=ggplot(data.frame(mean_cidx_list,alpha_list),aes(alpha_list,mean_cidx_list))+ geom_point(size=3)+xlab("alpha") + geom_line()
g1

g2=ggplot(data.frame(mean_acc_list,alpha_list),aes(alpha_list,mean_acc_list))+ geom_point(size=3)+xlab("alpha") + geom_line()
g2

g3=ggplot(data.frame(mean_err_list,alpha_list),aes(alpha_list,mean_err_list))+ geom_point(size=3)+xlab("alpha") + geom_line()
g3




library(glmnet)
library(GAGA)
library(mvtnorm)
set.seed(1234)

Nalpha=20
nfolds = 10


p_size = 20
C = 3
sample_size = 500
rate = 0.5 #Proportion of value zero in beta

R1 = 1
R2 = 5
lv = 0.2

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

alpha_list = seq(from=0.5,to=3,length=Nalpha)
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
    fit_gaga = GAGA(X_train,y_train,family = "multinomial",alpha = alpha)
    Ey = predict.GAGA(fit_gaga,X_test)
    tmp_pacc_list[index] = cal.w.acc(as.character(Ey),as.character(y_test))
    tmp_acc_list[index] = cal.w.acc(as.character(fit_gaga$beta!=0),as.character(beta_true!=0))
    tmp_err_list[index] = norm(fit_gaga$beta - beta_true, type = "f")/sqrt(length(beta_true))
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




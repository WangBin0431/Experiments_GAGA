library(glmnet)
library(GAGA)
set.seed(2022)
data(SparseExample)

X = SparseExample$x
y = SparseExample$y

Num = 100

X = cbind(1,X)
N = nrow(X)

train_size = round(N*0.8)
test_size = N - train_size

perr_glmnet = NULL
perr_GAGA = NULL

for(ii in 1:Num){
  #cat("iter is:",ii,"\n");
  ind = sample(1:N,N)
  
  X_train = X[ind[1:train_size],]
  X_test = X[-ind[1:train_size],]
  y_train = y[ind[1:train_size]]
  y_test = y[-ind[1:train_size]]
  
  #Estimation
  cvfit = cv.glmnet(X_train[,-1],y_train,family = "gaussian",type.measure = "mse",nfolds = 10,nlambda=100)
  Lmd = cvfit$lambda.min
  fit1 = glmnet(X_train[,-1],y_train,family = "gaussian",lambda = Lmd)
  Eb1 = fit1$beta
  
  fit2 = GAGA(as.matrix(X_train),y_train,alpha=1,itrNum=50,lamda_0 = 0.1,family = "gaussian")
  Eb2 = fit2$beta
  
  Ey1 = predict(fit1,newx=X_test[,-1])
  Ey2 = predict.GAGA(fit2,newx=as.matrix(X_test))
  
  perr_glmnet[ii] = norm(Ey1-y_test,type="2")
  perr_GAGA[ii] = norm(Ey2-y_test,type="2")
  
}

mean_perr_glmnet = mean(perr_glmnet)
mean_perr_GAGA = mean(perr_GAGA)

#Plotting
library(ggplot2)

PERR=c(perr_glmnet,perr_GAGA)
Algorithms=factor(c(rep('glmnet_LM',Num),rep('GAGA_LM',Num)),
                  levels=c('glmnet_LM','GAGA_LM'))

g1=ggplot(data.frame(PERR,Algorithms), aes(x=Algorithms, y=PERR,fill=Algorithms))+ylab("PERR") + geom_boxplot()#ERR box
g1



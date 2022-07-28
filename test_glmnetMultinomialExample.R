library(glmnet)
library(GAGA)
set.seed(2022)
data(MultinomialExample)

X = MultinomialExample$x
y = MultinomialExample$y
# X = scale(X, center = TRUE, scale = TRUE)

Num = 10

X = cbind(1,X)
N = nrow(X)

train_size = round(N*0.5)
test_size = N - train_size

pacc_glmnet = NULL
pacc_GAGA = NULL

for(ii in 1:Num){
  cat("iter is:",ii,"\n");
  ind = sample(1:N,N)
  
  X_train = X[ind[1:train_size],]
  X_test = X[-ind[1:train_size],]
  y_train = y[ind[1:train_size]]
  y_test = y[-ind[1:train_size]]
  
  #Estimation
  cvfit = cv.glmnet(X_train[,-1],y_train,family = "multinomial",type.measure = "class",nfolds = 10,nlambda=100)
  Lmd = cvfit$lambda.min
  fit1 = glmnet(X_train[,-1],y_train,family = "multinomial",lambda = Lmd)
  Eb1 = fit1$beta
  
  fit2 = GAGA(as.matrix(X_train),y_train,alpha=1,family = "multinomial")
  Eb2 = fit2$beta

  Ey1 = predict(fit1,newx=X_test[,-1],type = "class")
  Ey2 = predict.GAGA(fit2,newx=as.matrix(X_test))

  pacc_glmnet[ii] = cal.w.acc(as.character(Ey1),as.character(y_test))
  pacc_GAGA[ii] = cal.w.acc(as.character(Ey2),as.character(y_test))
  
}

mean_pacc_glmnet = mean(pacc_glmnet)
mean_pacc_GAGA = mean(pacc_GAGA)

#Plotting
library(ggplot2)

PACC=c(pacc_glmnet,pacc_GAGA)
Algorithms=factor(c(rep('glmnet_multinomial',Num),rep('GAGA_multinomial',Num)),
                  levels=c('glmnet_multinomial','GAGA_multinomial'))

g1=ggplot(data.frame(PACC,Algorithms), aes(x=Algorithms, y=PACC,fill=Algorithms))+ylab("PACC") + geom_boxplot()#ERR box
g1



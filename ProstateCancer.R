library (lasso2)
data(Prostate)
dim(Prostate)
library(ncvreg)
library(GAGA)

data_Pro=Prostate
X=model.matrix(lpsa~.,data_Pro)
y=Prostate$lpsa
samplesize=length(y)


set.seed(1234)
train_size=round(samplesize*0.9)
test_size=samplesize-train_size
expnum = 100


test_error_GAGA=NULL
test_error_GAGA_QR=NULL
test_error_LASSO=NULL
test_error_MCP=NULL
test_error_SCAD=NULL

for(iter in 1:expnum){
train = sample(1:samplesize,train_size)

##Adaptive LASSO
LSE=lm(y[train]~X[train,]-1)## Linear Regression to create the Adaptive Weights Vector
weight=abs(LSE$coefficients)^1# Using gamma = 1
XW=X[train,]%*%diag(weight)##

Nlambda=100##
# cvfit_ALASSO<- cv.ncvreg(XW, y[train],family="gaussian", penalty="lasso",nlambda =Nlambda, nfolds=10)
# fit_ALASSO = cvfit_ALASSO$fit
# tmp = fit_ALASSO$beta[,cvfit_ALASSO$min][-1]
# tmp = weight*tmp
# test_error_ALASSO[iter]=norm(y[-train]-X[-train,]%*%tmp,type='2')/test_size

cvfit_LASSO<- cv.ncvreg(X[train,-1], y[train],family="gaussian", penalty="lasso",nlambda =Nlambda, nfolds=10)
fit_LASSO = cvfit_LASSO$fit
tmp = fit_LASSO$beta[,cvfit_LASSO$min]
test_error_LASSO[iter]=norm(y[-train]-X[-train,]%*%tmp,type='2')/test_size


## SCAD 
cvfit_SCAD <- cv.ncvreg(X[train,-1], y[train],family="gaussian", penalty="SCAD",nlambda =Nlambda, nfolds=10)
fit_SCAD = cvfit_SCAD$fit

tmp = fit_SCAD$beta[,cvfit_SCAD$min]
test_error_SCAD[iter]=norm(y[-train]-X[-train,]%*%tmp,type='2')/test_size

## MCP
cvfit_MCP <- cv.ncvreg(X[train,-1], y[train],family="gaussian", penalty="MCP",nlambda =Nlambda, nfolds=10)
fit_MCP = cvfit_MCP$fit

tmp = fit_MCP$beta[,cvfit_MCP$min]
test_error_MCP[iter]=norm(y[-train]-X[-train,]%*%tmp,type='2')/test_size


##GaGa
mratio = 2
fit_gaga1 = GAGA(X[train,],y[train],alpha = mratio,family="gaussian")
EW = fit_gaga1$beta
test_error_GAGA[iter]=norm(y[-train]-X[-train,]%*%EW,type='2')/test_size

fit_gaga2 = GAGA(X[train,],y[train],alpha = mratio,family="gaussian",QR_flag = T)
EW2 = fit_gaga2$beta
test_error_GAGA_QR[iter]=norm(y[-train]-X[-train,]%*%EW2,type='2')/test_size

}

TEST_ERR=c(test_error_GAGA,test_error_GAGA_QR,test_error_LASSO,test_error_SCAD,test_error_MCP)
Algorithms=factor(c(rep('GAGA',expnum),rep('GAGA_QR',expnum),rep('LASSO_CV',expnum),rep('SCAD_CV',expnum),rep('MCP_CV',expnum)),
            levels=c('GAGA','GAGA_QR','LASSO_CV','SCAD_CV','MCP_CV'))

# Plotting
library(ggplot2)
TEST_ERR_Graph=data.frame(TEST_ERR,Algorithms)
g=ggplot(TEST_ERR_Graph, aes(x=Algorithms, y=TEST_ERR,fill=Algorithms))+ylab("TEST ERR") + geom_boxplot()
g






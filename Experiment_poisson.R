#Demo of Gaga dealing with Poisson regression
library(glmnet)
library(ncvreg)
library(mvtnorm)
library(GAGA)
set.seed(1234)

Nlambda=100
p_size = 50
sample_size = 500
test_size = 1000
rate = 0.5 #Proportion of value zero in beta
Num = 10 # Total number of experiments

R1 = 1/sqrt(p_size)
R2 = 5

ERR_MCP = NULL
ACC_MCP = NULL
ERR_SCAD = NULL
ACC_SCAD = NULL
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

  cov_mat=matrix(1:p_size*p_size,p_size,p_size) ##covariance matrix
  for(i in 1:p_size){for(j in 1:p_size) {if(i!=j){cov_mat[i,j]=0.0}else{cov_mat[i,j]=1}}}

  #Generate training samples
  X = R1*rmvnorm(n=sample_size, mean=rep(0, nrow(cov_mat)), sigma=cov_mat)
  X[1:sample_size,1]=1
  y = rpois(sample_size,lambda = as.vector(exp(X%*%beta_true)))
  y = as.vector(y)


  #Estimation
  cvfit = cv.glmnet(X[,-1], y, family = "poisson", type.measure = "deviance",nfolds = 10, lambda.min=0.00001,nlambda =Nlambda)
  Lmd = cvfit$lambda.min
  fitx = glmnet(X[,-1], y, family = "poisson",lambda = Lmd)
  Eb0 = rbind(fitx$a0, fitx$beta)
  
  cvfit <-
    cv.ncvreg(
      X[,-1],
      y,
      family = "poisson",
      penalty = "MCP",
      nlambda = Nlambda,
      nfolds = 10
    )
  Eb1 = cvfit$fit$beta[,cvfit$min]
  
  cvfit <-
    cv.ncvreg(
      X[,-1],
      y,
      family = "poisson",
      penalty = "SCAD",
      nlambda = Nlambda,
      nfolds = 10
    )
  Eb2 = cvfit$fit$beta[,cvfit$min]

  fit3 = GAGA(X, y, alpha = 2, family = "poisson")
  Eb3 = fit3$beta

  ERR_glmnet[ii] = norm(Eb0-beta_true,type="2")
  ERR_MCP[ii] = norm(Eb1-beta_true,type="2")
  ERR_SCAD[ii] = norm(Eb2-beta_true,type="2")
  ERR_GAGA[ii] = norm(Eb3-beta_true,type="2")

  ACC_glmnet[ii] = cal.w.acc(as.character(Eb0!=0),as.character(beta_true!=0))
  ACC_MCP[ii] = cal.w.acc(as.character(Eb1!=0),as.character(beta_true!=0))
  ACC_SCAD[ii] = cal.w.acc(as.character(Eb2!=0),as.character(beta_true!=0))
  ACC_GAGA[ii] = cal.w.acc(as.character(Eb3!=0),as.character(beta_true!=0))

}

mean_ERR_glmnet = mean(ERR_glmnet)
mean_ERR_MCP = mean(ERR_MCP)
mean_ERR_SCAD = mean(ERR_SCAD)
mean_ERR_GAGA = mean(ERR_GAGA)
mean_ACC_glmnet = mean(ACC_glmnet)
mean_ACC_MCP = mean(ACC_MCP)
mean_ACC_SCAD = mean(ACC_SCAD)
mean_ACC_GAGA = mean(ACC_GAGA)


#Plotting
library(ggplot2)

ERR=c(ERR_glmnet,ERR_MCP,ERR_SCAD,ERR_GAGA)
ACC=c(ACC_glmnet,ACC_MCP,ACC_SCAD,ACC_GAGA)
Algorithms=factor(c(rep('glmnet_poisson',Num),rep('MCP_poisson',Num),rep('SCAD_poisson',Num),rep('GAGA_poisson',Num)),
                  levels=c('glmnet_poisson','MCP_poisson','SCAD_poisson','GAGA_poisson'))
ERR_ACC=data.frame(ERR,ACC,Algorithms)

g1=ggplot(ERR_ACC, aes(x=Algorithms, y=ERR,fill=Algorithms))+ylab("ERR") + geom_boxplot()#ERR box
g1

g2=ggplot(ERR_ACC, aes(x=Algorithms, y=ACC,fill=Algorithms))+xlab("Algorithms")+ylab("ACC") + geom_boxplot()#ACC box
g2





library(ncvreg)
# data(Lung)
# X <- Lung$X
# y <- Lung$y
# 
# op <- par(mfrow=c(2,2))
# fit <- ncvsurv(X, y)
# plot(fit, main=expression(paste(gamma,"=",3)))
# fit <- ncvsurv(X, y, gamma=10)
# plot(fit, main=expression(paste(gamma,"=",10)))
# fit <- ncvsurv(X, y, gamma=1.5)
# plot(fit, main=expression(paste(gamma,"=",1.5)))
# fit <- ncvsurv(X, y, penalty="SCAD")
# plot(fit, main=expression(paste("SCAD, ",gamma,"=",3)))
# par(op)
# 
# fit <- ncvsurv(X,y)
# ll <- log(fit$lambda)
# op <- par(mfrow=c(2,1))
# plot(ll, BIC(fit), type="l", xlim=rev(range(ll)))
# lam <- fit$lambda[which.min(BIC(fit))]
# b <- coef(fit, lambda=lam)
# b[b!=0]
# plot(fit)
# abline(v=lam)
# par(op)
# 
# S <- predict(fit, X, type='survival', lambda=lam)
# plot(S, xlim=c(0,200))

data(Prostate)

cvfit <- cv.ncvreg(Prostate$X, Prostate$y)
plot(cvfit)
summary(cvfit)

fit <- cvfit$fit
plot(fit)
beta <- fit$beta[,cvfit$min]


X <- Prostate$X
y <- Prostate$y

cvfit <- cv.ncvreg(X, y, nfolds=length(y))

# Survival
data(Lung)
X <- Lung$X
y <- Lung$y

cvfit <- cv.ncvsurv(X, y)
summary(cvfit)
plot(cvfit)
plot(cvfit, type="rsq")
## ---- echo = FALSE-------------------------------------------------------
library(bayesm)
knitr::opts_chunk$set(fig.align = "center",
                      fig.height = 3.5,
                      warning = FALSE,
                      error = FALSE,
                      message = FALSE)

## ---- eval = FALSE-------------------------------------------------------
#  mydata <- list(y = y, X = X)
#  mymcmc <- list(R = 1e6, nprint = 0)
#  
#  out <- runireg(Data = mydata, Mcmc = mymcmc)

## ---- eval=FALSE---------------------------------------------------------
#  mydata2 <- list(y, X)
#  out <- runireg(Data = mydata2, Mcmc = mymcmc)

## ------------------------------------------------------------------------
set.seed(66)
R <- 2000
n <- 200
X <- cbind(rep(1,n), runif(n))
beta <- c(1,2)
sigsq <- 0.25
y <- X %*% beta + rnorm(n, sd = sqrt(sigsq))
out <- runireg(Data = list(y = y, X = X), Mcmc = list(R = R, nprint = 0))

summary(out$betadraw, tvalues = beta)

## ------------------------------------------------------------------------
data(cheese)
names(cheese) <- tolower(names(cheese))
str(cheese)

## ------------------------------------------------------------------------
options(digits=3)
cor(cheese$volume, cheese$price)
cor(cheese$volume, cheese$disp)

## ----fig.show = "hold"---------------------------------------------------
par(mfrow = c(1,2))

curve(dnorm(x,0,10), xlab = "", ylab = "", xlim = c(-30,30),
      main = expression(paste("Prior for ", beta[j])),
      col = "dodgerblue4")

nu  <- 3
ssq <- var(log(cheese$volume))
curve(nu*ssq/dchisq(x,nu), xlab = "", ylab = "", xlim = c(0,1),
      main = expression(paste("Prior for ", sigma^2)), 
      col = "darkred")

par(mfrow = c(1,1))

## ---- eval = FALSE-------------------------------------------------------
#  dat <- list(y = log(cheese$volume), X = model.matrix( ~ price + disp, data = cheese))
#  out <- runireg(Data = dat, Mcmc = list(R=1e4, nprint=1e3))

## ---- echo = FALSE-------------------------------------------------------
temp <- capture.output(
          {set.seed(1234);
           dat <- list(y = log(cheese$volume), X = model.matrix( ~ price + disp, data = cheese));
           out <- runireg(Data = dat, Mcmc = list(R=1e4, nprint=1e3))}, 
        file = NULL)

## ---- fig.show = "hold"--------------------------------------------------
B <- 1000+1 #burn in draws to discard 
R <- 10000

par(mfrow = c(1,2))
hist(out$betadraw[B:R,2], breaks = 30, 
     main = "Posterior Dist. of Price Coef.", 
     yaxt = "n", yaxs="i",
     xlab = "", ylab = "", 
     col = "dodgerblue4", border = "gray")
hist(out$sigmasqdraw[B:R], breaks = 30, 
     main = "Posterior Dist. of Sigma2", 
     yaxt = "n", yaxs="i",
     xlab = "", ylab = "", 
     col = "darkred", border = "gray")
par(mfrow = c(1,1))

## ------------------------------------------------------------------------
apply(out$betadraw[B:R,2:3], 2, mean)

## ------------------------------------------------------------------------
summary(out$betadraw)

## ------------------------------------------------------------------------
plot.bayesm.mat(out$betadraw[,2])

## ---- results = "hold"---------------------------------------------------
data(margarine)
str(margarine)
marg <- merge(margarine$choicePrice, margarine$demos, by = "hhid")

## ------------------------------------------------------------------------
y <- marg[,2]

## ------------------------------------------------------------------------
X1 <- createX(p=10, na=1, Xa=marg[,3:12], nd=NULL, Xd=NULL, base=1)
colnames(X1) <- c(names(marg[,3:11]), "price")
head(X1, n=10)

## ---- results = "hold"---------------------------------------------------
X2 <- createX(p=10, na=NULL, Xa=NULL, nd=2, Xd=as.matrix(marg[,c(13,16)]), base=1)
print(X2[1:10,1:9]); cat("\n")
print(X2[1:10,10:18])

## ---- eval = FALSE-------------------------------------------------------
#  X <- cbind(X1, X2[,10:ncol(X2)])
#  out <- rmnlIndepMetrop(Data = list(y=y, X=X, p=10),
#                         Mcmc = list(R=1e4, nprint=1e3))

## ---- echo = FALSE-------------------------------------------------------
temp <- capture.output(
          {set.seed(1234);
           X <- cbind(X1, X2[,10:ncol(X2)]);
           out <- rmnlIndepMetrop(Data = list(y=y, X=X, p=10), 
                                  Mcmc = list(R=1e4, nprint=1e3))}, 
        file = NULL)

## ------------------------------------------------------------------------
summary.bayesm.mat(out$betadraw[,10], names = "Price")

## ------------------------------------------------------------------------
plot.bayesm.mat(out$betadraw[,10], names = "Price")

## ------------------------------------------------------------------------
data(camera)
length(camera)
str(camera[[1]])
colnames(camera[[1]]$X)

## ---- eval = FALSE-------------------------------------------------------
#  data(margarine)
#  chpr <- margarine$choicePrice
#  chpr$hhid <- as.factor(chpr$hhid)
#  N <- nlevels(chpr$hhid)
#  dat <- vector(mode = "list", length = N)
#  for (i in 1:N) {
#    dat[[i]]$y <- chpr[chpr$hhid==levels(chpr$hhid)[i], "choice"]
#    dat[[i]]$X <- createX(p=10, na=1, Xa=chpr[chpr$hhid==levels(chpr$hhid)[i],3:12], nd=NULL, Xd=NULL)
#  }

## ------------------------------------------------------------------------
N <- length(camera)
dat <- matrix(NA, N*16, 2)
for (i in 1:length(camera)) {
  Ni <- length(camera[[i]]$y)
  dat[((i-1)*Ni+1):(i*Ni),1] <- i
  dat[((i-1)*Ni+1):(i*Ni),2] <- camera[[i]]$y
}
round(prop.table(table(dat[,2])), 3)

## ------------------------------------------------------------------------
round(prop.table(table(dat[,1], dat[,2])[41:50,], 1), 3)

## ------------------------------------------------------------------------
data  <- list(lgtdata = camera, p = 5)
prior <- list(ncomp = 1)
mcmc  <- list(R = 1e4, nprint = 0)

out <- rhierMnlRwMixture(Data = data, Prior = prior, Mcmc = mcmc)

## ------------------------------------------------------------------------
hist(apply(out$betadraw, 1:2, mean)[,9], col = "dodgerblue4", 
     xlab = "", ylab = "", yaxt="n", xlim = c(-4,6), breaks = 20,
     main = "Histogram of Posterior Means For Individual Wifi Coefs")


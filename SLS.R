set.seed(123)
list.of.packages <- c("foreign","SuperLearner")
new.packages <- list.of.packages[!(list.of.packages %in% installed.packages()[,"Package"])]
if(length(new.packages)) install.packages(new.packages, repos='http://cran.us.r-project.org')
library(SuperLearner)
library(foreign)
data <- read.csv("data.csv", sep=",")
attach(data)
SL.library <- c("SL.glm","SL.step","SL.glm.interaction")
n <- nrow(data)
nvar <- dim(data)[[2]]
Y <- data[,1]
A <- data[,2]
X <- data[,2:nvar]
W <- data[,3:nvar]
X1 <- X0 <- X
X1[,1] <- 1
X0[,1] <- 0
newdata <- rbind(X,X1,X0)
Q <- try(SuperLearner(Y = data[,1] ,X = X, SL.library=SL.library, family = "binomial", newX=newdata, method ="method.NNLS"), silent=TRUE)
Q <- as.data.frame(Q[[4]])
QAW <- Q[1:n,]
Q1W <- Q[((n+1):(2*n)),]
Q0W <- Q[((2*n+1):(3*n)),]
g <- suppressWarnings(SuperLearner(Y = data[,2], X = W, SL.library = SL.library, family = "binomial", method = "method.NNLS"))
ps <- g[[4]]
ps[ps<0.025] <- 0.025
ps[ps>0.975] <- 0.975
data <- cbind(data,QAW,Q1W,Q0W,ps,Y,A)
write.dta(data, "data2.dta")
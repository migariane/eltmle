*! version 1.5 Ensemble Learning Targeted Maximum Likelihodd by MALUQUE 05.APRIL.2017
***************************************************************************
**MIGUEL ANGEL LUQUE FERNANDEZ
**TMLE ALGORITHM IMPLEMENTATION IN STATA FOR BINARY OUTCOME AND TREATMENT 
**Improved AIPW with Super Learner (ensemble learning machine learning)
**Submitted to the IJE December 2016 (under review)
**IJE-2016-12-1473
**March 2017 
**This program requieres R to be installed in your computer (R-3.3)
****************************************************************************
capture program drop eltmle
program define eltmle
     syntax [varlist] [if] [pw] [, slaipw slaipwgbm slaipwbgam tmle tmlegbm tmlebgam aipw] 
	 marksample touse
	 local var `varlist' if `touse'
     local dir `c(pwd)'
	 cd "`dir'"
	 export delimited `var' using "data.csv", nolabel replace 
	 if "`slaipw'" == "" & "`slaipwgbm'" == "" & "`slaipwbgam'" == "" & "`tmlegbm'" == "" & "`tmlebgam'" == "" & "`aipw'" == ""{
		tmle `varlist'
	 }
	 else if "`tmlegbm'" == "tmlegbm" { 
		tmlegbm `varlist'
	 }
	 else if "`tmlebgam'" == "tmlebgam" {
		tmlebgam `varlist'
	 }
	 else if "`slaipw'" == "slaipw" { 
	 slaipw `varlist'
	 }
	 else if "`slaipwgbm'" == "slaipwgbm" {
		slaipwgbm `varlist'
	 }
	 else if "`slaipwbgam'" == "slaipwbgam" {
		slaipwbgam `varlist'
	 }
	  else if "`aipw'" == "aipw" {
		aipw `varlist'
	 }
end 

program tmle  
// Write R Code dependencies: foreign Surperlearner 
set more off
qui: file close _all
qui: file open rcode using SLS.R, write replace
qui: file write rcode ///
	`"set.seed(123)"' _newline ///
	`"list.of.packages <- c("foreign","SuperLearner")"' _newline ///
    `"new.packages <- list.of.packages[!(list.of.packages %in% installed.packages()[,"Package"])]"' _newline ///
    `"if(length(new.packages)) install.packages(new.packages)"' _newline ///
	`"library(SuperLearner)"' _newline ///
	`"library(foreign)"' _newline ///
	`"data <- read.csv("data.csv", sep=",")"' _newline ///
	`"attach(data)"' _newline ///
	`"SL.library <- c("SL.glm","SL.step","SL.glm.interaction")"' _newline ///
	`"n <- nrow(data)"' _newline ///
	`"nvar <- dim(data)[[2]]"' _newline ///
	`"Y <- data[,1]"' _newline ///
	`"A <- data[,2]"' _newline ///
	`"X <- data[,2:nvar]"' _newline ///
	`"W <- data[,3:nvar]"' _newline ///
	`"X1 <- X0 <- X"' _newline /// 
	`"X1[,1] <- 1"' _newline ///
	`"X0[,1] <- 0"' _newline ///
	`"newdata <- rbind(X,X1,X0)"' _newline /// 
	`"Q <- SuperLearner(Y = data[,1] ,X = X, SL.library=SL.library, family=binomial(), newX=newdata, method="method.NNLS")"' _newline ///
	`"Q <- as.data.frame(Q[[4]])"' _newline ///
	`"QAW <- Q[1:n,]"' _newline ///
	`"Q1W <- Q[((n+1):(2*n)),]"' _newline ///
	`"Q0W <- Q[((2*n+1):(3*n)),]"' _newline ///
	`"g <- SuperLearner(Y = data[,2], X = W, SL.library = SL.library, family = binomial(), method = "method.NNLS")"' _newline ///
	`"ps <- g[[4]]"' _newline ///
	`"ps[ps<0.025] <- 0.025"' _newline ///
	`"ps[ps>0.975] <- 0.975"' _newline ///
	`"data <- cbind(data,QAW,Q1W,Q0W,ps,Y,A)"' _newline ///
	`"write.dta(data, "data2.dta")"'  
qui: file close rcode

// Run R (you have to specify the path of your R executable file)
//shell "C:\Program Files\R\R-3.3.2\bin\x64\R.exe" CMD BATCH SLSTATA.R 
shell "/usr/local/bin/r" CMD BATCH SLS.R 

// Read Revised Data Back to Stata
clear
quietly: use "data2.dta", clear

// Q to logit scale
gen logQAW = log(QAW/(1 - QAW))
gen logQ1W = log(Q1W/(1 - Q1W))
gen logQ0W = log(Q0W/(1 - Q0W))

// Clever covariate HAW
gen double HAW = (A/ps) - ((1 - A)/(1 - ps))
gen double H1W = 1/ps
gen double H0W = -1/(1 - ps)

// Estimation of the substitution parameter (Epsilon)
qui: glm Y HAW, fam(binomial) offset(logQAW) noconstant
mat a= e(b)
gen epsilon = a[1,1]

// Targeted ATE, update from Q̅^0 (A,W) to Q̅^1 (A,W)
gen double  Qstar = exp(HAW*epsilon + logQAW)/(1 + exp(HAW*epsilon + logQAW))
gen double Q0star = exp(H0W*epsilon + logQ0W)/(1 + exp(H0W*epsilon + logQ0W))
gen double Q1star = exp(H1W*epsilon + logQ1W)/(1 + exp(H1W*epsilon + logQ1W))
summ Q1star Q0star ps

// Estimating the updated targeted ATE 
gen double ATE = (Q1star - Q0star)
qui: sum ATE
global ATEtmle = r(mean)

qui: sum Q1star
global Q1 = r(mean)
qui: sum Q0star
global Q0 = r(mean)
global RRtmle = $Q1/$Q0

// Statistical inference ATE and RR
// ATE
gen double IC = (HAW*(Y - QAW)) + (Q1W - Q0W) - $ATEtmle
qui: sum IC
global var = r(Var)
qui: count
global n = r(N)
global varICtmle = $var/$n
global pvalue = 2*(normalden(abs($ATEtmle/sqrt($varICtmle))))
global LCIa =  $ATEtmle -1.96*sqrt($varICtmle)
global UCIa =  $ATEtmle +1.96*sqrt($varICtmle)

// RR
global LCIr =  exp(log($RRtmle) -1.96*sqrt(($varICtmle)/log($RRtmle)))
global UCIr =  exp(log($RRtmle) +1.96*sqrt(($varICtmle)/log($RRtmle)))

di _newline
di "TMLE: Average Treatment Effect" _newline
di "ATE:" %9.4f $ATEtmle _col(5) "; SE:" %5.4f sqrt($varICtmle) _col(5) "; pvalue:" %5.4f $pvalue _col(5) "; 95%CI:(" %8.6f $LCIa ","  %8.6f $UCIa ")"

di _newline
di "TMLE: Relative Risk" _newline 
di "RR:" %9.4f $RRtmle _col(5) "; 95%CI:(" %6.4f $LCIr "," %6.4f $UCIr ")"

// Clean up
quietly: rm SLS.R
//quietly: rm SLS.Rout
quietly: rm data2.dta
quietly: rm data.csv 
end

///////////////////////////////////////

program tmlegbm 
// Write R Code dependencies: foreign Surperlearner 
set more off
qui: file close _all
qui: file open rcode using SLS.R, write replace
qui: file write rcode ///
	`"set.seed(123)"' _newline ///
	`"list.of.packages <- c("foreign","SuperLearner")"' _newline ///
    `"new.packages <- list.of.packages[!(list.of.packages %in% installed.packages()[,"Package"])]"' _newline ///
    `"if(length(new.packages)) install.packages(new.packages)"' _newline ///
	`"library(SuperLearner)"' _newline ///
	`"library(foreign)"' _newline ///
	`"data <- read.csv("data.csv", sep=",")"' _newline ///
	`"attach(data)"' _newline ///
	`"SL.library <- c("SL.glm","SL.step","SL.glm.interaction","SL.gbm")"' _newline ///
	`"n <- nrow(data)"' _newline ///
	`"nvar <- dim(data)[[2]]"' _newline ///
	`"Y <- data[,1]"' _newline ///
	`"A <- data[,2]"' _newline ///
	`"X <- data[,2:nvar]"' _newline ///
	`"W <- data[,3:nvar]"' _newline ///
	`"X1 <- X0 <- X"' _newline /// 
	`"X1[,1] <- 1"' _newline ///
	`"X0[,1] <- 0"' _newline ///
	`"newdata <- rbind(X,X1,X0)"' _newline /// 
	`"Q <- SuperLearner(Y = data[,1] ,X = X, SL.library=SL.library, family=binomial(), newX=newdata, method="method.NNLS")"' _newline ///
	`"Q <- as.data.frame(Q[[4]])"' _newline ///
	`"QAW <- Q[1:n,]"' _newline ///
	`"Q1W <- Q[((n+1):(2*n)),]"' _newline ///
	`"Q0W <- Q[((2*n+1):(3*n)),]"' _newline ///
	`"g <- SuperLearner(Y = data[,2], X = W, SL.library = SL.library, family = binomial(), method = "method.NNLS")"' _newline ///
	`"ps <- g[[4]]"' _newline ///
	`"ps[ps<0.025] <- 0.025"' _newline ///
	`"ps[ps>0.975] <- 0.975"' _newline ///
	`"data <- cbind(data,QAW,Q1W,Q0W,ps,Y,A)"' _newline ///
	`"write.dta(data, "data2.dta")"'  
qui: file close rcode
 
// Run R (you have to specify the path of your R executable file)
//shell "C:\Program Files\R\R-3.3.2\bin\x64\R.exe" CMD BATCH SLSTATA.R 
shell "/usr/local/bin/r" CMD BATCH SLS.R 

// Read Revised Data Back to Stata
clear
quietly: use "data2.dta", clear

// Q to logit scale
gen logQAW = log(QAW/(1 - QAW))
gen logQ1W = log(Q1W/(1 - Q1W))
gen logQ0W = log(Q0W/(1 - Q0W))

// Clever covariate HAW
gen double HAW = (A/ps) - ((1 - A)/(1 - ps))
gen double H1W = 1/ps
gen double H0W = -1/(1 - ps)

// Estimation of the substitution parameter (Epsilon)
qui: glm Y HAW, fam(binomial) offset(logQAW) noconstant
mat a= e(b)
gen epsilon = a[1,1]

// Targeted ATE, update from Q̅^0 (A,W) to Q̅^1 (A,W)
gen double  Qstar = exp(HAW*epsilon + logQAW)/(1 + exp(HAW*epsilon + logQAW))
gen double Q0star = exp(H0W*epsilon + logQ0W)/(1 + exp(H0W*epsilon + logQ0W))
gen double Q1star = exp(H1W*epsilon + logQ1W)/(1 + exp(H1W*epsilon + logQ1W))
summ Q1star Q0star ps

// Estimating the updated targeted ATE 
gen double ATE = (Q1star - Q0star)
qui: sum ATE
global ATEtmlegbm = r(mean)

qui: sum Q1star
global Q1 = r(mean)
qui: sum Q0star
global Q0 = r(mean)
global RRtmlegbm = $Q1/$Q0

// Statistical inference ATE and RR

// ATE
gen double IC = (HAW*(Y - QAW)) + (Q1W - Q0W) - $ATEtmlerf
qui: sum IC
global var = r(Var)
qui: count
global n = r(N)
global varICtmlegbm = $var/$n
global pvalue = 2*(normalden(abs($ATEtmlegbm / sqrt($varICtmlegbm))))
global LCIa =  $ATEtmlegbm -1.96*sqrt($varICtmlegbm)
global UCIa =  $ATEtmlegbm +1.96*sqrt($varICtmlegbm)

// RR
global LCIr =  exp(log($RRtmlegbm) -1.96*sqrt(($varICtmlegbm)/log($RRtmlegbm)))
global UCIr =  exp(log($RRtmlegbm) +1.96*sqrt(($varICtmlegbm)/log($RRtmlegbm)))

di _newline
di "TMLE + GBM: Average Treatment Effect" _newline
di "ATE:" %9.4f $ATEtmlegbm _col(5) "; SE:" %5.4f sqrt($varICtmlegbm) _col(5) "; pvalue:" %5.4f $pvalue _col(5) "; 95%CI:(" %8.6f $LCIa ","  %8.6f $UCIa ")"

di _newline
di "TMLE + GBM: Relative Risk" _newline 
di "RR:" %9.4f $RRtmlegbm _col(5) "; 95%CI:(" %6.4f $LCIr "," %6.4f $UCIr ")"

// Clean up
quietly: rm SLS.R
//quietly: rm SLS.Rout
quietly: rm data2.dta
quietly: rm data.csv 
end

///////////////////////////////////////

program tmlebgam 
// Write R Code dependencies: foreign Surperlearner 
set more off
qui: file close _all
qui: file open rcode using SLS.R, write replace
qui: file write rcode ///
	`"set.seed(123)"' _newline ///
	`"list.of.packages <- c("foreign","SuperLearner")"' _newline ///
    `"new.packages <- list.of.packages[!(list.of.packages %in% installed.packages()[,"Package"])]"' _newline ///
    `"if(length(new.packages)) install.packages(new.packages)"' _newline ///
	`"library(SuperLearner)"' _newline ///
	`"library(foreign)"' _newline ///
	`"data <- read.csv("data.csv", sep=",")"' _newline ///
	`"attach(data)"' _newline ///
	`"SL.library <- c("SL.glm","SL.step","SL.glm.interaction","SL.gam","SL.bayesglm")"' _newline ///
	`"n <- nrow(data)"' _newline ///
	`"nvar <- dim(data)[[2]]"' _newline ///
	`"Y <- data[,1]"' _newline ///
	`"A <- data[,2]"' _newline ///
	`"X <- data[,2:nvar]"' _newline ///
	`"W <- data[,3:nvar]"' _newline ///
	`"X1 <- X0 <- X"' _newline /// 
	`"X1[,1] <- 1"' _newline ///
	`"X0[,1] <- 0"' _newline ///
	`"newdata <- rbind(X,X1,X0)"' _newline /// 
	`"Q <- SuperLearner(Y = data[,1] ,X = X, SL.library=SL.library, family=binomial(), newX=newdata, method="method.NNLS")"' _newline ///
	`"Q <- as.data.frame(Q[[4]])"' _newline ///
	`"QAW <- Q[1:n,]"' _newline ///
	`"Q1W <- Q[((n+1):(2*n)),]"' _newline ///
	`"Q0W <- Q[((2*n+1):(3*n)),]"' _newline ///
	`"g <- SuperLearner(Y = data[,2], X = W, SL.library = SL.library, family = binomial(), method = "method.NNLS")"' _newline ///
	`"ps <- g[[4]]"' _newline ///
	`"ps[ps<0.025] <- 0.025"' _newline ///
	`"ps[ps>0.975] <- 0.975"' _newline ///
	`"data <- cbind(data,QAW,Q1W,Q0W,ps,Y,A)"' _newline ///
	`"write.dta(data, "data2.dta")"'  
qui: file close rcode
 
// Run R (you have to specify the path of your R executable file)
//shell "C:\Program Files\R\R-3.3.2\bin\x64\R.exe" CMD BATCH SLSTATA.R 
shell "/usr/local/bin/r" CMD BATCH SLS.R 

// Read Revised Data Back to Stata
clear
quietly: use "data2.dta", clear

// Q to logit scale
gen logQAW = log(QAW/(1 - QAW))
gen logQ1W = log(Q1W/(1 - Q1W))
gen logQ0W = log(Q0W/(1 - Q0W))

// Clever covariate HAW
gen double HAW = (A/ps) - ((1 - A)/(1 - ps))
gen double H1W = 1/ps
gen double H0W = -1/(1 - ps)

// Estimation of the substitution parameter (Epsilon)
qui: glm Y HAW, fam(binomial) offset(logQAW) noconstant
mat a= e(b)
gen epsilon = a[1,1]

// Targeted ATE, update from Q̅^0 (A,W) to Q̅^1 (A,W)
gen double  Qstar = exp(HAW*epsilon + logQAW)/(1 + exp(HAW*epsilon + logQAW))
gen double Q0star = exp(H0W*epsilon + logQ0W)/(1 + exp(H0W*epsilon + logQ0W))
gen double Q1star = exp(H1W*epsilon + logQ1W)/(1 + exp(H1W*epsilon + logQ1W))
summ Q1star Q0star ps

// Estimating the updated targeted ATE 
gen double ATE = (Q1star - Q0star)
qui: sum ATE
global ATEtmlebg = r(mean)

qui: sum Q1star
global Q1 = r(mean)
qui: sum Q0star
global Q0 = r(mean)
global RRtmlebg = $Q1/$Q0

// Statistical inference ATE and RR
// ATE
gen double IC = (HAW*(Y - QAW)) + (Q1W - Q0W) - $ATEtmlebg
qui: sum IC
global var = r(Var)
qui: count
global n = r(N)
global varICtmlebg = $var/$n
global pvalue = 2*(normalden(abs($ATEtmlebg/sqrt($varICtmlebg))))
global LCIa =  $ATEtmlebg -1.96*sqrt($varICtmlebg)
global UCIa =  $ATEtmlebg +1.96*sqrt($varICtmlebg)

// RR
global LCIr =  exp(log($RRtmlebg) -1.96*sqrt(($varICtmlebg)/log($RRtmlebg)))
global UCIr =  exp(log($RRtmlebg) +1.96*sqrt(($varICtmlebg)/log($RRtmlebg)))

di _newline
di "TMLE + Bayes GLM and GAM: Average Treatment Effect" _newline
di "ATE:" %9.4f $ATEtmlebg _col(5) "; SE:" %5.4f sqrt($varICtmlebg) _col(5) "; pvalue:" %5.4f $pvalue _col(5) "; 95%CI:(" %8.6f $LCIa ","  %8.6f $UCIa ")"

di _newline
di "TMLE + Bayes GLM and GAM: Relative Risk" _newline 
di "RR:" %9.4f $RRtmlebg _col(5) "; 95%CI:(" %6.4f $LCIr "," %6.4f $UCIr ")"

// Clean up
quietly: rm SLS.R
//quietly: rm SLS.Rout
quietly: rm data2.dta
quietly: rm data.csv 
end

////////////////////////////////////

program slaipw
// Write R Code dependencies: foreign Surperlearner 
set more off
qui: file close _all
qui: file open rcode using SLS.R, write replace
qui: file write rcode ///
	`"set.seed(123)"' _newline ///
	`"list.of.packages <- c("foreign","SuperLearner")"' _newline ///
    `"new.packages <- list.of.packages[!(list.of.packages %in% installed.packages()[,"Package"])]"' _newline ///
    `"if(length(new.packages)) install.packages(new.packages)"' _newline ///
	`"library(SuperLearner)"' _newline ///
	`"library(foreign)"' _newline ///
	`"data <- read.csv("data.csv", sep=",")"' _newline ///
	`"attach(data)"' _newline ///
	`"SL.library <- c("SL.glm","SL.step", "SL.glm.interaction")"' _newline ///
	`"n <- nrow(data)"' _newline ///
	`"nvar <- dim(data)[[2]]"' _newline ///
	`"Y <- data[,1]"' _newline ///
	`"A <- data[,2]"' _newline ///
	`"X <- data[,2:nvar]"' _newline ///
	`"W <- data[,3:nvar]"' _newline ///
	`"X1 <- X0 <- X"' _newline /// 
	`"X1[,1] <- 1"' _newline ///
	`"X0[,1] <- 0"' _newline ///
	`"newdata <- rbind(X,X1,X0)"' _newline /// 
	`"Q <- SuperLearner(Y = data[,1] ,X = X, SL.library=SL.library, family=binomial(), newX=newdata, method="method.NNLS")"' _newline ///
	`"Q <- as.data.frame(Q[[4]])"' _newline ///
	`"QAW <- Q[1:n,]"' _newline ///
	`"Q1W <- Q[((n+1):(2*n)),]"' _newline ///
	`"Q0W <- Q[((2*n+1):(3*n)),]"' _newline ///
	`"g <- SuperLearner(Y = data[,2], X = W, SL.library = SL.library, family = binomial(), method = "method.NNLS")"' _newline ///
	`"ps <- g[[4]]"'  _newline ///
	`"ps[ps<0.025] <- 0.025"' _newline ///
	`"ps[ps>0.975] <- 0.975"' _newline ///
	`"data <- cbind(data,QAW,Q1W,Q0W,ps,Y,A)"' _newline ///
	`"write.dta(data, "data2.dta")"'  
qui: file close rcode
 
// Run R (you have to specify the path of your R executable file)
//shell "C:\Program Files\R\R-3.3.2\bin\x64\R.exe" CMD BATCH SLSTATA.R 
shell "/usr/local/bin/r" CMD BATCH SLS.R 

// Read Revised Data Back to Stata
clear
quietly: use "data2.dta", clear

// IPTW
gen double HAW = (A/ps) - ((1 - A)/(1 - ps))
gen double H1W =  1/ps
gen double H0W =  -1/(1 - ps)

// AIPTW 
gen double ATE = HAW*(Y - QAW) + (Q1W - Q0W)
qui: sum ATE
global ATEslaipw = r(mean)

// Statistical inference ATE

// ATE
gen double IC = (HAW*(Y - QAW)) + (Q1W - Q0W) - $ATEslaipw
qui: sum IC
global var = r(Var)
qui: count
global n = r(N)
global varICslaipw = $var/$n
global pvalue = 2*(normalden(abs($ATEslaipw/sqrt($varICslaipw))))
global LCIa =  $ATEslaipw -1.96*sqrt($varICslaipw)
global UCIa =  $ATEslaipw +1.96*sqrt($varICslaipw)

// Augemented Q
gen double aQ1W = Q1W+(H1W*(Y-QAW))
gen double aQ0W = Q0W+(H0W*(Y-QAW))
sum aQ1W aQ0W ps

// RR
qui: sum aQ1W
global Q1 = r(mean)
qui: sum aQ0W
global Q0 = r(mean)
global RRslaipw = $Q1/$Q0

// RR
global LCIr =  exp(log($RRslaipw) -1.96*sqrt(($varICslaipw)/log($RRslaipw)))
global UCIr =  exp(log($RRslaipw) +1.96*sqrt(($varICslaipw)/log($RRslaipw)))

di _newline
di "AIPW ensemble learning: Average Treatment Effect" _newline
di "ATE:" %9.4f $ATEslaipw _col(5) "; SE:" %5.4f sqrt($varICslaipw) _col(5) "; pvalue:" %5.4f $pvalue _col(5) "; 95%CI:(" %8.6f $LCIa ","  %8.6f $UCIa ")"

di _newline
di "AIPW ensemble learning: Relative Risk" _newline 
di "RR:" %9.4f $RRslaipw _col(5) "; 95%CI:(" %6.4f $LCIr "," %6.4f $UCIr ")"

// Clean up
quietly: rm SLS.R
//quietly: rm SLS.Rout
quietly: rm data2.dta
quietly: rm data.csv 
end

///////////////////////

program slaipwgbm
// Write R Code dependencies: foreign Surperlearner 
set more off
qui: file close _all
qui: file open rcode using SLS.R, write replace
qui: file write rcode ///
	`"set.seed(123)"' _newline ///
	`"list.of.packages <- c("foreign","SuperLearner")"' _newline ///
    `"new.packages <- list.of.packages[!(list.of.packages %in% installed.packages()[,"Package"])]"' _newline ///
    `"if(length(new.packages)) install.packages(new.packages)"' _newline ///
	`"library(SuperLearner)"' _newline ///
	`"library(foreign)"' _newline ///
	`"data <- read.csv("data.csv", sep=",")"' _newline ///
	`"attach(data)"' _newline ///
	`"SL.library <- c("SL.glm","SL.step", "SL.glm.interaction","SL.randomForest")"' _newline ///
	`"n <- nrow(data)"' _newline ///
	`"nvar <- dim(data)[[2]]"' _newline ///
	`"Y <- data[,1]"' _newline ///
	`"A <- data[,2]"' _newline ///
	`"X <- data[,2:nvar]"' _newline ///
	`"W <- data[,3:nvar]"' _newline ///
	`"X1 <- X0 <- X"' _newline /// 
	`"X1[,1] <- 1"' _newline ///
	`"X0[,1] <- 0"' _newline ///
	`"newdata <- rbind(X,X1,X0)"' _newline /// 
	`"Q <- SuperLearner(Y = data[,1] ,X = X, SL.library=SL.library, family=binomial(), newX=newdata, method="method.NNLS")"' _newline ///
	`"Q <- as.data.frame(Q[[4]])"' _newline ///
	`"QAW <- Q[1:n,]"' _newline ///
	`"Q1W <- Q[((n+1):(2*n)),]"' _newline ///
	`"Q0W <- Q[((2*n+1):(3*n)),]"' _newline ///
	`"g <- SuperLearner(Y = data[,2], X = W, SL.library = SL.library, family = binomial(), method = "method.NNLS")"' _newline ///
	`"ps <- g[[4]]"'  _newline ///
	`"ps[ps<0.025] <- 0.025"' _newline ///
	`"ps[ps>0.975] <- 0.975"' _newline ///
	`"data <- cbind(data,QAW,Q1W,Q0W,ps,Y,A)"' _newline ///
	`"write.dta(data, "data2.dta")"'  
qui: file close rcode
 
// Run R (you have to specify the path of your R executable file)
//shell "C:\Program Files\R\R-3.3.2\bin\x64\R.exe" CMD BATCH SLSTATA.R 
shell "/usr/local/bin/r" CMD BATCH SLS.R 

// Read Revised Data Back to Stata
clear
quietly: use "data2.dta", clear

// IPTW
gen double HAW = (A/ps) - ((1 - A)/(1 - ps))
gen double H1W =  1/ps
gen double H0W =  -1/(1 - ps)

// AIPTW 
gen double ATE = HAW*(Y - QAW) + (Q1W - Q0W)
qui: sum ATE
global ATEslaipwgbm = r(mean)

// Augmented Q
gen double aQ1W = Q1W+(H1W*(Y-QAW))
gen double aQ0W = Q0W+(H0W*(Y-QAW))
sum aQ1W aQ0W ps

// RR
qui: sum aQ1W
global Q1 = r(mean)
qui: sum aQ0W
global Q0 = r(mean)
global RRslaipwgbm = $Q1/$Q0

// Statistical inference ATE 

// ATE
gen double IC = (HAW*(Y - QAW)) + (Q1W - Q0W) - $ATEslaipwrf
qui: sum IC
global var = r(Var)
qui: count
global n = r(N)
global varICslaipwgbm = $var/$n
global pvalue = 2*(normalden(abs($ATEslaipwgbm/sqrt($varICslaipwgbm))))
global LCIa =  $ATEslaipwgbm -1.96*sqrt($varICslaipwgbm)
global UCIa =  $ATEslaipwgbm +1.96*sqrt($varICslaipwgbm)

// RR
global LCIr =  exp(log($RRslaipwgbm) -1.96*sqrt(($varICslaipwgbm)/log($RRslaipwgbm)))
global UCIr =  exp(log($RRslaipwgbm) +1.96*sqrt(($varICslaipwgbm)/log($RRslaipwgbm)))

di _newline
di "AIPW + GMB: Average Treatment Effect" _newline
di "ATE:" %9.4f $ATEslaipwgbm _col(5) "; SE:" %5.4f sqrt($varICslaipwgbm) _col(5) "; pvalue:" %5.4f $pvalue _col(5) "; 95%CI:(" %8.6f $LCIa ","  %8.6f $UCIa ")"

di _newline
di "AIPW + GMB: Relative Risk" _newline 
di "RR:" %9.4f $RRslaipwgbm _col(5) "; 95%CI:(" %6.4f $LCIr "," %6.4f $UCIr ")"

// Clean up
quietly: rm SLS.R
//quietly: rm SLS.Rout
quietly: rm data2.dta
quietly: rm data.csv 
end

//////////////////////////////

program slaipwbgam
// Write R Code dependencies: foreign Surperlearner 
set more off
qui: file close _all
qui: file open rcode using SLS.R, write replace
qui: file write rcode ///
	`"set.seed(123)"' _newline ///
	`"list.of.packages <- c("foreign","SuperLearner")"' _newline ///
    `"new.packages <- list.of.packages[!(list.of.packages %in% installed.packages()[,"Package"])]"' _newline ///
    `"if(length(new.packages)) install.packages(new.packages)"' _newline ///
	`"library(SuperLearner)"' _newline ///
	`"library(foreign)"' _newline ///
	`"data <- read.csv("data.csv", sep=",")"' _newline ///
	`"attach(data)"' _newline ///
	`"SL.library <- c("SL.glm","SL.step","SL.bayesglm","SL.gam","SL.glm.interaction")"' _newline ///
	`"n <- nrow(data)"' _newline ///
	`"nvar <- dim(data)[[2]]"' _newline ///
	`"Y <- data[,1]"' _newline ///
	`"A <- data[,2]"' _newline ///
	`"X <- data[,2:nvar]"' _newline ///
	`"W <- data[,3:nvar]"' _newline ///
	`"X1 <- X0 <- X"' _newline /// 
	`"X1[,1] <- 1"' _newline ///
	`"X0[,1] <- 0"' _newline ///
	`"newdata <- rbind(X,X1,X0)"' _newline /// 
	`"Q <- SuperLearner(Y = data[,1] ,X = X, SL.library=SL.library, family=binomial(), newX=newdata, method="method.NNLS")"' _newline ///
	`"Q <- as.data.frame(Q[[4]])"' _newline ///
	`"QAW <- Q[1:n,]"' _newline ///
	`"Q1W <- Q[((n+1):(2*n)),]"' _newline ///
	`"Q0W <- Q[((2*n+1):(3*n)),]"' _newline ///
	`"g <- SuperLearner(Y = data[,2], X = W, SL.library = SL.library, family = binomial(), method = "method.NNLS")"' _newline ///
	`"ps <- g[[4]]"'  _newline ///
	`"ps[ps<0.025] <- 0.025"' _newline ///
	`"ps[ps>0.975] <- 0.975"' _newline ///
	`"data <- cbind(data,QAW,Q1W,Q0W,ps,Y,A)"' _newline ///
	`"write.dta(data, "data2.dta")"'  
qui: file close rcode
 
// Run R (you have to specify the path of your R executable file)
//shell "C:\Program Files\R\R-3.3.2\bin\x64\R.exe" CMD BATCH SLSTATA.R 
shell "/usr/local/bin/r" CMD BATCH SLS.R 

// Read Revised Data Back to Stata
clear
quietly: use "data2.dta", clear

// IPTW
gen double HAW = (A/ps) - ((1 - A)/(1 - ps))
gen double H1W =  1/ps
gen double H0W =  -1/(1 - ps)

// AIPTW 
gen double ATE = HAW*(Y - QAW) + (Q1W - Q0W)
qui: sum ATE
global ATEslaipwbg = r(mean)

// Augmented Q
gen double aQ1W = Q1W+(H1W*(Y-QAW))
gen double aQ0W = Q0W+(H0W*(Y-QAW))
sum aQ1W aQ0W ps

// RR
qui: sum aQ1W
global Q1 = r(mean)
qui: sum aQ0W
global Q0 = r(mean)
global RRslaipwbg = $Q1/$Q0

// Statistical inference ATE

// ATE
gen double IC = (HAW*(Y - QAW)) + (Q1W - Q0W) - $ATEslaipwbg
qui: sum IC
global var = r(Var)
qui: count
global n = r(N)
global varICslaipwbg = $var/$n
global pvalue = 2*(normalden(abs($ATEslaipwbg/sqrt($varICslaipwbg))))
global LCIa =  $ATEslaipwbg -1.96*sqrt($varICslaipwbg)
global UCIa =  $ATEslaipwbg +1.96*sqrt($varICslaipwbg)

// RR
global LCIr = exp(log($RRslaipwbg) -1.96*sqrt(($varICslaipwbg)/log($RRslaipwbg)))
global UCIr = exp(log($RRslaipwbg) +1.96*sqrt(($varICslaipwbg)/log($RRslaipwbg)))

di _newline
di "AIPW Bayes GLM and GAM: Average Treatment Effect" _newline
di "ATE:" %9.4f $ATEslaipwbg _col(5) "; SE:" %5.4f sqrt($varICslaipwbg) _col(5) "; pvalue:" %5.4f $pvalue _col(5) "; 95%CI:(" %8.6f $LCIa ","  %8.6f $UCIa ")"

di _newline
di "AIPW Bayes GLM and GAM: Relative Risk" _newline 
di "RR:" %9.4f $RRslaipwbg _col(5) "; 95%CI:(" %6.4f $LCIr "," %6.4f $UCIr ")"

// Clean up
quietly: rm SLS.R
//quietly: rm SLS.Rout
quietly: rm data2.dta
quietly: rm data.csv 
end

program aipw
// Write R Code dependencies: foreign Surperlearner 
set more off
qui: file close _all
qui: file open rcode using SLS.R, write replace
qui: file write rcode ///
	`"set.seed(123)"' _newline ///
	`"list.of.packages <- c("foreign")"' _newline ///
    `"new.packages <- list.of.packages[!(list.of.packages %in% installed.packages()[,"Package"])]"' _newline ///
    `"if(length(new.packages)) install.packages(new.packages)"' _newline ///
	`"library(foreign)"' _newline ///
	`"data <- read.csv("data.csv", sep=",")"' _newline ///
	`"attach(data)"' _newline ///
	`"nvar <- dim(data)[[2]]"' _newline ///
	`"Y <- data[,1]"' _newline ///
	`"A <- data[,2]"' _newline ///
	`"X <- data[,2:nvar]"' _newline ///
	`"W <- data[,3:nvar]"' _newline ///
	`"nvarX <- dim(X)[[2]]"' _newline ///
	`"xnam <- paste0("w", 2:nvarX-1)"' _newline ///
	`"fmlaY <- as.formula(paste("Y ~ A +", paste(names(W), collapse= "+")))"' _newline ///
	`"fmlaA <- as.formula(paste("A ~ ", paste(names(W), collapse= "+")))"' _newline ///
	`"Q <- glm(fmlaY, family=binomial, data=data)"' _newline ///
	`"QAW <- predict(Q, type="response")"' _newline ///
	`"Q1W = predict(Q, newdata = data.frame(A = 1, W), type="response")"' _newline ///
	`"Q0W = predict(Q, newdata = data.frame(A = 0, W), type="response")"' _newline ///
	`"g <- glm(fmlaA, family=binomial, data=data)"' _newline ///
	`"ps <- predict(g, type="response")"' _newline ///
	`"ps[ps<0.025] <- 0.025"' _newline ///
	`"ps[ps>0.975] <- 0.975"' _newline ///
	`"data <- cbind(data,QAW,Q1W,Q0W,ps,Y,A)"' _newline ///
	`"write.dta(data, "data2.dta")"'  
qui: file close rcode
 
// Run R (you have to specify the path of your R executable file)
//shell "C:\Program Files\R\R-3.3.2\bin\x64\R.exe" CMD BATCH SLSTATA.R 
shell "/usr/local/bin/r" CMD BATCH SLS.R 

// Read Revised Data Back to Stata
clear
quietly: use "data2.dta", clear

// IPTW
gen double HAW = (A/ps) - ((1 - A)/(1 - ps))
gen double H1W =  1/ps
gen double H0W =  -1/(1 - ps)

// AIPTW 
gen double ATE = (HAW*(Y - QAW)) + (Q1W - Q0W)
qui: sum ATE
global ATE = r(mean)

// Augmented Q
gen double aQ1W = Q1W+(H1W*(Y-QAW))
gen double aQ0W = Q0W+(H0W*(Y-QAW))
sum aQ1W aQ0W ps

// RR
qui: sum aQ1W
global Q1 = r(mean)
qui: sum aQ0W
global Q0 = r(mean)
global RR = $Q1/$Q0

// Statistical inference ATE

// ATE
gen double IC = (HAW*(Y - QAW)) + (Q1W - Q0W) - $ATE
qui: sum IC
global var = r(Var)
qui: count
global n = r(N)
global varIC = $var/$n
global pvalue = 2*(normalden(abs($ATE/sqrt($varIC))))
global LCIa =  $ATE -1.96*sqrt($varIC)
global UCIa =  $ATE +1.96*sqrt($varIC)

// RR
global LCIr =  exp(log($RR) -1.96*sqrt(($varIC)/log($RR)))
global UCIr =  exp(log($RR) +1.96*sqrt(($varIC)/log($RR)))

di _newline
di "AIPW: Average Treatment Effect" _newline
di "ATE:" %9.4f $ATE _col(5) "; SE:" %5.4f sqrt($varIC) _col(5) "; pvalue:" %5.4f $pvalue _col(5) "; 95%CI:(" %8.6f $LCIa ","  %8.6f $UCIa ")"

di _newline
di "AIPW: Relative Risk" _newline 
di "RR:" %9.4f $RR _col(5) "; 95%CI:(" %6.4f $LCIr "," %6.4f $UCIr ")"

// Clean up
quietly: rm SLS.R
//quietly: rm SLS.Rout
quietly: rm data2.dta
quietly: rm data.csv 
end

//program drop eltmle tmle tmlebgam tmlegbm slaipw slaipwgbm slaipwbgam aipw

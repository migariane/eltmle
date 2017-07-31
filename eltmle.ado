*! version 2.2  30.Jun.2017
*! ELTMLE: Stata module for Ensemble Learning Targeted Maximum Likelihood Estimation
*! by Miguel Angel Luque-Fernandez [cre,aut]
*! Bug reports:
*! miguel-angel.luque at lshtm.ac.uk

/*
Copyright (c) 2017  <Miguel Angel Luque-Fernandez>

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in
all copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN
THE SOFTWARE.
*/

***************************************************************************
** MIGUEL ANGEL LUQUE FERNANDEZ
** mluquefe@hsph.havard.edu // miguel-angel.luque@lshtm.ac.uk
** TMLE ALGORITHM IMPLEMENTATION IN STATA FOR BINARY OR CONTINUOUS 
** OUTCOME AND BINARY TREATMENT FOR WINDOWS USERS 
** Improved AIPTW with Super Learner (Ensemble Learning)
** This program requires R to be installed in your computer
** JUNE 2017
****************************************************************************

* Improved Imfluence curve estimation for the causal relative risk
* Improved display including potential outcomes and weights to check balance
* Included estimation for continuous outcomes 
* Added additive causal effects for continuous outcome
* Fixed risk difference and causal relative risk for SLAIPW and SLAIPWBGAM 
* Checked consistency of confidence intervals for RR

capture program drop eltmle
program define eltmle
     syntax [varlist] [if] [pw] [, slaipw slaipwbgam tmle tmlebgam] 
	 version 13.2
	 marksample touse
	 local var `varlist' if `touse'
	 tokenize `var'
	 local yvar = "`1'" 
     global flag = cond(`yvar'<=1,1,0)
	 qui sum `yvar'
	 global b = `r(max)'
	 global a = `r(min)'
	 qui replace `yvar' = (`yvar' - `r(min)') / (`r(max)' - `r(min)') if `yvar'>1
     local dir `c(pwd)'
	 cd "`dir'"
	 qui export delimited `var' using "data.csv", nolabel replace 
	 if "`slaipw'" == "" & "`slaipwbgam'" == "" & "`tmlebgam'" == "" {
		tmle `varlist'
	 }
	 else if "`tmlebgam'" == "tmlebgam" {
		tmlebgam `varlist'
	 }
	 else if "`slaipw'" == "slaipw" { 
	    slaipw `varlist'
	 }
	 else if "`slaipwbgam'" == "slaipwbgam" {
		slaipwbgam `varlist'
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
        `"if(length(new.packages)) install.packages(new.packages, repos='http://cran.us.r-project.org')"' _newline ///
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
        `"Q <- try(SuperLearner(Y = data[,1] ,X = X, SL.library=SL.library, family=binomial(), newX=newdata, method="method.NNLS"), silent=TRUE)"' _newline ///
        `"Q <- as.data.frame(Q[[4]])"' _newline ///
        `"QAW <- Q[1:n,]"' _newline ///
        `"Q1W <- Q[((n+1):(2*n)),]"' _newline ///
        `"Q0W <- Q[((2*n+1):(3*n)),]"' _newline ///
        `"g <- suppressWarnings(SuperLearner(Y = data[,2], X = W, SL.library = SL.library, family = binomial(), method = "method.NNLS"))"' _newline ///
        `"ps <- g[[4]]"' _newline ///
        `"ps[ps<0.025] <- 0.025"' _newline ///
        `"ps[ps>0.975] <- 0.975"' _newline ///
        `"data <- cbind(data,QAW,Q1W,Q0W,ps,Y,A)"' _newline ///
        `"write.dta(data, "data2.dta")"'  
qui: file close rcode

// Write bacth file to find R.exe path and R version
set more off
qui: file close _all
qui: file open bat using setup.bat, write replace
qui: file write bat ///
`"@echo off"' _newline ///
`"SET PATHROOT=C:\Program Files\R\"' _newline ///
`"echo Locating path of R..."' _newline ///
`"echo."' _newline ///
`"if not exist "%PATHROOT%" goto:NO_R"' _newline ///
`"for /f "delims=" %%r in (' dir /b "%PATHROOT%R*" ') do ("' _newline ///
        `"echo Found %%r"' _newline ///
        `"echo shell "%PATHROOT%%%r\bin\x64\R.exe" CMD BATCH SLS.R > runr.do"' _newline ///
        `"echo All set!"' _newline ///  
        `"goto:DONE"' _newline ///
`")"' _newline ///
`":NO_R"' _newline ///
`"echo R is not installed in your system."' _newline ///
`"echo."' _newline ///
`"echo Download it from https://cran.r-project.org/bin/windows/base/"' _newline ///
`"echo Install it and re-run this script"' _newline ///
`":DONE"' _newline ///
`"echo."' _newline ///
`"pause"'
qui: file close bat

//Run batch
shell setup.bat 
//Run R 
do runr.do

// Read Revised Data Back to Stata
clear
quietly: use "data2.dta", clear

// Q to logit scale
gen logQAW = log(QAW / (1 - QAW))
gen logQ1W = log(Q1W / (1 - Q1W))
gen logQ0W = log(Q0W / (1 - Q0W))
 
// Clever covariate HAW
gen double HAW = (A / ps) - ((1 - A) / (1 - ps))
gen double H1W = 1 / ps
gen double H0W = -1 / (1 - ps)

// Estimation of the substitution parameter (Epsilon)
qui glm Y HAW, fam(binomial) offset(logQAW) robust noconstant
mat a= e(b)
gen double epsilon = a[1,1]

// Targeted ATE, update from Q̅^0 (A,W) to Q̅^1 (A,W)

gen double  Qstar = exp(HAW*epsilon + logQAW)/(1 + exp(HAW*epsilon + logQAW))
gen double Q0star = exp(H0W*epsilon + logQ0W)/(1 + exp(H0W*epsilon + logQ0W))
gen double Q1star = exp(H1W*epsilon + logQ1W)/(1 + exp(H1W*epsilon + logQ1W))

global cf = $b - $a

gen double POM1 = cond($flag==1,Q1star,Q1star*($cf),.)
gen double POM0 = cond($flag==1,Q0star,Q0star*($cf),.)
gen double   PO = cond($flag==1,Qstar,Qstar*($cf),.)

gen    WT = HAW
gen    PS = ps
summ   POM1 POM0 WT PS

// Estimating the updated targeted ATE binary outcome
gen double ATE = cond($flag==1,(Q1star - Q0star),(Q1star - Q0star)*($cf),.)
qui sum ATE
global ATEtmle = r(mean)

// ATE continuous outcome 
gen double ATEci = (Q1star - Q0star)
qui sum ATEci 
global ATEci = r(mean)

// Relative risk
qui sum Q1star
global Q1 = r(mean)
qui sum Q0star
global Q0 = r(mean)
global RRtmle = $Q1/$Q0

// Statistical inference ATE 
gen double IC = cond($flag==1,HAW*(Y - Qstar) + (Q1star - Q0star) - $ATEci,(HAW*(Y - Qstar) + (Q1star - Q0star) - $ATEci)*$cf,.)

qui sum IC
global var = r(Var)
qui count
global n = r(N)
global varICtmle = $var/$n

global pvalue = cond($flag==1,2*(normalden(abs($ATEci/sqrt($varICtmle)))),2*(normalden(abs($ATEtmle/sqrt($varICtmle)))),.)

global LCIa =  cond($flag==1,$ATEci -1.96*sqrt($varICtmle),$ATEtmle -1.96*sqrt($varICtmle),.)
global UCIa =  cond($flag==1,$ATEci +1.96*sqrt($varICtmle),$ATEtmle +1.96*sqrt($varICtmle),.)

// Statistical inference RR
gen double ICrr = 1/Q1star*((A/ps)*(Y - QAW) + (Q1W) - Q1star) - 1/Q0star*((1-A)/(1-ps)*(Y - QAW) + (Q0W) - Q0star)
qui sum ICrr
global varr = r(Var)
global varICrr = $varr/$n

global LCIr =  exp(log($RRtmle) - 1.96*sqrt($varICrr))
global UCIr =  exp(log($RRtmle) + 1.96*sqrt($varICrr))

// IC for Additive Risk differences
gen double ICrd = HAW*(Y - Qstar) + (Q1star - Q0star) - $ATEci
qui sum ICrd
global varrd = r(Var)
global varICrd = $varrd/$n
global LCIard =  $ATEci -1.96*sqrt($varICrd)
global UCIard =  $ATEci +1.96*sqrt($varICrd)

// Display Results 
local bin  ""ACE (Risk Differences):" %10.4f $ATEtmle _col(5) "; SE:" %10.5f sqrt($varICtmle) _col(5) "; p-value:" %7.4f $pvalue _col(5) "; 95%CI:("  %5.4f $LCIa ","   %7.4f $UCIa ")""
local cont ""ACE (Additive Effect):" %10.4f $ATEtmle _col(5) "; Estimated Variance:" %10.4f $varICtmle _col(5) "; p-value:" %7.4f $pvalue _col(5) "; 95%CI:("  %8.2f $LCIa ","  %9.2f $UCIa ")""
local contrd  ""ACE (Risk Differences):" %10.4f $ATEci _col(5) "; SE:" %10.5f sqrt($varICrd) _col(5) "; p-value:" %7.4f $pvalue _col(5) "; 95%CI:("  %5.4f $LCIard ","  %7.4f $UCIard ")""

if $flag==1 {
di _newline
di "TMLE: Average Causal Effect" _newline
di `bin'
}
else if $flag!=1{
di _newline
di `cont'
di _newline
di `contrd'
}

local rrbin ""RR:" %9.4f $RRtmle _col(5) "; 95%CI:(" %5.4f $LCIr "," %7.4f $UCIr ")""
if $flag==1 {
di _newline
di "TMLE: Causal Relative Risk" _newline 
display `rrbin'
}
else if $flag!=1{
}

drop ATEci ICrr ICrd logQAW logQ1W logQ0W HAW H1W H0W QAW Q1W Q0W Qstar Q1star Q0star ps Y A epsilon

label var POM1 "Potential Outcome Y(1)"
label var POM0 "Potential Otucome Y(0)"
label var ATE "Average Treatment Effect"
label var WT "Inverse Probability Treatment Weights"
label var IC "Variance ATE"
label var PS "Propensity Score"

// Clean up
quietly: rm SLS.R
//quietly: rm SLS.Rout
quietly: rm data2.dta
quietly: rm data.csv
quietly: rm .RData
quietly: rm runr.do
quietly: rm setup.bat 
end

///////////////////////////////////////

program tmlebgam 
// Write R Code dependencies: foreign Surperlearner 
set more off
qui: file close _all
qui: file open rcode using SLS.R, write replace
qui: file write rcode ///
	`"set.seed(123)"' _newline ///
	`"list.of.packages <- c("foreign","SuperLearner","gam","arm")"' _newline ///
    `"new.packages <- list.of.packages[!(list.of.packages %in% installed.packages()[,"Package"])]"' _newline ///
    `"if(length(new.packages)) install.packages(new.packages, repos='http://cran.us.r-project.org')"' _newline ///
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
qui file close rcode

// Write bacth file to find R.exe path and R version
set more off
qui: file close _all
qui: file open bat using setup.bat, write replace
qui: file write bat ///
`"@echo off"' _newline ///
`"SET PATHROOT=C:\Program Files\R\"' _newline ///
`"echo Locating path of R..."' _newline ///
`"echo."' _newline ///
`"if not exist "%PATHROOT%" goto:NO_R"' _newline ///
`"for /f "delims=" %%r in (' dir /b "%PATHROOT%R*" ') do ("' _newline ///
        `"echo Found %%r"' _newline ///
        `"echo shell "%PATHROOT%%%r\bin\x64\R.exe" CMD BATCH SLS.R > runr.do"' _newline ///
        `"echo All set!"' _newline ///  
        `"goto:DONE"' _newline ///
`")"' _newline ///
`":NO_R"' _newline ///
`"echo R is not installed in your system."' _newline ///
`"echo."' _newline ///
`"echo Download it from https://cran.r-project.org/bin/windows/base/"' _newline ///
`"echo Install it and re-run this script"' _newline ///
`":DONE"' _newline ///
`"echo."' _newline ///
`"pause"'
qui: file close bat

//Run batch
shell setup.bat 
//Run R 
do runr.do

// Read Revised Data Back to Stata
clear
quietly: use "data2.dta", clear

// Q to logit scale
gen logQAW = log(QAW / (1 - QAW))
gen logQ1W = log(Q1W / (1 - Q1W))
gen logQ0W = log(Q0W / (1 - Q0W))
 
// Clever covariate HAW
gen double HAW = (A / ps) - ((1 - A) / (1 - ps))
gen double H1W = 1 / ps
gen double H0W = -1 / (1 - ps)

// Estimation of the substitution parameter (Epsilon)
qui glm Y HAW, fam(binomial) offset(logQAW) robust noconstant
mat a= e(b)
gen double epsilon = a[1,1]

// Targeted ATE, update from Q̅^0 (A,W) to Q̅^1 (A,W)

gen double  Qstar = exp(HAW*epsilon + logQAW)/(1 + exp(HAW*epsilon + logQAW))
gen double Q0star = exp(H0W*epsilon + logQ0W)/(1 + exp(H0W*epsilon + logQ0W))
gen double Q1star = exp(H1W*epsilon + logQ1W)/(1 + exp(H1W*epsilon + logQ1W))

global cfbg = $b - $a

gen double POM1 = cond($flag==1,Q1star,Q1star*($cfbg),.)
gen double POM0 = cond($flag==1,Q0star,Q0star*($cfbg),.)
gen double   PO = cond($flag==1,Qstar,Qstar*($cfbg),.)

gen    WT = HAW
gen    PS = ps
summ   POM1 POM0 WT PS

// Estimating the updated targeted ATE binary outcome
gen double ATE = cond($flag==1,(Q1star - Q0star),(Q1star - Q0star)*($cfbg),.)
qui sum ATE
global ATEtmlebg = r(mean)

// ATE continuous outcome 
gen double ATEci = (Q1star - Q0star)
qui sum ATEci 
global ATEcibg = r(mean)

// Relative risk
qui sum Q1star
global Q1bg = r(mean)
qui sum Q0star
global Q0bg = r(mean)
global RRtmlebg = $Q1bg/$Q0bg

// Statistical inference ATE 
gen double IC = cond($flag==1,HAW*(Y - Qstar) + (Q1star - Q0star) - $ATEcibg,(HAW*(Y - Qstar) + (Q1star - Q0star) - $ATEcibg)*$cfbg,.)

qui sum IC
global varbg = r(Var)
qui count
global n = r(N)
global varICtmlebg = $varbg/$n

global pvaluebg = cond($flag==1,2*(normalden(abs($ATEcibg/sqrt($varICtmlebg)))),2*(normalden(abs($ATEtmlebg/sqrt($varICtmlebg)))),.)

global LCIabg =  cond($flag==1,$ATEcibg -1.96*sqrt($varICtmlebg),$ATEtmlebg -1.96*sqrt($varICtmlebg),.)
global UCIabg =  cond($flag==1,$ATEcibg +1.96*sqrt($varICtmlebg),$ATEtmlebg +1.96*sqrt($varICtmlebg),.)

// Statistical inference RR
gen double ICrr = 1/Q1star*((A/ps)*(Y - QAW) + (Q1W) - Q1star) - 1/Q0star*((1-A)/(1-ps)*(Y - QAW) + (Q0W) - Q0star)
qui sum ICrr
global varrbg = r(Var)
global varICrrbg = $varrbg/$n

global LCIrbg =  exp(log($RRtmlebg) - 1.96*sqrt($varICrrbg))
global UCIrbg =  exp(log($RRtmlebg) + 1.96*sqrt($varICrrbg))

// IC for Additive Risk differences
gen double ICrd = HAW*(Y - Qstar) + (Q1star - Q0star) - $ATEcibg
qui sum ICrd
global varrdbg = r(Var)
global varICrdbg = $varrdbg/$n
global LCIardbg =  $ATEcibg -1.96*sqrt($varICrdbg)
global UCIardbg =  $ATEcibg +1.96*sqrt($varICrdbg)

// Display Results 
local bin  ""ACE (Risk Differences):" %10.4f $ATEtmlebg _col(5) "; SE:" %10.5f sqrt($varICtmlebg) _col(5) "; p-value:" %7.4f $pvaluebg _col(5) "; 95%CI:("  %5.4f $LCIabg ","   %7.4f $UCIabg ")""
local cont ""ACE (Additive Effect):" %10.4f $ATEtmlebg _col(5) "; Estimated Variance:" %10.4f $varICtmlebg _col(5) "; p-value:" %7.4f $pvaluebg _col(5) "; 95%CI:("  %8.2f $LCIabg ","  %9.2f $UCIabg ")""
local contrd  ""ACE (Risk Differences):" %10.4f $ATEcibg _col(5) "; SE:" %10.5f sqrt($varICrdbg) _col(5) "; p-value:" %7.4f $pvaluebg _col(5) "; 95%CI:("  %5.4f $LCIardbg ","  %7.4f $UCIardbg ")""

if $flag==1 {
di _newline
di "TMLE: Average Causal Effect" _newline
di `bin'
}
else if $flag!=1{
di _newline
di `cont'
di _newline
di `contrd'
}

local rrbin ""RR:" %9.4f $RRtmlebg _col(5) "; 95%CI:(" %5.4f $LCIrbg "," %7.4f $UCIrbg ")""
if $flag==1 {
di _newline
di "TMLE: Causal Relative Risk" _newline 
display `rrbin'
}
else if $flag!=1{
}

drop ATEci ICrr ICrd logQAW logQ1W logQ0W HAW H1W H0W QAW Q1W Q0W Qstar Q1star Q0star ps Y A epsilon

label var POM1 "Potential Outcome Y(1)"
label var POM0 "Potential Otucome Y(0)"
label var ATE "Average Treatment Effect"
label var WT "Inverse Probability Treatment Weights"
labe var IC "Variance ATE"
labe var PS "Propensity Score"

// Clean up
quietly: rm SLS.R
//quietly: rm SLS.Rout
quietly: rm data2.dta
quietly: rm data.csv
quietly: rm .RData
quietly: rm runr.do
quietly: rm setup.bat 
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
    `"if(length(new.packages)) install.packages(new.packages, repos='http://cran.us.r-project.org')"' _newline ///
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
 
// Write bacth file to find R.exe path and R version
set more off
qui: file close _all
qui: file open bat using setup.bat, write replace
qui: file write bat ///
`"@echo off"' _newline ///
`"SET PATHROOT=C:\Program Files\R\"' _newline ///
`"echo Locating path of R..."' _newline ///
`"echo."' _newline ///
`"if not exist "%PATHROOT%" goto:NO_R"' _newline ///
`"for /f "delims=" %%r in (' dir /b "%PATHROOT%R*" ') do ("' _newline ///
        `"echo Found %%r"' _newline ///
        `"echo shell "%PATHROOT%%%r\bin\x64\R.exe" CMD BATCH SLS.R > runr.do"' _newline ///
        `"echo All set!"' _newline ///  
        `"goto:DONE"' _newline ///
`")"' _newline ///
`":NO_R"' _newline ///
`"echo R is not installed in your system."' _newline ///
`"echo."' _newline ///
`"echo Download it from https://cran.r-project.org/bin/windows/base/"' _newline ///
`"echo Install it and re-run this script"' _newline ///
`":DONE"' _newline ///
`"echo."' _newline ///
`"pause"'
qui: file close bat

//Run batch
shell setup.bat 
//Run R 
do runr.do

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
global varaipw = r(Var)
qui: count
global n = r(N)
global varICslaipw = $varaipw/$n
global pvalueslaipw = 2*(normalden(abs($ATEslaipw/sqrt($varICslaipw))))
global LCIaslaipw =  $ATEslaipw -1.96*sqrt($varICslaipw)
global UCIaslaipw =  $ATEslaipw +1.96*sqrt($varICslaipw)

// Augemented Q
gen double aQ1W = Q1W+(H1W*(Y-QAW))
gen double aQ0W = Q0W+(H0W*(Y-QAW))
gen  POM1 = aQ1W
gen  POM0 = aQ0W
gen    WT = HAW
gen    PS = ps
summ POM1 POM0 WT PS

// RR
qui: sum aQ1W
global Q1slaipw = r(mean)
qui: sum aQ0W
global Q0slaipw = r(mean)
global RRslaipw = $Q1slaipw/$Q0slaipw

// RR
gen double ICrr = ((A/ps)*(Y - QAW) + (Q1W) - aQ1W) - ((1-A)/(1-ps)*(Y - QAW) + (Q0W) - aQ0W)
qui sum ICrr
global varraipw = r(Var)
global varICrrslaipw = $varraipw/$n

global LCIrslaipw =  exp(log($RRslaipw) - 1.96*sqrt($varICrrslaipw))
global UCIrslaipw =  exp(log($RRslaipw) + 1.96*sqrt($varICrrslaipw))

di _newline
di "AIPW ensemble learning: Average Treatment Effect" _newline
di "ATE:" %9.4f $ATEslaipw _col(5) "; SE:" %5.4f sqrt($varICslaipw) _col(5) "; p-value:" %5.4f $pvalueslaipw _col(5) "; 95%CI:(" %5.4f $LCIaslaipw ","  %7.4f $UCIaslaipw ")"

di _newline
di "AIPW ensemble learning: Relative Risk" _newline 
di "RR:" %9.4f $RRslaipw _col(5) "; 95%CI:(" %5.4f $LCIrslaipw "," %7.4f $UCIrslaipw ")"

drop HAW H1W H0W aQ1W aQ0W ps Y A ICrr

label var POM1 "Potential Outcome Y(1)"
label var POM0 "Potential Otucome Y(0)"
label var ATE "Average Treatment Effect"
label var WT "Inverse Probability Treatment Weights"
label var IC "Variance ATE"
label var PS "Propensity Score"

// Clean up
quietly: rm SLS.R
//quietly: rm SLS.Rout
quietly: rm data2.dta
quietly: rm data.csv
quietly: rm .RData
quietly: rm runr.do
quietly: rm setup.bat 
end

///////////////////////

program slaipwbgam
// Write R Code dependencies: foreign Surperlearner 
set more off
qui: file close _all
qui: file open rcode using SLS.R, write replace
qui: file write rcode ///
	`"set.seed(123)"' _newline ///
	`"list.of.packages <- c("foreign","SuperLearner","gam","arm")"' _newline ///
    `"new.packages <- list.of.packages[!(list.of.packages %in% installed.packages()[,"Package"])]"' _newline ///
    `"if(length(new.packages)) install.packages(new.packages, repos='http://cran.us.r-project.org')"' _newline ///
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
 
// Write bacth file to find R.exe path and R version
set more off
qui: file close _all
qui: file open bat using setup.bat, write replace
qui: file write bat ///
`"@echo off"' _newline ///
`"SET PATHROOT=C:\Program Files\R\"' _newline ///
`"echo Locating path of R..."' _newline ///
`"echo."' _newline ///
`"if not exist "%PATHROOT%" goto:NO_R"' _newline ///
`"for /f "delims=" %%r in (' dir /b "%PATHROOT%R*" ') do ("' _newline ///
        `"echo Found %%r"' _newline ///
        `"echo shell "%PATHROOT%%%r\bin\x64\R.exe" CMD BATCH SLS.R > runr.do"' _newline ///
        `"echo All set!"' _newline ///  
        `"goto:DONE"' _newline ///
`")"' _newline ///
`":NO_R"' _newline ///
`"echo R is not installed in your system."' _newline ///
`"echo."' _newline ///
`"echo Download it from https://cran.r-project.org/bin/windows/base/"' _newline ///
`"echo Install it and re-run this script"' _newline ///
`":DONE"' _newline ///
`"echo."' _newline ///
`"pause"'
qui: file close bat

//Run batch
shell setup.bat 
//Run R 
do runr.do

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
gen  POM1 = aQ1W
gen  POM0 = aQ0W
gen    WT = HAW
gen    PS = ps
summ POM1 POM0 WT PS

// RR
qui: sum aQ1W
global Q1slaipwbg = r(mean)
qui: sum aQ0W
global Q0slaipwbg = r(mean)
global RRslaipwbg = $Q1slaipwbg/$Q0slaipwbg

// Statistical inference ATE

// ATE
gen double IC = (HAW*(Y - QAW)) + (Q1W - Q0W) - $ATEslaipwbg
qui: sum IC
global varbg = r(Var)
qui: count
global n = r(N)
global varICslaipwbg = $varbg/$n
global pvaluebg = 2*(normalden(abs($ATEslaipwbg/sqrt($varICslaipwbg))))
global LCIaaipwbg =  $ATEslaipwbg -1.96*sqrt($varICslaipwbg)
global UCIaaipwbg =  $ATEslaipwbg +1.96*sqrt($varICslaipwbg)

// RR
gen double ICrr = ((A/ps)*(Y - QAW) + (Q1W) - aQ1W) - ((1-A)/(1-ps)*(Y - QAW) + (Q0W) - aQ0W)
qui sum ICrr
global varrbg = r(Var)
global varICrrbg = $varrbg/$n

global LCIraipwbg =  exp(log($RRslaipwbg) - 1.96*sqrt($varICrrbg))
global UCIraipwbg =  exp(log($RRslaipwbg) + 1.96*sqrt($varICrrbg))

di _newline
di "AIPW Bayes GLM and GAM: Average Treatment Effect" _newline
di "ATE:" %9.4f $ATEslaipwbg _col(5) "; SE:" %5.4f sqrt($varICslaipwbg) _col(5) "; p-value:" %5.4f $pvaluebg _col(5) "; 95%CI:(" %5.4f $LCIaaipwbg ","  %7.4f $UCIaaipwbg ")"

di _newline
di "AIPW Bayes GLM and GAM: Relative Risk" _newline 
di "RR:" %9.4f $RRslaipwbg _col(5) "; 95%CI:(" %5.4f $LCIraipwbg "," %7.4f $UCIraipwbg ")"

drop HAW H1W H0W aQ1W aQ0W ps Y A ICrr

label var POM1 "Potential Outcome Y(1)"
label var POM0 "Potential Otucome Y(0)"
label var ATE "Average Treatment Effect"
label var WT "Inverse Probability Treatment Weights"
label var IC "Variance ATE"
label var PS "Propensity Score"

// Clean up
quietly: rm SLS.R
//quietly: rm SLS.Rout
quietly: rm data2.dta
quietly: rm data.csv
quietly: rm .RData
quietly: rm runr.do
quietly: rm setup.bat 
end
////////////////////////////////////////

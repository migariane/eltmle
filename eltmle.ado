*! version 3.1.5  29.October.2024
*! ELTMLE: Stata module for Ensemble Learning Targeted Maximum Likelihood Estimation
*! by Miguel Angel Luque-Fernandez [cre,aut]
*! and Matthew J. Smith [aut]
*! and Camille Maringe [aut]
*! Bug reports:
*! miguel-angel.luque at lshtm.ac.uk
*! matt.smith at lshtm.ac.uk
*! camille.maringe at lshtm.ac.uk


/*
Copyright (c) 2024  <Miguel Angel Luque-Fernandez> & <Matthew J. Smith> & <Camille Maringe> 

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
FITNESS FOR A PARTICULAR PURPOSE AND NON INFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN
THE SOFTWARE.
*/

***************************************************************************
** MIGUEL ANGEL LUQUE FERNANDEZ
** miguel-angel.luque at lshtm.ac.uk 
** mluquefe at ugr.es
** TMLE ALGORITHM IMPLEMENTATION IN STATA FOR BINARY OR CONTINUOUS
** OUTCOME AND BINARY TREATMENT FOR MAC and WINDOWS USERS
** This program requires R to be installed in your computer
** June 2024
****************************************************************************

* October 2018:
* Improved the output including potential outcomes and propensity score
* Included estimation for continuous outcomes
* Included marginal odds ratio
* Improved estimation of the clever covariate for both H1W and H0W
* Included Influence curve (IC) estimation for the marginal OR
* Improved IC estimation
* Update globals to locals where possible
* Just one ado file for both Mac and Windows users
* Included additive effect for continuous outcomes
* Fixed ATE 95%CI for additive risk difference

* February 2019:
* Included HAW as a sampling weight in MLE for targeted step (gain in efficiency) for the ATE

* July 2019:
* Updated as a rclass program: returning scalars for ATE, ATE 95%CI, ATE SE, CRR, MOR, CRR SEs, and MOR SEs
* Improved the output display

* November 2020:
* * Keep initial dataset

* January 2021:
* Added bal option to visually display positivity violations

* June 2021:
* Complete case (Listwise) analysis

* July 2022:
* Improved display of the results Stata like format 

* January 2023:
* Improved display of the CRR and MOR results to Stata like format

* February 2023:
* Added bal option
* Added elements option 

* March 2023:
* Updated 95%CI for LogOR to match R-TMLE implementation

* March 2024:
* Removed a bug that gave different results when using the 'bal' option

* July 2024:
* Added an option (cvtmle) to cross-validate of the initial outcome
* Added an option (cvtmleglsrf) to cross-validate the initial outcome and include random forest algorithm for prediction
* Added an option (cvfolds(#)) for the user to choose the number of folds during cross-validation
* Added an option (seed()) for the user to specify the starting seed and obtain the same results during cross-validation
* Added an option (cvres) to show the fold-specific results during cross-validation
* Updated cross-validation programs to get covariate balance table using the 'bal' option

* September 2024:
* Removed the bug that was causing the user to manually close the Windows shell. Used 'rscript' to locate where R is stored.

* October 2024:
* Seed option is no longer necessary but remains in the syntax for reproducibility.
* Added a progress bar for cross-validation.


capture program drop eltmle
program define eltmle
	syntax varlist(min=3) [if] [pw] [, tmle tmlebgam tmleglsrf bal elements cvtmle cvtmleglsrf cvtmlebgam cvfolds(int 10) seed(numlist)]
        version 13.2

	 // Drop missing values	only for variables in the varlist
		 foreach v of varlist `varlist' {
			qui drop if missing(`v')
         }
	 
	 // Confirm that exposure is binary
		local var  "`varlist'" 
		tokenize `var'
		local xvar = "`2'"
		quietly tabulate `xvar'
		if r(r)!=2  di in red "`xvar' is not binary. Please either convert this variable to be binary or use a binary variable."
		if r(r)!=2 exit(0)

		
     // Drop the elements if they have been defined already
				capture drop _merge
				capture drop foldid
				capture drop rowid
				capture drop _d1
				capture drop _d0
				capture drop _QAW
				capture drop _Q1W
				capture drop _Q0W
				capture drop _Q1star
				capture drop _Qa1star
				capture drop _Q0star
				capture drop _Qa0star
				capture drop _ATE
				capture drop _IC
				capture drop _Y
				capture drop _A
				capture drop _POM1
				capture drop _POM0
				capture drop _ps
				capture drop _cin
				capture drop _ipw
				capture drop d1A
				capture drop x1pointsa
				capture drop d0A
				capture drop x0pointsa


	// Specify whether the user is doing Cross-Validation
	if "`cvtmle'" != "" | "`cvtmleglsrf'" != "" | "`cvtmlebgam'" != "" {
		qui export delimited using "fulldata.csv", nolabel replace

		global folds = `cvfolds'
		capture drop foldid
		
		if "`seed'" != "" {
			set seed `seed'
			qui xtile foldid = uniform() , nq(`cvfolds') 
			}
		else if "`seed'" == "" {
			qui xtile foldid = uniform() , nq(`cvfolds') 
			}
		* global setseed = `seed'		// set seed `seed'
		*set seed $setseed
		*qui xtile foldid = uniform() , nq(`cvfolds') 
		//marksample touse
		global variablelist `varlist'
		local var  "`varlist'" // if `touse'
		
		* Create a row identifier
		capture drop rowid
		gen rowid = _n
		
		// Identify continuous/binary outcome
		tokenize `var'
		local yvar = "`1'"
		macro shift 
		local rest "`*'"
		//global flag = cond(`yvar'<=1 & `yvar'>=0,1,0)   // 1 if binary, 0 if continuous.
		qui tabulate `yvar'
		global flag = cond(r(r)==2,1,0)
		*di "Flag = $flag"
				
		if $flag == 0 {  
			*di "Continuous outcome"
			capture drop ytempvar
			gen ytempvar = `yvar'
			qui sum ytempvar
			global b = r(max)
			global a = r(min)
			replace ytempvar = (ytempvar - r(min)) / (r(max) - r(min))  // if `yvar'>1
			global varusedforcv "ytempvar `rest' foldid rowid"
		} 
		else if $flag == 1 {
			*di "Binary outcome"
			qui sum `yvar'
			global b = r(max)
			global a = r(min)
			global varusedforcv "`var' foldid rowid"
		}		
		
		// Directory and data sets
		local dir `c(pwd)'
		qui cd "`dir'"
		//tempfile cvdata
		qui save "cvdata.dta", replace
		qui export delimited $varusedforcv using "cvdata.csv", nolabel replace
	}
	else if "`cvtmle'" == "" & "`cvtmleglsrf'" == "" & "`cvtmlebgam'" == "" {
		qui export delimited using "fulldata.csv", nolabel replace
		marksample touse
		global variablelist `varlist'
		local var `varlist' if `touse'
		tokenize `var'
		local yvar = "`1'"
		global flag = cond(`yvar'<=1,1,0) 
		*di "Flag = $flag"
		qui sum `yvar'
		global b = `r(max)'
		global a = `r(min)'
		qui replace `yvar' = (`yvar' - `r(min)') / (`r(max)' - `r(min)') // if `yvar'>1
		local dir `c(pwd)'
		cd "`dir'"
		tempfile data
		qui save "`data'.dta", replace
		qui export delimited `var' using "data.csv", nolabel replace
	}


	// Create global macro to keep the elements of the TMLE
		if "`elements'" != "" {
			global keepvars = 1
		}
		else if "`elements'" == "" {
			global keepvars = 0
		}

	// Create global macro to show the covariate balance table 
		if "`bal'" != "" {
			global covbalancetable = 1
		}
		else if "`bal'" == "" {
			global covbalancetable = 0
		}

	// Create global macro to show fold-specific results
		if "`cvres'" != "" {
			global cvresultstable = 1
		}
		else if "`cvres'" == "" {
			global cvresultstable = 0
		}	
		
		// 	// Create global macro to show progress of CVTMLE
		// 		if "`progress'" != "" {
		// 			global progress = 1
		// 		}
		// 		else if "`cvres'" == "" {
		// 			global progress = 0
		// 		}	


	// Select the program based on the defined options
		if "`cvtmle'" == "cvtmle" & "`cvfolds'" == "" {
			cvtmle `varlist'
		}
		else if "`cvtmle'" == "cvtmle" & "`cvfolds'" != "" {
			cvtmle `varlist'
		}
		else if "`cvtmleglsrf'" == "cvtmleglsrf" & "`cvfolds'" == "" {
			if $flag == 0 {	
				cvtmleglsrfcont `varlist' 
			}
			if $flag == 1 {
				cvtmleglsrfbin `varlist' 
			}
		}
		else if "`cvtmleglsrf'" == "cvtmleglsrf" & "`cvfolds'" != "" {
			if $flag == 0 {	
				cvtmleglsrfcont `varlist' 
			}
			if $flag == 1 { 
				cvtmleglsrfbin `varlist' 
			}
		}
		else if "`cvtmlebgam'" == "cvtmlebgam" & "`cvfolds'" == "" {
			cvtmlebgam `varlist'
		}
		else if "`cvtmlebgam'" == "cvtmlebgam" & "`cvfolds'" != "" {
			cvtmlebgam `varlist'
		}
		else if "`tmlebgam'" == "" & "`tmleglsrf'" == "" & "`bal'" == ""{
			tmle `varlist'
		}
		else if "`tmle'" != "" & "`tmlebgam'" != "" {
			di as error "Both tmle and tmlebgam are specified. Please specify only tmle or tmlebgam, but not both."
		}
		else if "`tmle'" != "" & "`tmleglsrf'" != "" {
			di as error "Both tmle and tmleglsrf are specified. Please specify only tmle or tmleglsrf, but not both."
		}
		else if "`tmlebgam'" != "" & "`tmleglsrf'" != "" {
			di as error "Both tmlebgam and tmleglsrf are specified. Please specify only tmlebgam or tmleglsrf, but not both."
		}
		else if "`tmle'" == "tmle" & "`bal'" == "bal" {
			tmlebal `varlist'
			}
		else if "`tmlebgam'" == "tmlebgam" & "`bal'" == "bal" {
			tmlebgambal `varlist'
			}
		else if "`tmleglsrf'" == "tmleglsrf" & "`bal'" == "bal" {
			tmleglsrfbal `varlist'
			}
		else if "`tmlebgam'" == "tmlebgam" {
			tmlebgam `varlist'
			}
		else if "`tmleglsrf'" == "tmleglsrf" {
			tmleglsrf `varlist'
			}
		else if "`bal'" == "bal" {
			tmlebal `varlist'
			}
end





program tmle, rclass
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
        `"fulldata <- read.csv("fulldata.csv", sep=",")"' _newline ///
        `"attach(data)"' _newline ///
        `"SL.library <- c("SL.glm","SL.step","SL.glm.interaction")"' _newline ///
        `"n <- nrow(data)"' _newline ///
        `"nvar <- dim(data)[[2]]"' _newline ///
        `"Y <- data[,1]"' _newline ///
        `"A <- data[,2]"' _newline ///
        `"X <- data[,2:nvar]"' _newline ///
        `"W <- as.data.frame(data[,3:nvar])"' _newline ///
        `"X1 <- X0 <- X"' _newline ///
        `"X1[,1] <- 1"' _newline ///
        `"X0[,1] <- 0"' _newline ///
        `"newdata <- rbind(X,X1,X0)"' _newline ///
        `"Q <- try(SuperLearner(Y = data[,1] ,X = X, SL.library=SL.library, family = "binomial", newX=newdata, method ="method.NNLS"), silent=TRUE)"' _newline ///
        `"Q <- as.data.frame(Q[[4]])"' _newline ///
        `"QAW <- Q[1:n,]"' _newline ///
        `"Q1W <- Q[((n+1):(2*n)),]"' _newline ///
        `"Q0W <- Q[((2*n+1):(3*n)),]"' _newline ///
        `"g <- suppressWarnings(SuperLearner(Y = data[,2], X = W, SL.library = SL.library, family = "binomial", method = "method.NNLS"))"' _newline ///
        `"ps <- g[[4]]"' _newline ///
        `"ps[ps<0.025] <- 0.025"' _newline ///
		`"ps[ps>0.975] <- 0.975"' _newline ///
		`"data <- cbind(fulldata,QAW,Q1W,Q0W,ps,Y,A)"' _newline ///
		`"write.dta(data, "data2.dta")"'
qui: file close rcode


	if "`c(os)'" == "MacOSX" {
	qui shell "/usr/local/bin/r" CMD BATCH SLS.R
	}
	else{
	// Write batch file to find R.exe path and R version
		set more off
		qui: file close _all
		qui: file open bat using setup.bat, write replace
		qui: file write bat ///
		`"@echo off"' _newline ///
		`"SET PATHROOT=C:\Program Files\R\"' _newline ///
		`"echo Locating path of R..."' _newline ///
		`"echo."' _newline ///
		`"if not exist "%PATHROOT%" goto:NO_R"' _newline ///
		`":NO_R"' _newline ///
		`"echo R is not installed in your system."' _newline ///
		`"echo."' _newline ///
		`"echo Download it from https://cran.r-project.org/bin/windows/base/"' _newline ///
		`"echo Install it and re-run this script"' _newline ///
		`":DONE"' _newline ///
		`"echo."' _newline ///
		`"pause"'
		qui: file close bat
		
		// 		`"for /f "delims=" %%r in (' dir /b "%PATHROOT%R*" ') do ("' _newline ///
		// 				`"echo Found %%r"' _newline ///
		// 				`"echo ! "%PATHROOT%%%r\bin\x64\R.exe" CMD BATCH SLS.R > runr.do"' _newline ///
		// 				`"echo All set!"' _newline ///
		// 				`"goto:DONE"' _newline ///
		// 		`")"' _newline ///
		
		// 		//Run batch
		// 		! setup.bat
		// 		//Run R
		// 		do runr.do
		
		* Check you have 'rscript' command installed.
		qui: net install rscript, from("https://raw.githubusercontent.com/reifjulian/rscript/master") replace
		
		* Run SLS.R in Windows - 'rscript' command will search for common folders to find R.exe
		qui: rscript using SLS.R
	}

// Read Revised Data Back to Stata
	clear
	quietly: use "data2.dta", clear
	qui cap drop X__000000
	tempvar logQAW logQ1W logQ0W HAW H1W H0W eps1 eps2 eps ATE ICrr ICor

// Q to logit scale
	gen `logQAW' = log(QAW / (1 - QAW))
	gen `logQ1W' = log(Q1W / (1 - Q1W))
	gen `logQ0W' = log(Q0W / (1 - Q0W))

// Clever covariate HAW
	gen  `HAW' = (A / ps) - ((1 - A) / (1 - ps))
	gen  `H1W' = A / ps
	gen  `H0W' = (1 - A) / (1 - ps)

// Estimation of the substitution parameter (Epsilon)
	qui glm Y `H1W' `H0W', fam(binomial) offset(`logQAW') robust noconstant
	mat a= e(b)
	gen `eps1' = a[1,1]
	gen `eps2' = a[1,2]

	qui glm Y `HAW', fam(binomial) offset(`logQAW') robust noconstant
	mat a= e(b)
	gen `eps' = a[1,1]
	
	
// Targeted ATE, update from Q̅^0 (A,W) to Q̅^1 (A,W)
	gen double Qa0star = exp(`H0W'*`eps' + `logQ0W')/(1 + exp(`H0W'*`eps' + `logQ0W'))
	gen double Qa1star = exp(`H1W'*`eps' + `logQ1W')/(1 + exp(`H1W'*`eps' + `logQ1W'))

	gen double Q0star = exp(`H0W'*`eps2' + `logQ0W')/(1 + exp(`H0W'*`eps2' + `logQ0W'))
	gen double Q1star = exp(`H1W'*`eps1' + `logQ1W')/(1 + exp(`H1W'*`eps1' + `logQ1W'))

	//gen double cin = ($b - $a) 

	gen double POM1 = cond($flag == 1, Qa1star, (Qa1star * ( $b - $a )) + $a , .)
	gen double POM0 = cond($flag == 1, Qa0star, (Qa0star * ( $b - $a )) + $a , .)

	sum POM1 POM0 ps
	di as text " "

// Estimating the updated targeted ATE binary outcome
	gen double ATE = cond($flag == 1, (Qa1star - Qa0star), ((Qa1star * ( $b - $a )) + $a ) - ((Qa0star * ( $b - $a )) + $a ) , .)
	qui sum ATE
	return scalar ATEtmle = r(mean)

// Relative risk
	qui sum Q1star
	local Q1 = r(mean)
	qui sum Q0star
	local Q0 = r(mean)

// Relative risk and Odds ratio
	local RRtmle = `Q1'/`Q0'
	local logRRtmle = log(`Q1') - log(`Q0')
	local ORtmle = (`Q1' * (1 - `Q0')) / ((1 - `Q1') * `Q0')

// Statistical inference (Efficient Influence Curve)
	gen d1 = cond($flag == 1,(A * (Y - Q1star) / ps) + Q1star - `Q1',(A * (Y - Qa1star) / ps) + Qa1star - `Q1' ,.)
	gen d0 = cond($flag == 1,(1 - A) * (Y - Q0star) / (1 - ps) + Q0star - `Q0',(1 - A) * (Y - Qa0star) / (1 - ps) + Qa0star - `Q0' ,.)
	gen IC = cond($flag == 1,(d1 - d0), (d1 * ($b - $a ) + $a ) - (d0 * ($b - $a ) + $a ), .)
	qui sum IC
	return scalar ATE_SE_tmle = sqrt(r(Var)/r(N))

// Statistical inference ATE
	return scalar ATE_pvalue =  2 * (normalden(abs(return(ATEtmle) / (return(ATE_SE_tmle)))))
	return scalar ATE_LCIa   =  return(ATEtmle) - 1.96 * return(ATE_SE_tmle)
	return scalar ATE_UCIa   =  return(ATEtmle) + 1.96 * return(ATE_SE_tmle)

// Statistical inference RR
	gen `ICrr' = (1/`Q1' * d1) + ((1/`Q0') * d0)
	qui sum `ICrr'
	local varICrr = r(Var)/r(N)

	local LCIrr =  exp(`logRRtmle' - 1.96 * sqrt(`varICrr'))
	local UCIrr =  exp(`logRRtmle' + 1.96 * sqrt(`varICrr'))

// Statistical inference OR (We do not use the log of the OR)
	local ORtmle = (`Q1'  / (1 - `Q1')) / (`Q0' / ((1 - `Q0')))
	gen `ICor' = (1/ (`Q1' *(1 - `Q1')) * d1) - (1/ (`Q0' *(1 - `Q0')) * d0) // Partial derivatives (use wolfram to check: d/dx log(x/(1-x))+ log(y/(1-y)))
	qui sum `ICor'
	local varICor = r(Var)/r(N)
	local LCIOr = exp(log(`ORtmle') - 1.96 * sqrt(`varICor'))
	local UCIOr = exp(log(`ORtmle') + 1.96 * sqrt(`varICor'))
	
	/*
	// Statistical inference OR
		gen `ICor' = ((1 - `Q0') / `Q0' / (1 - `Q1')^2) * d1 - (`Q1' / (1 - `Q1') / `Q0'^2) * d0
		qui sum `ICor'
		local varICor = r(Var)/r(N)

		local LCIOr =  `ORtmle' - 1.96 * sqrt(`varICor')
		local UCIOr =  `ORtmle' + 1.96 * sqrt(`varICor')
	*/
	
	
// Display Results of ATE
	return scalar CRR = `RRtmle'
	return scalar SE_log_CRR  = sqrt(`varICrr')
	return scalar MOR = `ORtmle'
	return scalar SE_log_MOR  = sqrt(`varICor')

	if $flag==1 {
	disp as text "{hline 63}"
	di "         {c |}" "    ATE         SE     P-value           95% CI"
	disp as text "{hline 63}"
	disp as text "TMLE:    {c |}" %7.4f as result return(ATEtmle) "    " %7.4f as result return(ATE_SE_tmle) "     " %7.4f as result return(ATE_pvalue) as text "     (" %7.4f as result return(ATE_LCIa) as text ","  %7.4f as result return(ATE_UCIa) as text " )"
	disp as text "{hline 63}"
	disp as text " "
	}
	else if $flag!=1{
	disp as text "{hline 63}"
	di "         {c |}" "    ATE         SE     P-value           95% CI"
	disp as text "{hline 63}"
	disp as text "TMLE:    {c |}" %7.1f as result return(ATEtmle) "    " %7.1f as result return(ATE_SE_tmle) "     " %7.4f as result return(ATE_pvalue) as text "     (" %7.1f as result return(ATE_LCIa) as text ","  %7.1f as result return(ATE_UCIa) as text " )"
	disp as text "{hline 63}"
	disp as text " "
	}

	local rrbin ""CRR: "%4.2f `RRtmle'  "; 95%CI:("%3.2f `LCIrr' ", "%3.2f `UCIrr' ")""
	local orbin ""MOR: "%4.2f `ORtmle'  "; 95%CI:("%3.2f `LCIOr' ", "%3.2f `UCIOr' ")""

	* Display the results of CRR and MOR
	disp as text "{hline 51}"
	di "                           Estimate          95% CI"
	disp as text "{hline 51}"
	disp as text "Causal Risk Ratio:      " "{c |}      "  %04.2f as result `RRtmle' as text "     (" %03.2f as result `LCIrr' as text ","  %03.2f as result `UCIrr' as text ")"
	disp as text "Marginal Odds Ratio:    " "{c |}      "  %04.2f as result `ORtmle' as text "     (" %03.2f as result `LCIOr' as text "," %03.2f as result `UCIOr' as text ")"
	disp as text "{hline 51}"

// Drop the variables if the elements option is not specified
		if $keepvars == 0 {
			drop d1 d0
			drop QAW Q1W Q0W
			drop Q1star Qa1star Q0star Qa0star
			drop ATE IC
			drop Y A
			drop POM1 POM0 ps
			capture drop cin
			capture drop ytempvar
		}

// Rename and label the variables if the elements option *is* specified
		if $keepvars == 1 {
			* Label the variables
				lab var d1 "Parameter for the influence curve"
				lab var d0 "Parameter for the influence curve"
				lab var QAW "Initial prediction of the outcome"
				lab var Q1W "Initial prediction of the outcome for A = 1"
				lab var Q0W "Initial prediction of the outcome for A = 0"
				lab var Q1star "Update of the initial prediction for A = 1"
				lab var Qa1star "Update of the initial prediction for A = 1"
				lab var Q0star "Update of the initial prediction for A = 0"
				lab var Qa0star "Update of the initial prediction for A = 0"
				lab var A "Exposure/Treatment"
				lab var Y "Outcome"
				lab var ATE "Average Treatment Effect"
				lab var IC "Influence Curve"
				lab var Q1star "Update of initial plug-in estimate for A=1"
				lab var Q0star "Update of initial plug-in estimate for A=0"
				lab var POM1 "Potential Outcome Y(1)"
				lab var POM0 "Potential Otucome Y(0)"
				lab var ps "Propensity Score"
				//lab var cin "Range of Y"
				capture drop ytempvar
			* Rename the elements variables
				foreach var of varlist d1 d0 QAW Q1W Q0W Q1star Qa1star Q0star Qa0star ATE IC  Y A  POM1 POM0 ps {
					rename `var' _`var'
				}
		}


	// Clean up
	quietly: rm SLS.R
	// quietly: rm SLS.Rout
	quietly: rm data2.dta
	quietly: rm data.csv
	quietly: rm fulldata.csv
	// quietly: rm .RData
	quietly: memory clean
end









program cvtmle, rclass
	*di "Starting cvtmle..."

preserve
	forval fold = 1/$folds {
		
		if `fold' == 1 _dots 0, title (Performing $folds folds) 
		
		* Load the data to be used for the folds
		qui use "cvdata.dta", clear
		
		* Create a variable indicating the fold number
		capture drop foldnumber
		qui gen foldnumber = `fold'
		
		* Count number of observations
		qui count
		local nobstotal = r(N)
		
		* Export to csv file
		qui: export delimited $varusedforcv foldnumber rowid using "cvdata.csv", nolabel replace
		
		* Rcode for the SuperLearner
		qui: set more off
		qui: file close _all
		qui: file open rcode using SLS.R, write replace
		qui: file write rcode ///
			`"set.seed(123)"' _newline ///
			`"list.of.packages <- c("foreign","SuperLearner","dplyr")"' _newline ///
			`"new.packages <- list.of.packages[!(list.of.packages %in% installed.packages()[,"Package"])]"' _newline ///
			`"if(length(new.packages)) install.packages(new.packages, repos='http://cran.us.r-project.org')"' _newline ///
			`"library(SuperLearner)"' _newline ///
			`"library(foreign)"' _newline ///
			`"library(dplyr)"' _newline ///
			`"data <- read.csv("cvdata.csv", sep=",")"' _newline ///
			`"attach(data)"' _newline ///
			`"SL.library <- c("SL.glm","SL.step","SL.glm.interaction")"' _newline ///
			`"nvar <- dim(data)[[2]]"' _newline ///
			`"f <- foldnumber[1]"' _newline ///
			`"data <- data %>% select(-c("foldnumber"))"' _newline ///
			`"tdata <- data[foldid!=f,]"' _newline ///		
			`"vdataAa <- data[foldid==f,]"' _newline ///
			`"vdataA1 <- vdataAa; vdataA1[,2] <- 1"' _newline ///
			`"vdataA0 <- vdataAa; vdataA0[,2] <- 0"' _newline ///
			`"Vn <- nrow(vdataAa)"' _newline ///
			`"Tn <- nrow(tdata)"' _newline ///
			`"Y <- tdata[,1]"' _newline ///
			`"A <- tdata[,2]"' _newline ///
			`"X <- tdata %>% select(-c(1,    "foldid", "rowid"))"' _newline ///
			`"W <- tdata %>% select(-c(1, 2, "foldid", "rowid"))"' _newline ///
			`"tset <- try(SuperLearner(Y = tdata[,1], X = X, SL.library=SL.library, family = "binomial", method ="method.NNLS"), silent=TRUE)"' _newline ///
			`"QAW <- predict(tset, newdata=vdataAa)[[1]]"' _newline ///
			`"Q1W <- predict(tset, newdata=vdataA1)[[1]]"' _newline ///
			`"Q0W <- predict(tset, newdata=vdataA0)[[1]]"' _newline ///
			`"QAWt <- predict(tset)[[1]]"' _newline ///
			`"data <- cbind(vdataAa,QAW,Q1W,Q0W)"' _newline ///
			`"write.dta(data, "data2.dta")"' 
		qui: file close rcode
		
		* Specify the operating system
	if "`c(os)'" == "MacOSX" {
	qui shell "/usr/local/bin/r" CMD BATCH SLS.R
	}
	else{
	// Write batch file to find R.exe path and R version
		set more off
		qui: file close _all
		qui: file open bat using setup.bat, write replace
		qui: file write bat ///
		`"@echo off"' _newline ///
		`"SET PATHROOT=C:\Program Files\R\"' _newline ///
		`"echo Locating path of R..."' _newline ///
		`"echo."' _newline ///
		`"if not exist "%PATHROOT%" goto:NO_R"' _newline ///
		`":NO_R"' _newline ///
		`"echo R is not installed in your system."' _newline ///
		`"echo."' _newline ///
		`"echo Download it from https://cran.r-project.org/bin/windows/base/"' _newline ///
		`"echo Install it and re-run this script"' _newline ///
		`":DONE"' _newline ///
		`"echo."' _newline ///
		`"pause"'
		qui: file close bat
		
		// 		`"for /f "delims=" %%r in (' dir /b "%PATHROOT%R*" ') do ("' _newline ///
		// 				`"echo Found %%r"' _newline ///
		// 				`"echo ! "%PATHROOT%%%r\bin\x64\R.exe" CMD BATCH SLS.R > runr.do"' _newline ///
		// 				`"echo All set!"' _newline ///
		// 				`"goto:DONE"' _newline ///
		// 		`")"' _newline ///
		
		// 		//Run batch
		// 		! setup.bat
		// 		//Run R
		// 		do runr.do
		
		* Check you have 'rscript' command installed.
		qui: net install rscript, from("https://raw.githubusercontent.com/reifjulian/rscript/master") replace
		
		* Run SLS.R in Windows - 'rscript' command will search for common folders to find R.exe
		qui: rscript using SLS.R
	}
		
	
	* Load the data created by R
		quietly: use "data2.dta", clear
		local folddataname "folddata`fold'"
		qui save "`folddataname'.dta", replace
		
	*if $progress == 1 di "Done fold `fold'"
	
	_dots `fold' 0
		
	}
	
	
	* Combine the initial estimate into one data set
		* Append the data sets
		qui: use "folddata1.dta", clear
		forval fold = 2/$folds {
			qui: append using "folddata`fold'"
		}

	
	* SuperLearner for the propensity score
		* Export to csv file
		qui export delimited using "cvdata.csv", nolabel replace
		
		quietly: rm data2.dta
		
	// Write R Code dependencies: foreign Superlearner
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
			`"data <- read.csv("cvdata.csv", sep=",")"' _newline ///
			`"attach(data)"' _newline ///
			`"SL.library <- c("SL.glm","SL.step","SL.glm.interaction")"' _newline ///
			`"n <- nrow(data)"' _newline ///
			`"dt <- data"' _newline ///
			`"dt\$foldid <- NULL"' _newline ///
			`"dt\$QAW <- NULL"' _newline ///
			`"dt\$Q1W <- NULL"' _newline ///
			`"dt\$Q0W <- NULL"' _newline ///
			`"dt\$rowid <- NULL"' _newline ///
			`"nvar <- dim(dt)[[2]]"' _newline ///
			`"Y <- data[,1]"' _newline ///
			`"A <- data[,2]"' _newline ///
			`"X <- data[,2:nvar]"' _newline ///
			`"W <- as.data.frame(data[,3:nvar])"' _newline ///
			`"g <- suppressWarnings(SuperLearner(Y = data[,2], X = W, SL.library = SL.library, family = "binomial", method = "method.NNLS"))"' _newline ///
			`"ps <- g[[4]]"' _newline ///
			`"ps[ps<0.025] <- 0.025"' _newline ///
			`"ps[ps>0.975] <- 0.975"' _newline ///
			`"data <- cbind(data,ps,Y,A)"' _newline ///
			`"write.dta(data, "data2.dta")"'
	qui: file close rcode


	
	* Specify the operating system
	if "`c(os)'" == "MacOSX" {
	qui shell "/usr/local/bin/r" CMD BATCH SLS.R
	}
	else{
	// Write batch file to find R.exe path and R version
		set more off
		qui: file close _all
		qui: file open bat using setup.bat, write replace
		qui: file write bat ///
		`"@echo off"' _newline ///
		`"SET PATHROOT=C:\Program Files\R\"' _newline ///
		`"echo Locating path of R..."' _newline ///
		`"echo."' _newline ///
		`"if not exist "%PATHROOT%" goto:NO_R"' _newline ///
		`":NO_R"' _newline ///
		`"echo R is not installed in your system."' _newline ///
		`"echo."' _newline ///
		`"echo Download it from https://cran.r-project.org/bin/windows/base/"' _newline ///
		`"echo Install it and re-run this script"' _newline ///
		`":DONE"' _newline ///
		`"echo."' _newline ///
		`"pause"'
		qui: file close bat
		
		// 		`"for /f "delims=" %%r in (' dir /b "%PATHROOT%R*" ') do ("' _newline ///
		// 				`"echo Found %%r"' _newline ///
		// 				`"echo ! "%PATHROOT%%%r\bin\x64\R.exe" CMD BATCH SLS.R > runr.do"' _newline ///
		// 				`"echo All set!"' _newline ///
		// 				`"goto:DONE"' _newline ///
		// 		`")"' _newline ///
		
		// 		//Run batch
		// 		! setup.bat
		// 		//Run R
		// 		do runr.do
		
		* Check you have 'rscript' command installed.
		qui: net install rscript, from("https://raw.githubusercontent.com/reifjulian/rscript/master") replace
		
		* Run SLS.R in Windows - 'rscript' command will search for common folders to find R.exe
		qui: rscript using SLS.R
	}
		
		
// Read Revised Data Back to Stata
	qui clear
	quietly: use "data2.dta", clear
	*import delimited "data2.csv", case(preserve) clear 
	qui cap drop X__000000
	tempvar logQAW logQ1W logQ0W HAW H1W H0W eps1 eps2 eps ATE ICrr ICor

	
// Q to logit scale
	gen `logQAW' = log(QAW / (1 - QAW))
	gen `logQ1W' = log(Q1W / (1 - Q1W))
	gen `logQ0W' = log(Q0W / (1 - Q0W))


// Clever covariate HAW
	gen  `HAW' = (A / ps) - ((1 - A) / (1 - ps))
	gen  `H1W' = A / ps
	gen  `H0W' = (1 - A) / (1 - ps)

	
// Estimation of the substitution parameter (Epsilon)
	qui glm Y `H1W' `H0W', fam(binomial) offset(`logQAW') robust noconstant
	mat a= e(b)
	gen `eps1' = a[1,1]
	gen `eps2' = a[1,2]

	qui glm Y `HAW', fam(binomial) offset(`logQAW') robust noconstant
	mat a= e(b)
	gen `eps' = a[1,1]
	
	
// Targeted ATE, update from Q̅^0 (A,W) to Q̅^1 (A,W)
	gen double Qa0star = exp(`H0W'*`eps' + `logQ0W')/(1 + exp(`H0W'*`eps' + `logQ0W'))
	gen double Qa1star = exp(`H1W'*`eps' + `logQ1W')/(1 + exp(`H1W'*`eps' + `logQ1W'))

	gen double Q0star = exp(`H0W'*`eps2' + `logQ0W')/(1 + exp(`H0W'*`eps2' + `logQ0W'))
	gen double Q1star = exp(`H1W'*`eps1' + `logQ1W')/(1 + exp(`H1W'*`eps1' + `logQ1W'))

	// gen double cin = ($b - $a) 

	gen double POM1 = cond($flag == 1, Qa1star, (Qa1star * ($b - $a )) + $a , .)
	gen double POM0 = cond($flag == 1, Qa0star, (Qa0star * ($b - $a )) + $a , .)

	di as text " "
	sum POM1 POM0 ps
	di as text " "

// Estimating the updated targeted ATE binary outcome
	gen double ATE = cond($flag == 1, (Qa1star - Qa0star), ((Qa1star * ($b - $a )) + $a ) - ((Qa0star * ($b - $a )) + $a ) , .)
	qui sum ATE
	return scalar ATEtmle = r(mean)

// Relative risk
	qui sum Q1star
	local Q1 = r(mean)
	qui sum Q0star
	local Q0 = r(mean)

// Relative risk and Odds ratio
	local RRtmle = `Q1'/`Q0'
	local logRRtmle = log(`Q1') - log(`Q0')
	local ORtmle = (`Q1' * (1 - `Q0')) / ((1 - `Q1') * `Q0')


// Statistical inference (Efficient Influence Curve)
	gen d1 = cond($flag == 1,(A * (Y - Q1star) / ps) + Q1star - `Q1',(A * (Y - Qa1star) / ps) + Qa1star - `Q1' ,.)
	gen d0 = cond($flag == 1,(1 - A) * (Y - Q0star) / (1 - ps) + Q0star - `Q0',(1 - A) * (Y - Qa0star) / (1 - ps) + Qa0star - `Q0' ,.)
	gen IC = cond($flag == 1,(d1 - d0), ((d1 * ($b - $a )) + $a ) - ((d0 * ($b - $a )) + $a ) , .)
	qui sum IC
	return scalar ATE_SE_tmle = sqrt(r(Var)/r(N))

// Statistical inference ATE
	return scalar ATE_pvalue =  2 * (normalden(abs(return(ATEtmle) / (return(ATE_SE_tmle)))))
	return scalar ATE_LCIa   =  return(ATEtmle) - 1.96 * return(ATE_SE_tmle)
	return scalar ATE_UCIa   =  return(ATEtmle) + 1.96 * return(ATE_SE_tmle)

// Statistical inference RR
	gen `ICrr' = (1/`Q1' * d1) + ((1/`Q0') * d0)
	qui sum `ICrr'
	local varICrr = r(Var)/r(N)
	
	local LCIrr =  exp(`logRRtmle' - 1.96 * sqrt(`varICrr'))
	local UCIrr =  exp(`logRRtmle' + 1.96 * sqrt(`varICrr'))
	
// Statistical inference OR (We do not use the log of the OR)
	local ORtmle = (`Q1' * (1 - `Q0')) / ((1 - `Q1') * `Q0')
	gen `ICor' = (1/ (`Q1' *(1 - `Q1')) * d1) - (1/ (`Q0' *(1 - `Q0')) * d0) // Partial derivatives (use wolfram to check: d/dx log(x/(1-x))+ log(y/(1-y)))
	qui sum `ICor'
	local varICor = r(Var)/r(N)
	local LCIOr = exp(log(`ORtmle') - 1.96 * sqrt(`varICor'))
	local UCIOr = exp(log(`ORtmle') + 1.96 * sqrt(`varICor'))
	*di sqrt(`varICor')
	
// Display Results of ATE
	return scalar CRR = `RRtmle'
	return scalar SE_log_CRR  = sqrt(`varICrr')
	return scalar MOR = `ORtmle'
	return scalar SE_log_MOR  = sqrt(`varICor')

	if $flag==1 {
	disp as text "{hline 63}"
	di "         {c |}" "    ATE         SE     P-value           95% CI"
	disp as text "{hline 63}"
	disp as text "TMLE:    {c |}" %7.4f as result return(ATEtmle) "    " %7.4f as result return(ATE_SE_tmle) "     " %7.4f as result return(ATE_pvalue) as text "     (" %7.4f as result return(ATE_LCIa) as text ","  %7.4f as result return(ATE_UCIa) as text " )"
	disp as text "{hline 63}"
	disp as text " "
	}
	else if $flag!=1{
	disp as text "{hline 63}"
	di "         {c |}" "    ATE         SE     P-value           95% CI"
	disp as text "{hline 63}"
	disp as text "TMLE:    {c |}" %7.1f as result return(ATEtmle) "    " %7.1f as result return(ATE_SE_tmle) "     " %7.4f as result return(ATE_pvalue) as text "     (" %7.1f as result return(ATE_LCIa) as text ","  %7.1f as result return(ATE_UCIa) as text " )"
	disp as text "{hline 63}"
	disp as text " "
	}

	local rrbin ""CRR: "%4.2f `RRtmle'  "; 95%CI:("%3.2f `LCIrr' ", "%3.2f `UCIrr' ")""
	local orbin ""MOR: "%4.2f `ORtmle'  "; 95%CI:("%3.2f `LCIOr' ", "%3.2f `UCIOr' ")""

	* Display the results of CRR and MOR
	disp as text "{hline 51}"
	di "                           Estimate          95% CI"
	disp as text "{hline 51}"
	disp as text "Causal Risk Ratio:      " "{c |}      "  %04.2f as result `RRtmle' as text "     (" %03.2f as result `LCIrr' as text ","  %03.2f as result `UCIrr' as text ")"
	disp as text "Marginal Odds Ratio:    " "{c |}      "  %04.2f as result `ORtmle' as text "     (" %03.2f as result `LCIOr' as text "," %03.2f as result `UCIOr' as text ")"
	disp as text "{hline 51}"

	
// Display covariate balance table

	if $covbalancetable == 0 {
		
			* Create macros from the varlist
				tokenize $variablelist
				local outcome `1'
				macro shift
				local exposure `1'
				macro shift
				local varlist `*'

			* Create the IPW weights
				*di "Capture dropping _ipw 1"
				qui capture drop _ipw
				qui gen _ipw = .
				qui replace _ipw = (`exposure'==1) / ps if `exposure'==1
				qui replace _ipw = (`exposure'==0) / (1- ps) if `exposure'==0
	}


	if $covbalancetable == 1 {

		* Create macros from the varlist
			tokenize $variablelist
			local outcome `1'
			macro shift
			local exposure `1'
			macro shift
			local varlist `*'

		* Create the IPW weights
			*di "Capture dropping _ipw 1"
			qui capture drop _ipw
			qui gen _ipw = .
			qui replace _ipw = (`exposure'==1) / ps if `exposure'==1
			qui replace _ipw = (`exposure'==0) / (1- ps) if `exposure'==0

		* Layout of the results
			di as text "{hline 67}"
			di as text "                 Standardised Differences            Variance ratio"
			di as text "                          Raw    Weighted           Raw    Weighted"
			di as text "{hline 67}"

			
		* Calculate the covariate balance
			foreach var in `varlist' {
					di as text "`var'"

					* Raw SMD
					qui summarize `var' if `exposure'==1
					local m1 = r(mean)
					local v1 = r(Var)
					qui summarize `var' if `exposure'==0
					local m0 = r(mean)
					local v0 = r(Var)
					* Calculate the Standardised "mean difference"
					local rSMD = (`m1' - `m0') / sqrt( (`v1' + `v0') /2 )

					* Weighted SMD
					qui summarize `var' [iw=_ipw] if `exposure'==1
					local m1 = r(mean)
					local v1 = r(Var)
					qui summarize `var' [iw=_ipw] if `exposure'==0
					local m0 = r(mean)
					local v0 = r(Var)
					* Calculate the Standardised "mean difference"
					local wSMD = (`m1' - `m0') / sqrt( (`v1' + `v0') /2 )

					* Raw VR
					qui sum `var' if `exposure' ==1
					local v1 = r(Var)
					qui sum `var' if `exposure' ==0
					local v0 = r(Var)
					* Calculate the variance ratio
					local rVR = `v1' / `v0'

					* Weighted VR
					qui sum `var' [iw=_ipw] if `exposure' ==1
					local v1 = r(Var)
					qui sum `var' [iw=_ipw] if `exposure' ==0
					local v0 = r(Var)
					* Calculate the variance ratio
					local wVR `v1' / `v0'

					* Display the values
					di as text "                    " %9.7g as result `rSMD' as text "   " %9.7g as result `wSMD' as text "     " %9.7g as result `rVR' as text "   " %9.7g as result `wVR'
			}
			di as text "{hline 67}"

	}
	
	
// Create density plot to check positivity assumption
	if $covbalancetable == 1 {
	*preserve
	*sort A
	*qui by A: summarize ps
	qui kdensity ps if A==1, generate(x1pointsa d1A) nograph n(10000)
	qui kdensity ps if A==0, generate(x0pointsa d0A) nograph n(10000)
	label variable d1A "A = 1"
	label variable d0A "A = 0"
	twoway (line d0A x0pointsa) || ///
			(line d1A x1pointsa), ///
			xtitle("Propensity score") ///
			ytitle("Density") ///
			graphregion(color(white)) bgcolor(white) plotregion(fcolor(white)) ///
			legend(order(2 "Treated" 1 "Not treated") position(11) cols(1) ring(0) ///
				region(style(none)))
	*restore
	}


// Drop the variables if the elements option is not specified
		if $keepvars == 0 {
			drop d1 d0
			drop QAW Q1W Q0W
			drop Q1star Qa1star Q0star Qa0star
			drop ATE IC
			rename Y _Y
			rename A _A
			drop _Y _A
			*di "Capture dropping _ipw 2"
			capture drop _ipw
			drop POM1 POM0 ps
			capture drop cin
			drop foldid rowid
			capture drop x1pointsa x0pointsa d1A d0A
			capture drop ytempvar
		}

// Rename and label the variables if the elements option *is* specified
		if $keepvars == 1 {
			
			capture drop x1pointsa x0pointsa d1A d0A
			
			* Label the variables
				lab var d1 "Parameter for the influence curve"
				lab var d0 "Parameter for the influence curve"
				lab var QAW "Initial prediction of the outcome"
				lab var Q1W "Initial prediction of the outcome for A = 1"
				lab var Q0W "Initial prediction of the outcome for A = 0"
				lab var Q1star "Update of the initial prediction for A = 1"
				lab var Qa1star "Update of the initial prediction for A = 1"
				lab var Q0star "Update of the initial prediction for A = 0"
				lab var Qa0star "Update of the initial prediction for A = 0"
				lab var A "Exposure/Treatment"
				lab var Y "Outcome"
				lab var ATE "Average Treatment Effect"
				lab var IC "Influence Curve"
				lab var Q1star "Update of initial plug-in estimate for A=1"
				lab var Q0star "Update of initial plug-in estimate for A=0"
				lab var POM1 "Potential Outcome Y(1)"
				lab var POM0 "Potential Otucome Y(0)"
				lab var ps "Propensity Score"
				//lab var cin "Range of Y"
				lab var _ipw "Inverse probability of treatment weights"
				capture drop ytempvar
			* Rename the elements variables
				foreach var of varlist d1 d0 QAW Q1W Q0W Q1star Qa1star Q0star Qa0star ATE IC Y A  POM1 POM0 ps {
					rename `var' _`var'
				}
				*di "Keeping _ipw 1"
				keep rowid _d1 _d0 _QAW _Q1W _Q0W _Q1star _Qa1star _Q0star _Qa0star _ATE _IC _Y _A _POM1 _POM0 _ps _ipw 
				qui save "elementsdata.dta", replace
		}


restore

	// Drop the temporary variable for Y
	capture drop ytempvar

	// Clean up
	quietly: rm SLS.R
	// quietly: rm SLS.Rout
	quietly: rm data2.dta
	quietly: rm cvdata.csv
	quietly: rm cvdata.dta
	quietly: rm fulldata.csv
	// quietly: rm .RData
	quietly: memory clean

* Note the completion of the program
	//di "Success"



	// Merge in the elements
		if $keepvars == 0 {
			qui drop foldid rowid
		}

	// Merge in the elements
		if $keepvars == 1 {
			qui merge 1:1 rowid using "elementsdata.dta"
		}
	
end





program cvtmleglsrfbin, rclass

di "Starting cvtmleglsrf... please be patient, additional algorithms are computationally intensive."

preserve

	forval fold = 1/$folds {
		
		if `fold' == 1 _dots 0, title (Performing $folds folds) 
		
		* Load the data to be used for the folds
		qui use "cvdata.dta", clear
		
		* Create a variable indicating the fold number
		capture drop foldnumber
		qui gen foldnumber = `fold'
		
		* Count number of observations
		qui count
		local nobstotal = r(N)
		
		* Export to csv file
		qui: export delimited $varusedforcv foldnumber rowid using "cvdata.csv", nolabel replace
		
		* Rcode for the SuperLearner
		qui: set more off
		qui: file close _all
		qui: file open rcode using SLS.R, write replace
		qui: file write rcode ///
			`"set.seed(123)"' _newline ///
			`"list.of.packages <- c("foreign","SuperLearner","dplyr","glmnet","randomForest")"' _newline ///
			`"new.packages <- list.of.packages[!(list.of.packages %in% installed.packages()[,"Package"])]"' _newline ///
			`"if(length(new.packages)) install.packages(new.packages, repos='http://cran.us.r-project.org')"' _newline ///
			`"library(SuperLearner)"' _newline ///
			`"library(foreign)"' _newline ///
			`"library(dplyr)"' _newline ///
			`"library(randomForest)"' _newline ///
			`"library(glmnet)"' _newline ///
			`"data <- read.csv("cvdata.csv", sep=",")"' _newline ///
			`"attach(data)"' _newline ///
			`"SL.library <- c("SL.glm","SL.step","SL.glm.interaction","SL.gam","SL.glmnet","SL.randomForest")"' _newline ///
			`"nvar <- dim(data)[[2]]"' _newline ///
			`"f <- foldnumber[1]"' _newline ///
			`"data <- data %>% select(-c("foldnumber"))"' _newline ///
			`"tdata <- data[foldid!=f,]"' _newline ///		
			`"vdataAa <- data[foldid==f,]"' _newline ///
			`"vdataA1 <- vdataAa; vdataA1[,2] <- 1"' _newline ///
			`"vdataA0 <- vdataAa; vdataA0[,2] <- 0"' _newline ///
			`"Vn <- nrow(vdataAa)"' _newline ///
			`"Tn <- nrow(tdata)"' _newline ///
			`"Y <- tdata[,1]"' _newline ///
			`"A <- tdata[,2]"' _newline ///
			`"X <- tdata %>% select(-c(1,    "foldid", "rowid"))"' _newline ///
			`"W <- tdata %>% select(-c(1, 2, "foldid", "rowid"))"' _newline ///
			`"tset <- try(SuperLearner(Y = tdata[,1], X = X, SL.library=SL.library, family = "binomial", method ="method.NNLS"), silent=TRUE)"' _newline ///
			`"QAW <- predict(tset, newdata=vdataAa)[[1]]"' _newline ///
			`"Q1W <- predict(tset, newdata=vdataA1)[[1]]"' _newline ///
			`"Q0W <- predict(tset, newdata=vdataA0)[[1]]"' _newline ///
			`"QAWt <- predict(tset)[[1]]"' _newline ///
			`"data <- cbind(vdataAa,QAW,Q1W,Q0W)"' _newline ///
			`"write.dta(data, "data2.dta")"' 
		qui: file close rcode

		
		
		* Specify the operating system
	if "`c(os)'" == "MacOSX" {
		qui shell "/usr/local/bin/r" CMD BATCH SLS.R
	}
	else{
	// Write batch file to find R.exe path and R version
		set more off
		qui: file close _all
		qui: file open bat using setup.bat, write replace
		qui: file write bat ///
		`"@echo off"' _newline ///
		`"SET PATHROOT=C:\Program Files\R\"' _newline ///
		`"echo Locating path of R..."' _newline ///
		`"echo."' _newline ///
		`"if not exist "%PATHROOT%" goto:NO_R"' _newline ///
		`":NO_R"' _newline ///
		`"echo R is not installed in your system."' _newline ///
		`"echo."' _newline ///
		`"echo Download it from https://cran.r-project.org/bin/windows/base/"' _newline ///
		`"echo Install it and re-run this script"' _newline ///
		`":DONE"' _newline ///
		`"echo."' _newline ///
		`"pause"'
		qui: file close bat
		
		// 		`"for /f "delims=" %%r in (' dir /b "%PATHROOT%R*" ') do ("' _newline ///
		// 				`"echo Found %%r"' _newline ///
		// 				`"echo ! "%PATHROOT%%%r\bin\x64\R.exe" CMD BATCH SLS.R > runr.do"' _newline ///
		// 				`"echo All set!"' _newline ///
		// 				`"goto:DONE"' _newline ///
		// 		`")"' _newline ///
		
		// 		//Run batch
		// 		! setup.bat
		// 		//Run R
		// 		do runr.do
		
		* Check you have 'rscript' command installed.
		qui: net install rscript, from("https://raw.githubusercontent.com/reifjulian/rscript/master") replace
		
		* Run SLS.R in Windows - 'rscript' command will search for common folders to find R.exe
		qui: rscript using SLS.R
	}
		
	
	* Load the data created by R
		quietly: use "data2.dta", clear
		local folddataname "folddata`fold'"
		qui save "`folddataname'.dta", replace
		
	* if $progress == 1 di "Done fold `fold'"
	
	_dots `fold' 0
		
	}
	
	
	
* Combine the initial estimate into one data set
	* Append the data sets
	qui use "folddata1.dta", clear
	forval fold = 2/$folds {
		qui append using "folddata`fold'"
	}

	
* SuperLearner for the propensity score
	* Export to csv file
	qui export delimited using "cvdata.csv", nolabel replace
	
	quietly: rm data2.dta
		
	// Write R Code dependencies: foreign Superlearner
	set more off
	qui: file close _all
	qui: file open rcode using SLS.R, write replace
	qui: file write rcode ///
			`"set.seed(123)"' _newline ///
			`"list.of.packages <- c("foreign","SuperLearner","dplyr","glmnet","randomForest")"' _newline ///
			`"new.packages <- list.of.packages[!(list.of.packages %in% installed.packages()[,"Package"])]"' _newline ///
			`"if(length(new.packages)) install.packages(new.packages, repos='http://cran.us.r-project.org')"' _newline ///
			`"library(SuperLearner)"' _newline ///
			`"library(foreign)"' _newline ///
			`"library(dplyr)"' _newline ///
			`"library(glmnet)"' _newline ///
			`"library(randomForest)"' _newline ///
			`"data <- read.csv("cvdata.csv", sep=",")"' _newline ///
			`"attach(data)"' _newline ///
			`"SL.library <- c("SL.glm","SL.step","SL.glm.interaction","SL.gam","SL.glmnet","SL.randomForest")"' _newline ///
			`"n <- nrow(data)"' _newline ///
			`"dt <- data"' _newline ///
			`"dt <- dt %>% select(-c("foldid"))"' _newline ///
			`"dt <- dt %>% select(-c("rowid"))"' _newline ///
			`"dt <- dt %>% select(-c("QAW"))"' _newline ///
			`"dt <- dt %>% select(-c("Q1W"))"' _newline ///
			`"dt <- dt %>% select(-c("Q0W"))"' _newline ///
			`"nvar <- dim(dt)[[2]]"' _newline ///
			`"Y <- data[,1]"' _newline ///
			`"A <- data[,2]"' _newline ///
			`"X <- data[,2:nvar]"' _newline ///
			`"W <- as.data.frame(data[,3:nvar])"' _newline ///
			`"g <- suppressWarnings(SuperLearner(Y = data[,2], X = W, SL.library = SL.library, family = "binomial", method = "method.NNLS"))"' _newline ///
			`"ps <- g[[4]]"' _newline ///
			`"ps[ps<0.025] <- 0.025"' _newline ///
			`"ps[ps>0.975] <- 0.975"' _newline ///
			`"data <- cbind(data,ps,Y,A)"' _newline ///
			`"write.dta(data, "data2.dta")"'
	qui: file close rcode

	* Specify the operating system
	if "`c(os)'" == "MacOSX" {
		qui shell "/usr/local/bin/r" CMD BATCH SLS.R
	}
	else{
	// Write batch file to find R.exe path and R version
		set more off
		qui: file close _all
		qui: file open bat using setup.bat, write replace
		qui: file write bat ///
		`"@echo off"' _newline ///
		`"SET PATHROOT=C:\Program Files\R\"' _newline ///
		`"echo Locating path of R..."' _newline ///
		`"echo."' _newline ///
		`"if not exist "%PATHROOT%" goto:NO_R"' _newline ///
		`":NO_R"' _newline ///
		`"echo R is not installed in your system."' _newline ///
		`"echo."' _newline ///
		`"echo Download it from https://cran.r-project.org/bin/windows/base/"' _newline ///
		`"echo Install it and re-run this script"' _newline ///
		`":DONE"' _newline ///
		`"echo."' _newline ///
		`"pause"'
		qui: file close bat
		
		// 		`"for /f "delims=" %%r in (' dir /b "%PATHROOT%R*" ') do ("' _newline ///
		// 				`"echo Found %%r"' _newline ///
		// 				`"echo ! "%PATHROOT%%%r\bin\x64\R.exe" CMD BATCH SLS.R > runr.do"' _newline ///
		// 				`"echo All set!"' _newline ///
		// 				`"goto:DONE"' _newline ///
		// 		`")"' _newline ///
		
		// 		//Run batch
		// 		! setup.bat
		// 		//Run R
		// 		do runr.do
		
		* Check you have 'rscript' command installed.
		qui: net install rscript, from("https://raw.githubusercontent.com/reifjulian/rscript/master") replace
		
		* Run SLS.R in Windows - 'rscript' command will search for common folders to find R.exe
		qui: rscript using SLS.R
	}
		
		
// Read Revised Data Back to Stata
	qui clear
	quietly: use "data2.dta", clear
	*import delimited "data2.csv", case(preserve) clear 
	qui cap drop X__000000
	tempvar logQAW logQ1W logQ0W HAW H1W H0W eps1 eps2 eps ATE ICrr ICor

// Q to logit scale
	gen `logQAW' = log(QAW / (1 - QAW))
	gen `logQ1W' = log(Q1W / (1 - Q1W))
	gen `logQ0W' = log(Q0W / (1 - Q0W))

// Clever covariate HAW
	gen  `HAW' = (A / ps) - ((1 - A) / (1 - ps))
	gen  `H1W' = A / ps
	gen  `H0W' = (1 - A) / (1 - ps)

// Estimation of the substitution parameter (Epsilon)
	qui glm Y `H1W' `H0W', fam(binomial) offset(`logQAW') robust noconstant
	mat a= e(b)
	gen `eps1' = a[1,1]
	gen `eps2' = a[1,2]

	qui glm Y `HAW', fam(binomial) offset(`logQAW') robust noconstant
	mat a= e(b)
	gen `eps' = a[1,1]
	
	
// Targeted ATE, update from Q̅^0 (A,W) to Q̅^1 (A,W)
	gen double Qa0star = exp(`H0W'*`eps' + `logQ0W')/(1 + exp(`H0W'*`eps' + `logQ0W'))
	gen double Qa1star = exp(`H1W'*`eps' + `logQ1W')/(1 + exp(`H1W'*`eps' + `logQ1W'))

	gen double Q0star = exp(`H0W'*`eps2' + `logQ0W')/(1 + exp(`H0W'*`eps2' + `logQ0W'))
	gen double Q1star = exp(`H1W'*`eps1' + `logQ1W')/(1 + exp(`H1W'*`eps1' + `logQ1W'))

	// gen double cin = ($b - $a)

	gen double POM1 = cond($flag == 1, Qa1star, (Qa1star * ($b - $a )) + $a , .)
	gen double POM0 = cond($flag == 1, Qa0star, (Qa0star * ($b - $a )) + $a , .)
	
	di as text " "
	sum POM1 POM0 ps
	di as text " "

// Estimating the updated targeted ATE binary outcome
	gen double ATE = cond($flag == 1, (Qa1star - Qa0star), ((Qa1star * ($b - $a )) + $a ) - ((Qa0star * ($b - $a )) + $a ) , .)
	qui sum ATE
	return scalar ATEtmle = r(mean)

// Relative risk
	qui sum Q1star
	local Q1 = r(mean)
	qui sum Q0star
	local Q0 = r(mean)

// Relative risk and Odds ratio
	local RRtmle = `Q1'/`Q0'
	local logRRtmle = log(`Q1') - log(`Q0')
	local ORtmle = (`Q1' * (1 - `Q0')) / ((1 - `Q1') * `Q0')


// Statistical inference (Efficient Influence Curve)
	gen d1 = cond($flag == 1,(A * (Y - Q1star) / ps) + Q1star - `Q1',(A * (Y - Qa1star) / ps) + Qa1star - `Q1' ,.)
	gen d0 = cond($flag == 1,(1 - A) * (Y - Q0star) / (1 - ps) + Q0star - `Q0',(1 - A) * (Y - Qa0star) / (1 - ps) + Qa0star - `Q0' ,.)
	gen IC = cond($flag == 1,(d1 - d0), ((d1 * ($b - $a )) + $a ) - ((d0 * ($b - $a )) + $a ) , .)
	qui sum IC
	return scalar ATE_SE_tmle = sqrt(r(Var)/r(N))

// Statistical inference ATE
	return scalar ATE_pvalue =  2 * (normalden(abs(return(ATEtmle) / (return(ATE_SE_tmle)))))
	return scalar ATE_LCIa   =  return(ATEtmle) - 1.96 * return(ATE_SE_tmle)
	return scalar ATE_UCIa   =  return(ATEtmle) + 1.96 * return(ATE_SE_tmle)

// Statistical inference RR
	gen `ICrr' = (1/`Q1' * d1) + ((1/`Q0') * d0)
	qui sum `ICrr'
	local varICrr = r(Var)/r(N)
	
	local LCIrr =  exp(`logRRtmle' - 1.96 * sqrt(`varICrr'))
	local UCIrr =  exp(`logRRtmle' + 1.96 * sqrt(`varICrr'))
	
// Statistical inference OR (We do not use the log of the OR)
	local ORtmle = (`Q1' * (1 - `Q0')) / ((1 - `Q1') * `Q0')
	gen `ICor' = (1/ (`Q1' *(1 - `Q1')) * d1) - (1/ (`Q0' *(1 - `Q0')) * d0) // Partial derivatives (use wolfram to check: d/dx log(x/(1-x))+ log(y/(1-y)))
	qui sum `ICor'
	local varICor = r(Var)/r(N)
	local LCIOr = exp(log(`ORtmle') - 1.96 * sqrt(`varICor'))
	local UCIOr = exp(log(`ORtmle') + 1.96 * sqrt(`varICor'))
	*di sqrt(`varICor')
	
// Display Results of ATE
	return scalar CRR = `RRtmle'
	return scalar SE_log_CRR  = sqrt(`varICrr')
	return scalar MOR = `ORtmle'
	return scalar SE_log_MOR  = sqrt(`varICor')

	if $flag==1 {
		disp as text "{hline 63}"
		di "         {c |}" "    ATE         SE     P-value           95% CI"
		disp as text "{hline 63}"
		disp as text "TMLE:    {c |}" %7.4f as result return(ATEtmle) "    " %7.4f as result return(ATE_SE_tmle) "     " %7.4f as result return(ATE_pvalue) as text "     (" %7.4f as result return(ATE_LCIa) as text ","  %7.4f as result return(ATE_UCIa) as text " )"
		disp as text "{hline 63}"
		disp as text " "
	}
	else if $flag!=1{
		disp as text "{hline 63}"
		di "         {c |}" "    ATE         SE     P-value           95% CI"
		disp as text "{hline 63}"
		disp as text "TMLE:    {c |}" %7.1f as result return(ATEtmle) "    " %7.1f as result return(ATE_SE_tmle) "     " %7.4f as result return(ATE_pvalue) as text "     (" %7.1f as result return(ATE_LCIa) as text ","  %7.1f as result return(ATE_UCIa) as text " )"
		disp as text "{hline 63}"
		disp as text " "
	}

	local rrbin ""CRR: "%4.2f `RRtmle'  "; 95%CI:("%3.2f `LCIrr' ", "%3.2f `UCIrr' ")""
	local orbin ""MOR: "%4.2f `ORtmle'  "; 95%CI:("%3.2f `LCIOr' ", "%3.2f `UCIOr' ")""

	* Display the results of CRR and MOR
	disp as text "{hline 51}"
	di "                           Estimate          95% CI"
	disp as text "{hline 51}"
	disp as text "Causal Risk Ratio:      " "{c |}      "  %04.2f as result `RRtmle' as text "     (" %03.2f as result `LCIrr' as text ","  %03.2f as result `UCIrr' as text ")"
	disp as text "Marginal Odds Ratio:    " "{c |}      "  %04.2f as result `ORtmle' as text "     (" %03.2f as result `LCIOr' as text "," %03.2f as result `UCIOr' as text ")"
	disp as text "{hline 51}"

	
// Display covariate balance table

	if $covbalancetable == 0 {
		
			* Create macros from the varlist
				tokenize $variablelist
				local outcome `1'
				macro shift
				local exposure `1'
				macro shift
				local varlist `*'

			* Create the IPW weights
				*di "Capture dropping _ipw 1"
				qui capture drop _ipw
				qui gen _ipw = .
				qui replace _ipw = (`exposure'==1) / ps if `exposure'==1
				qui replace _ipw = (`exposure'==0) / (1- ps) if `exposure'==0
	}


 if $covbalancetable == 1 {
 	
		* Create macros from the varlist
			tokenize $variablelist
			local outcome `1'
			macro shift
			local exposure `1'
			macro shift
			local varlist `*'

		* Create the IPW weights
			capture drop _ipw
			qui gen _ipw = .
			qui replace _ipw = (`exposure'==1) / ps if `exposure'==1
			qui replace _ipw = (`exposure'==0) / (1- ps) if `exposure'==0

		* Layout of the results
			di as text "{hline 67}"
			di as text "                 Standardised Differences            Variance ratio"
			di as text "                          Raw    Weighted           Raw    Weighted"
			di as text "{hline 67}"


		* Calculate the covariate balance
			foreach var in `varlist' {
					di as text "`var'"

					* Raw SMD
					qui summarize `var' if `exposure'==1
					local m1 = r(mean)
					local v1 = r(Var)
					qui summarize `var' if `exposure'==0
					local m0 = r(mean)
					local v0 = r(Var)
					* Calculate the Standardised "mean difference"
					local rSMD = (`m1' - `m0') / sqrt( (`v1' + `v0') /2 )

					* Weighted SMD
					qui summarize `var' [iw=_ipw] if `exposure'==1
					local m1 = r(mean)
					local v1 = r(Var)
					qui summarize `var' [iw=_ipw] if `exposure'==0
					local m0 = r(mean)
					local v0 = r(Var)
					* Calculate the Standardised "mean difference"
					local wSMD = (`m1' - `m0') / sqrt( (`v1' + `v0') /2 )

					* Raw VR
					qui sum `var' if `exposure' ==1
					local v1 = r(Var)
					qui sum `var' if `exposure' ==0
					local v0 = r(Var)
					* Calculate the variance ratio
					local rVR = `v1' / `v0'

					* Weighted VR
					qui sum `var' [iw=_ipw] if `exposure' ==1
					local v1 = r(Var)
					qui sum `var' [iw=_ipw] if `exposure' ==0
					local v0 = r(Var)
					* Calculate the variance ratio
					local wVR `v1' / `v0'

					* Display the values
					di as text "                    " %9.7g as result `rSMD' as text "   " %9.7g as result `wSMD' as text "     " %9.7g as result `rVR' as text "   " %9.7g as result `wVR'
			}
			di as text "{hline 67}"
	
	}
	

// Create density plot to check positivity assumption
	if $covbalancetable == 1 {
	*preserve
	*sort A
	*qui by A: summarize ps
	qui kdensity ps if A==1, generate(x1pointsa d1A) nograph n(10000)
	qui kdensity ps if A==0, generate(x0pointsa d0A) nograph n(10000)
	label variable d1A "A = 1"
	label variable d0A "A = 0"
	twoway (line d0A x0pointsa) || ///
			(line d1A x1pointsa), ///
			xtitle("Propensity score") ///
			ytitle("Density") ///
			graphregion(color(white)) bgcolor(white) plotregion(fcolor(white)) ///
			legend(order(2 "Treated" 1 "Not treated") position(11) cols(1) ring(0) ///
				region(style(none)))
	*restore
	}
	
	
// Drop the variables if the elements option is not specified
		if $keepvars == 0 {
			drop d1 d0
			drop QAW Q1W Q0W
			drop Q1star Qa1star Q0star Qa0star
			drop ATE IC
			rename Y _Y
			rename A _A
			drop _Y _A
			drop _ipw
			drop POM1 POM0 ps
			*drop cin
			drop foldid rowid
			capture drop x1pointsa x0pointsa d1A d0A
		}

// Rename and label the variables if the elements option *is* specified
		if $keepvars == 1 {
			
			capture drop x1pointsa x0pointsa d1A d0A
			
			* Label the variables
				lab var d1 "Parameter for the influence curve"
				lab var d0 "Parameter for the influence curve"
				lab var QAW "Initial prediction of the outcome"
				lab var Q1W "Initial prediction of the outcome for A = 1"
				lab var Q0W "Initial prediction of the outcome for A = 0"
				lab var Q1star "Update of the initial prediction for A = 1"
				lab var Qa1star "Update of the initial prediction for A = 1"
				lab var Q0star "Update of the initial prediction for A = 0"
				lab var Qa0star "Update of the initial prediction for A = 0"
				lab var A "Exposure/Treatment"
				lab var Y "Outcome"
				lab var ATE "Average Treatment Effect"
				lab var IC "Influence Curve"
				lab var Q1star "Update of initial plug-in estimate for A=1"
				lab var Q0star "Update of initial plug-in estimate for A=0"
				lab var POM1 "Potential Outcome Y(1)"
				lab var POM0 "Potential Otucome Y(0)"
				lab var ps "Propensity Score"
				*lab var cin "Range of Y"
				lab var _ipw "Inverse probability of treatment weights"
			* Rename the elements variables
				foreach var of varlist d1 d0 QAW Q1W Q0W Q1star Qa1star Q0star Qa0star ATE IC  Y A  POM1 POM0 ps {
					rename `var' _`var'
				}
				keep rowid _d1 _d0 _QAW _Q1W _Q0W _Q1star _Qa1star _Q0star _Qa0star _ATE _IC _Y _A _POM1 _POM0 _ps _ipw
				qui save "elementsdata.dta", replace
		}


restore
		
	// Drop the temporary variable for Y
	capture drop ytempvar

	// Clean up
	quietly: rm SLS.R
	// quietly: rm SLS.Rout
	quietly: rm data2.dta
	quietly: rm cvdata.csv
	quietly: rm cvdata.dta
	quietly: rm fulldata.csv
	// quietly: rm .RData
	quietly: memory clean

* Note the completion of the program
	*di "Success"
	



	// Merge in the elements
		if $keepvars == 0 {
			qui drop foldid rowid
		}

	// Merge in the elements
		if $keepvars == 1 {
			qui merge 1:1 rowid using "elementsdata.dta"
		}


end





program cvtmleglsrfcont, rclass

di "Starting cvtmleglsrf... please be patient, additional algorithms are computationally intensive."

preserve

	forval fold = 1/$folds {
		
		if `fold' == 1 _dots 0, title (Performing $folds folds) 
		
		* Load the data to be used for the folds
		qui use "cvdata.dta", clear
		
		* Create a variable indicating the fold number
		capture drop foldnumber
		qui gen foldnumber = `fold'
		
		* Count number of observations
		qui count
		local nobstotal = r(N)
		
		* Export to csv file
		qui: export delimited $varusedforcv foldnumber rowid using "cvdata.csv", nolabel replace
		
		* Rcode for the SuperLearner
		qui: set more off
		qui: file close _all
		qui: file open rcode using SLS.R, write replace
		qui: file write rcode ///
			`"set.seed(123)"' _newline ///
			`"list.of.packages <- c("foreign","SuperLearner","dplyr","randomForest")"' _newline ///
			`"new.packages <- list.of.packages[!(list.of.packages %in% installed.packages()[,"Package"])]"' _newline ///
			`"if(length(new.packages)) install.packages(new.packages, repos='http://cran.us.r-project.org')"' _newline ///
			`"library(SuperLearner)"' _newline ///
			`"library(foreign)"' _newline ///
			`"library(dplyr)"' _newline ///
			`"library(randomForest)"' _newline ///
			`"data <- read.csv("cvdata.csv", sep=",")"' _newline ///
			`"attach(data)"' _newline ///
			`"SL.library <- c("SL.glm","SL.step","SL.glm.interaction","SL.gam","SL.randomForest")"' _newline ///
			`"nvar <- dim(data)[[2]]"' _newline ///
			`"f <- foldnumber[1]"' _newline ///
			`"data <- data %>% select(-c("foldnumber"))"' _newline ///
			`"tdata <- data[foldid!=f,]"' _newline ///		
			`"vdataAa <- data[foldid==f,]"' _newline ///
			`"vdataA1 <- vdataAa; vdataA1[,2] <- 1"' _newline ///
			`"vdataA0 <- vdataAa; vdataA0[,2] <- 0"' _newline ///
			`"Vn <- nrow(vdataAa)"' _newline ///
			`"Tn <- nrow(tdata)"' _newline ///
			`"Y <- tdata[,1]"' _newline ///
			`"A <- tdata[,2]"' _newline ///
			`"X <- tdata %>% select(-c(1,    "foldid", "rowid"))"' _newline ///
			`"W <- tdata %>% select(-c(1, 2, "foldid", "rowid"))"' _newline ///
			`"tset <- try(SuperLearner(Y = tdata[,1], X = X, SL.library=SL.library, family = "binomial", method ="method.NNLS"), silent=TRUE)"' _newline ///
			`"QAW <- predict(tset, newdata=vdataAa)[[1]]"' _newline ///
			`"Q1W <- predict(tset, newdata=vdataA1)[[1]]"' _newline ///
			`"Q0W <- predict(tset, newdata=vdataA0)[[1]]"' _newline ///
			`"QAWt <- predict(tset)[[1]]"' _newline ///
			`"data <- cbind(vdataAa,QAW,Q1W,Q0W)"' _newline ///
			`"write.dta(data, "data2.dta")"' 
		qui: file close rcode

		
		
		* Specify the operating system
	if "`c(os)'" == "MacOSX" {
		qui shell "/usr/local/bin/r" CMD BATCH SLS.R
	}
	else{
	// Write batch file to find R.exe path and R version
		set more off
		qui: file close _all
		qui: file open bat using setup.bat, write replace
		qui: file write bat ///
		`"@echo off"' _newline ///
		`"SET PATHROOT=C:\Program Files\R\"' _newline ///
		`"echo Locating path of R..."' _newline ///
		`"echo."' _newline ///
		`"if not exist "%PATHROOT%" goto:NO_R"' _newline ///
		`":NO_R"' _newline ///
		`"echo R is not installed in your system."' _newline ///
		`"echo."' _newline ///
		`"echo Download it from https://cran.r-project.org/bin/windows/base/"' _newline ///
		`"echo Install it and re-run this script"' _newline ///
		`":DONE"' _newline ///
		`"echo."' _newline ///
		`"pause"'
		qui: file close bat
		
		// 		`"for /f "delims=" %%r in (' dir /b "%PATHROOT%R*" ') do ("' _newline ///
		// 				`"echo Found %%r"' _newline ///
		// 				`"echo ! "%PATHROOT%%%r\bin\x64\R.exe" CMD BATCH SLS.R > runr.do"' _newline ///
		// 				`"echo All set!"' _newline ///
		// 				`"goto:DONE"' _newline ///
		// 		`")"' _newline ///
		
		// 		//Run batch
		// 		! setup.bat
		// 		//Run R
		// 		do runr.do
		
		* Check you have 'rscript' command installed.
		qui: net install rscript, from("https://raw.githubusercontent.com/reifjulian/rscript/master") replace
		
		* Run SLS.R in Windows - 'rscript' command will search for common folders to find R.exe
		qui: rscript using SLS.R
	}
		
	
	* Load the data created by R
		qui clear
		quietly: use "data2.dta", clear
		local folddataname "folddata`fold'"
		qui save "`folddataname'.dta", replace
		
	*if $progress == 1 di "Done fold `fold'"
	
	_dots `fold' 0
		
	}
	
	
* Combine the initial estimate into one data set
	* Append the data sets
	qui use "folddata1.dta", clear
	forval fold = 2/$folds {
		append using "folddata`fold'"
	}

	
* SuperLearner for the propensity score
	* Export to csv file
	qui export delimited using "cvdata.csv", nolabel replace
	
	quietly: rm data2.dta
		
	// Write R Code dependencies: foreign Superlearner
	set more off
	qui: file close _all
	qui: file open rcode using SLS.R, write replace
	qui: file write rcode ///
			`"set.seed(123)"' _newline ///
			`"list.of.packages <- c("foreign","SuperLearner","dplyr","glmnet","randomForest")"' _newline ///
			`"new.packages <- list.of.packages[!(list.of.packages %in% installed.packages()[,"Package"])]"' _newline ///
			`"if(length(new.packages)) install.packages(new.packages, repos='http://cran.us.r-project.org')"' _newline ///
			`"library(SuperLearner)"' _newline ///
			`"library(foreign)"' _newline ///
			`"library(dplyr)"' _newline ///
			`"library(glmnet)"' _newline ///
			`"library(randomForest)"' _newline ///
			`"data <- read.csv("cvdata.csv", sep=",")"' _newline ///
			`"attach(data)"' _newline ///
			`"SL.library <- c("SL.glm","SL.step","SL.glm.interaction","SL.gam","SL.glmnet","SL.randomForest")"' _newline ///
			`"n <- nrow(data)"' _newline ///
			`"dt <- data"' _newline ///
			`"dt <- dt %>% select(-c("foldid"))"' _newline ///
			`"dt <- dt %>% select(-c("rowid"))"' _newline ///
			`"dt <- dt %>% select(-c("QAW"))"' _newline ///
			`"dt <- dt %>% select(-c("Q1W"))"' _newline ///
			`"dt <- dt %>% select(-c("Q0W"))"' _newline ///
			`"nvar <- dim(dt)[[2]]"' _newline ///
			`"Y <- data[,1]"' _newline ///
			`"A <- data[,2]"' _newline ///
			`"X <- data[,2:nvar]"' _newline ///
			`"W <- as.data.frame(data[,3:nvar])"' _newline ///
			`"g <- suppressWarnings(SuperLearner(Y = data[,2], X = W, SL.library = SL.library, family = "binomial", method = "method.NNLS"))"' _newline ///
			`"ps <- g[[4]]"' _newline ///
			`"ps[ps<0.025] <- 0.025"' _newline ///
			`"ps[ps>0.975] <- 0.975"' _newline ///
			`"data <- cbind(data,ps,Y,A)"' _newline ///
			`"write.dta(data, "data2.dta")"'
	qui: file close rcode

	* Specify the operating system
	if "`c(os)'" == "MacOSX" {
		qui shell "/usr/local/bin/r" CMD BATCH SLS.R
	}
	else{
	// Write batch file to find R.exe path and R version
		set more off
		qui: file close _all
		qui: file open bat using setup.bat, write replace
		qui: file write bat ///
		`"@echo off"' _newline ///
		`"SET PATHROOT=C:\Program Files\R\"' _newline ///
		`"echo Locating path of R..."' _newline ///
		`"echo."' _newline ///
		`"if not exist "%PATHROOT%" goto:NO_R"' _newline ///
		`":NO_R"' _newline ///
		`"echo R is not installed in your system."' _newline ///
		`"echo."' _newline ///
		`"echo Download it from https://cran.r-project.org/bin/windows/base/"' _newline ///
		`"echo Install it and re-run this script"' _newline ///
		`":DONE"' _newline ///
		`"echo."' _newline ///
		`"pause"'
		qui: file close bat
		
		// 		`"for /f "delims=" %%r in (' dir /b "%PATHROOT%R*" ') do ("' _newline ///
		// 				`"echo Found %%r"' _newline ///
		// 				`"echo ! "%PATHROOT%%%r\bin\x64\R.exe" CMD BATCH SLS.R > runr.do"' _newline ///
		// 				`"echo All set!"' _newline ///
		// 				`"goto:DONE"' _newline ///
		// 		`")"' _newline ///
		
		// 		//Run batch
		// 		! setup.bat
		// 		//Run R
		// 		do runr.do
		
		* Check you have 'rscript' command installed.
		qui: net install rscript, from("https://raw.githubusercontent.com/reifjulian/rscript/master") replace
		
		* Run SLS.R in Windows - 'rscript' command will search for common folders to find R.exe
		qui: rscript using SLS.R
	}
		
		
// Read Revised Data Back to Stata
	qui clear
	quietly: use "data2.dta", clear
	*import delimited "data2.csv", case(preserve) clear 
	qui cap drop X__000000
	tempvar logQAW logQ1W logQ0W HAW H1W H0W eps1 eps2 eps ATE ICrr ICor

// Q to logit scale
	gen `logQAW' = log(QAW / (1 - QAW))
	gen `logQ1W' = log(Q1W / (1 - Q1W))
	gen `logQ0W' = log(Q0W / (1 - Q0W))

// Clever covariate HAW
	gen  `HAW' = (A / ps) - ((1 - A) / (1 - ps))
	gen  `H1W' = A / ps
	gen  `H0W' = (1 - A) / (1 - ps)

// Estimation of the substitution parameter (Epsilon)
	qui glm Y `H1W' `H0W', fam(binomial) offset(`logQAW') robust noconstant
	mat a= e(b)
	gen `eps1' = a[1,1]
	gen `eps2' = a[1,2]

	qui glm Y `HAW', fam(binomial) offset(`logQAW') robust noconstant
	mat a= e(b)
	gen `eps' = a[1,1]
	
	
// Targeted ATE, update from Q̅^0 (A,W) to Q̅^1 (A,W)
	gen double Qa0star = exp(`H0W'*`eps' + `logQ0W')/(1 + exp(`H0W'*`eps' + `logQ0W'))
	gen double Qa1star = exp(`H1W'*`eps' + `logQ1W')/(1 + exp(`H1W'*`eps' + `logQ1W'))

	gen double Q0star = exp(`H0W'*`eps2' + `logQ0W')/(1 + exp(`H0W'*`eps2' + `logQ0W'))
	gen double Q1star = exp(`H1W'*`eps1' + `logQ1W')/(1 + exp(`H1W'*`eps1' + `logQ1W'))

	// gen double cin = ($b - $a)

	gen double POM1 = cond($flag == 1, Qa1star, (Qa1star * ($b - $a )) + $a , .)
	gen double POM0 = cond($flag == 1, Qa0star, (Qa0star * ($b - $a )) + $a , .)

	di as text " "
	sum POM1 POM0 ps
	di as text " "

// Estimating the updated targeted ATE binary outcome
	gen double ATE = cond($flag == 1, (Qa1star - Qa0star), ((Qa1star * ($b - $a )) + $a ) - ((Qa0star * ($b - $a )) + $a ) , .)
	qui sum ATE
	return scalar ATEtmle = r(mean)

// Relative risk
	qui sum Q1star
	local Q1 = r(mean)
	qui sum Q0star
	local Q0 = r(mean)

// Relative risk and Odds ratio
	local RRtmle = `Q1'/`Q0'
	local logRRtmle = log(`Q1') - log(`Q0')
	local ORtmle = (`Q1' * (1 - `Q0')) / ((1 - `Q1') * `Q0')


// Statistical inference (Efficient Influence Curve)
	gen d1 = cond($flag == 1,(A * (Y - Q1star) / ps) + Q1star - `Q1',(A * (Y - Qa1star) / ps) + Qa1star - `Q1' ,.)
	gen d0 = cond($flag == 1,(1 - A) * (Y - Q0star) / (1 - ps) + Q0star - `Q0',(1 - A) * (Y - Qa0star) / (1 - ps) + Qa0star - `Q0' ,.)
	gen IC = cond($flag == 1,(d1 - d0), ((d1 * ($b - $a )) + $a ) - ((d0 * ($b - $a )) + $a ) , .)
	qui sum IC
	return scalar ATE_SE_tmle = sqrt(r(Var)/r(N))

// Statistical inference ATE
	return scalar ATE_pvalue =  2 * (normalden(abs(return(ATEtmle) / (return(ATE_SE_tmle)))))
	return scalar ATE_LCIa   =  return(ATEtmle) - 1.96 * return(ATE_SE_tmle)
	return scalar ATE_UCIa   =  return(ATEtmle) + 1.96 * return(ATE_SE_tmle)

// Statistical inference RR
	gen `ICrr' = (1/`Q1' * d1) + ((1/`Q0') * d0)
	qui sum `ICrr'
	local varICrr = r(Var)/r(N)
	
	local LCIrr =  exp(`logRRtmle' - 1.96 * sqrt(`varICrr'))
	local UCIrr =  exp(`logRRtmle' + 1.96 * sqrt(`varICrr'))
	
// Statistical inference OR (We do not use the log of the OR)
	local ORtmle = (`Q1' * (1 - `Q0')) / ((1 - `Q1') * `Q0')
	gen `ICor' = (1/ (`Q1' *(1 - `Q1')) * d1) - (1/ (`Q0' *(1 - `Q0')) * d0) // Partial derivatives (use wolfram to check: d/dx log(x/(1-x))+ log(y/(1-y)))
	qui sum `ICor'
	local varICor = r(Var)/r(N)
	local LCIOr = exp(log(`ORtmle') - 1.96 * sqrt(`varICor'))
	local UCIOr = exp(log(`ORtmle') + 1.96 * sqrt(`varICor'))
	*di sqrt(`varICor')
	
// Display Results of ATE
	return scalar CRR = `RRtmle'
	return scalar SE_log_CRR  = sqrt(`varICrr')
	return scalar MOR = `ORtmle'
	return scalar SE_log_MOR  = sqrt(`varICor')

	if $flag==1 {
	disp as text "{hline 63}"
	di "         {c |}" "    ATE         SE     P-value           95% CI"
	disp as text "{hline 63}"
	disp as text "TMLE:    {c |}" %7.4f as result return(ATEtmle) "    " %7.4f as result return(ATE_SE_tmle) "     " %7.4f as result return(ATE_pvalue) as text "     (" %7.4f as result return(ATE_LCIa) as text ","  %7.4f as result return(ATE_UCIa) as text " )"
	disp as text "{hline 63}"
	disp as text " "
	}
	else if $flag!=1{
	disp as text "{hline 63}"
	di "         {c |}" "    ATE         SE     P-value           95% CI"
	disp as text "{hline 63}"
	disp as text "TMLE:    {c |}" %7.1f as result return(ATEtmle) "    " %7.1f as result return(ATE_SE_tmle) "     " %7.4f as result return(ATE_pvalue) as text "     (" %7.1f as result return(ATE_LCIa) as text ","  %7.1f as result return(ATE_UCIa) as text " )"
	disp as text "{hline 63}"
	disp as text " "
	}

	local rrbin ""CRR: "%4.2f `RRtmle'  "; 95%CI:("%3.2f `LCIrr' ", "%3.2f `UCIrr' ")""
	local orbin ""MOR: "%4.2f `ORtmle'  "; 95%CI:("%3.2f `LCIOr' ", "%3.2f `UCIOr' ")""

	* Display the results of CRR and MOR
	disp as text "{hline 51}"
	di "                           Estimate          95% CI"
	disp as text "{hline 51}"
	disp as text "Causal Risk Ratio:      " "{c |}      "  %04.2f as result `RRtmle' as text "     (" %03.2f as result `LCIrr' as text ","  %03.2f as result `UCIrr' as text ")"
	disp as text "Marginal Odds Ratio:    " "{c |}      "  %04.2f as result `ORtmle' as text "     (" %03.2f as result `LCIOr' as text "," %03.2f as result `UCIOr' as text ")"
	disp as text "{hline 51}"

	
// Display covariate balance table

	if $covbalancetable == 0 {
		
			* Create macros from the varlist
				tokenize $variablelist
				local outcome `1'
				macro shift
				local exposure `1'
				macro shift
				local varlist `*'

			* Create the IPW weights
				*di "Capture dropping _ipw 1"
				qui capture drop _ipw
				qui gen _ipw = .
				qui replace _ipw = (`exposure'==1) / ps if `exposure'==1
				qui replace _ipw = (`exposure'==0) / (1- ps) if `exposure'==0
	}


 if $covbalancetable == 1 {
 	
		* Create macros from the varlist
			tokenize $variablelist
			local outcome `1'
			macro shift
			local exposure `1'
			macro shift
			local varlist `*'

		* Create the IPW weights
			capture drop _ipw
			qui gen _ipw = .
			qui replace _ipw = (`exposure'==1) / ps if `exposure'==1
			qui replace _ipw = (`exposure'==0) / (1- ps) if `exposure'==0

		* Layout of the results
			di as text "{hline 67}"
			di as text "                 Standardised Differences            Variance ratio"
			di as text "                          Raw    Weighted           Raw    Weighted"
			di as text "{hline 67}"


		* Calculate the covariate balance
			foreach var in `varlist' {
					di as text "`var'"

					* Raw SMD
					qui summarize `var' if `exposure'==1
					local m1 = r(mean)
					local v1 = r(Var)
					qui summarize `var' if `exposure'==0
					local m0 = r(mean)
					local v0 = r(Var)
					* Calculate the Standardised "mean difference"
					local rSMD = (`m1' - `m0') / sqrt( (`v1' + `v0') /2 )

					* Weighted SMD
					qui summarize `var' [iw=_ipw] if `exposure'==1
					local m1 = r(mean)
					local v1 = r(Var)
					qui summarize `var' [iw=_ipw] if `exposure'==0
					local m0 = r(mean)
					local v0 = r(Var)
					* Calculate the Standardised "mean difference"
					local wSMD = (`m1' - `m0') / sqrt( (`v1' + `v0') /2 )

					* Raw VR
					qui sum `var' if `exposure' ==1
					local v1 = r(Var)
					qui sum `var' if `exposure' ==0
					local v0 = r(Var)
					* Calculate the variance ratio
					local rVR = `v1' / `v0'

					* Weighted VR
					qui sum `var' [iw=_ipw] if `exposure' ==1
					local v1 = r(Var)
					qui sum `var' [iw=_ipw] if `exposure' ==0
					local v0 = r(Var)
					* Calculate the variance ratio
					local wVR `v1' / `v0'

					* Display the values
					di as text "                    " %9.7g as result `rSMD' as text "   " %9.7g as result `wSMD' as text "     " %9.7g as result `rVR' as text "   " %9.7g as result `wVR'
			}
			di as text "{hline 67}"
	
	}
	

// Create density plot to check positivity assumption
	if $covbalancetable == 1 {
	*preserve
	*sort A
	*qui by A: summarize ps
	qui kdensity ps if A==1, generate(x1pointsa d1A) nograph n(10000)
	qui kdensity ps if A==0, generate(x0pointsa d0A) nograph n(10000)
	label variable d1A "A = 1"
	label variable d0A "A = 0"
	twoway (line d0A x0pointsa) || ///
			(line d1A x1pointsa), ///
			xtitle("Propensity score") ///
			ytitle("Density") ///
			graphregion(color(white)) bgcolor(white) plotregion(fcolor(white)) ///
			legend(order(2 "Treated" 1 "Not treated") position(11) cols(1) ring(0) ///
				region(style(none)))
	*restore
	}
	
	
// Drop the variables if the elements option is not specified
		if $keepvars == 0 {
			drop d1 d0
			drop QAW Q1W Q0W
			drop Q1star Qa1star Q0star Qa0star
			drop ATE IC
			rename Y _Y
			rename A _A
			drop _Y _A
			drop _ipw
			drop POM1 POM0 ps
			*drop cin
			drop foldid rowid
			capture drop x1pointsa x0pointsa d1A d0A
		}

// Rename and label the variables if the elements option *is* specified
		if $keepvars == 1 {
			
			capture drop x1pointsa x0pointsa d1A d0A
			
			* Label the variables
				lab var d1 "Parameter for the influence curve"
				lab var d0 "Parameter for the influence curve"
				lab var QAW "Initial prediction of the outcome"
				lab var Q1W "Initial prediction of the outcome for A = 1"
				lab var Q0W "Initial prediction of the outcome for A = 0"
				lab var Q1star "Update of the initial prediction for A = 1"
				lab var Qa1star "Update of the initial prediction for A = 1"
				lab var Q0star "Update of the initial prediction for A = 0"
				lab var Qa0star "Update of the initial prediction for A = 0"
				lab var A "Exposure/Treatment"
				lab var Y "Outcome"
				lab var ATE "Average Treatment Effect"
				lab var IC "Influence Curve"
				lab var Q1star "Update of initial plug-in estimate for A=1"
				lab var Q0star "Update of initial plug-in estimate for A=0"
				lab var POM1 "Potential Outcome Y(1)"
				lab var POM0 "Potential Otucome Y(0)"
				lab var ps "Propensity Score"
				*lab var cin "Range of Y"
				lab var _ipw "Inverse probability of treatment weights"
			* Rename the elements variables
				foreach var of varlist d1 d0 QAW Q1W Q0W Q1star Qa1star Q0star Qa0star ATE IC  Y A  POM1 POM0 ps {
					rename `var' _`var'
				}
				keep rowid _d1 _d0 _QAW _Q1W _Q0W _Q1star _Qa1star _Q0star _Qa0star _ATE _IC _Y _A _POM1 _POM0 _ps _ipw
				qui save "elementsdata.dta", replace
		}


restore
		
	// Drop the temporary variable for Y
	capture drop ytempvar

	// Clean up
	quietly: rm SLS.R
	// quietly: rm SLS.Rout
	quietly: rm data2.dta
	quietly: rm cvdata.csv
	quietly: rm cvdata.dta
	quietly: rm fulldata.csv
	// quietly: rm .RData
	quietly: memory clean

* Note the completion of the program
	*di "Success"
	



	// Merge in the elements
		if $keepvars == 0 {
			qui drop foldid rowid
		}

	// Merge in the elements
		if $keepvars == 1 {
			qui merge 1:1 rowid using "elementsdata.dta"
		}


end






program cvtmlebgam, rclass

di "Starting cvtmlebgam... please be patient, additional algorithms are computationally intensive."

preserve

	forval fold = 1/$folds {
		
		if `fold' == 1 _dots 0, title (Performing $folds folds) 
		
		* Load the data to be used for the folds
		qui use "cvdata.dta", clear
		
		* Create a variable indicating the fold number
		capture drop foldnumber
		qui gen foldnumber = `fold'
		
		* Count number of observations
		qui count
		local nobstotal = r(N)
		
		* Export to csv file
		qui: export delimited $varusedforcv foldnumber rowid using "cvdata.csv", nolabel replace
		
		* Rcode for the SuperLearner
		qui: set more off
		qui: file close _all
		qui: file open rcode using SLS.R, write replace
		qui: file write rcode ///
			`"set.seed(123)"' _newline ///
			`"list.of.packages <- c("foreign","SuperLearner","dplyr","gam","arm")"' _newline ///
			`"new.packages <- list.of.packages[!(list.of.packages %in% installed.packages()[,"Package"])]"' _newline ///
			`"if(length(new.packages)) install.packages(new.packages, repos='http://cran.us.r-project.org')"' _newline ///
			`"library(SuperLearner)"' _newline ///
			`"library(foreign)"' _newline ///
			`"library(dplyr)"' _newline ///
			`"library(gam)"' _newline ///			
			`"library(arm)"' _newline ///
			`"data <- read.csv("cvdata.csv", sep=",")"' _newline ///
			`"attach(data)"' _newline ///
			`"SL.library <- c("SL.glm","SL.step","SL.glm.interaction","SL.gam","SL.bayesglm")"' _newline ///
			`"nvar <- dim(data)[[2]]"' _newline ///
			`"f <- foldnumber[1]"' _newline ///
			`"data <- subset(data, select = -(foldnumber))"' _newline ///
			`"tdata <- data[foldid!=f,]"' _newline ///		
			`"vdataAa <- data[foldid==f,]"' _newline ///
			`"vdataA1 <- vdataAa; vdataA1[,2] <- 1"' _newline ///
			`"vdataA0 <- vdataAa; vdataA0[,2] <- 0"' _newline ///
			`"Vn <- nrow(vdataAa)"' _newline ///
			`"Tn <- nrow(tdata)"' _newline ///
			`"Y <- tdata[,1]"' _newline ///
			`"A <- tdata[,2]"' _newline ///
			`"X <- subset(tdata, select = -c(1, foldid, rowid))"' _newline ///
			`"W <- subset(tdata, select = -c(1, 2, foldid, rowid))"' _newline ///
			`"tset <- try(SuperLearner(Y = tdata[,1], X = X, SL.library=SL.library, family = "binomial", method ="method.NNLS"), silent=TRUE)"' _newline ///
			`"QAW <- predict(tset, newdata=vdataAa)[[1]]"' _newline ///
			`"Q1W <- predict(tset, newdata=vdataA1)[[1]]"' _newline ///
			`"Q0W <- predict(tset, newdata=vdataA0)[[1]]"' _newline ///
			`"QAWt <- predict(tset)[[1]]"' _newline ///
			`"data <- cbind(vdataAa,QAW,Q1W,Q0W)"' _newline ///
			`"write.dta(data, "data2.dta")"' 
		qui: file close rcode

		
		
		* Specify the operating system
	if "`c(os)'" == "MacOSX" {
		qui shell "/usr/local/bin/r" CMD BATCH SLS.R
	}
	else{
	// Write batch file to find R.exe path and R version
		set more off
		qui: file close _all
		qui: file open bat using setup.bat, write replace
		qui: file write bat ///
		`"@echo off"' _newline ///
		`"SET PATHROOT=C:\Program Files\R\"' _newline ///
		`"echo Locating path of R..."' _newline ///
		`"echo."' _newline ///
		`"if not exist "%PATHROOT%" goto:NO_R"' _newline ///
		`":NO_R"' _newline ///
		`"echo R is not installed in your system."' _newline ///
		`"echo."' _newline ///
		`"echo Download it from https://cran.r-project.org/bin/windows/base/"' _newline ///
		`"echo Install it and re-run this script"' _newline ///
		`":DONE"' _newline ///
		`"echo."' _newline ///
		`"pause"'
		qui: file close bat
		
		// 		`"for /f "delims=" %%r in (' dir /b "%PATHROOT%R*" ') do ("' _newline ///
		// 				`"echo Found %%r"' _newline ///
		// 				`"echo ! "%PATHROOT%%%r\bin\x64\R.exe" CMD BATCH SLS.R > runr.do"' _newline ///
		// 				`"echo All set!"' _newline ///
		// 				`"goto:DONE"' _newline ///
		// 		`")"' _newline ///
		
		// 		//Run batch
		// 		! setup.bat
		// 		//Run R
		// 		do runr.do
		
		* Check you have 'rscript' command installed.
		qui: net install rscript, from("https://raw.githubusercontent.com/reifjulian/rscript/master") replace
		
		* Run SLS.R in Windows - 'rscript' command will search for common folders to find R.exe
		qui: rscript using SLS.R
	}
		
	
	* Load the data created by R
		clear
		quietly: use "data2.dta", clear
		local folddataname "folddata`fold'"
		qui save "`folddataname'.dta", replace
		
	*if $progress == 1 di "Done fold `fold'"
	
	_dots `fold' 0
		
	}
	
	
* Combine the initial estimate into one data set
	* Append the data sets
	qui use "folddata1.dta", clear
	forval fold = 2/$folds {
		qui append using "folddata`fold'"
	}

	
* SuperLearner for the propensity score
	* Export to csv file
	qui export delimited using "cvdata.csv", nolabel replace
	
	quietly: rm data2.dta
		
	// Write R Code dependencies: foreign Superlearner
	set more off
	qui: file close _all
	qui: file open rcode using SLS.R, write replace
	qui: file write rcode ///
			`"set.seed(123)"' _newline ///
			`"list.of.packages <- c("foreign","SuperLearner","dplyr","gam","arm")"' _newline ///
			`"new.packages <- list.of.packages[!(list.of.packages %in% installed.packages()[,"Package"])]"' _newline ///
			`"if(length(new.packages)) install.packages(new.packages, repos='http://cran.us.r-project.org')"' _newline ///
			`"library(SuperLearner)"' _newline ///
			`"library(foreign)"' _newline ///
			`"library(dplyr)"' _newline ///
			`"library(gam)"' _newline ///			
			`"library(arm)"' _newline ///			
			`"data <- read.csv("cvdata.csv", sep=",")"' _newline ///
			`"attach(data)"' _newline ///
			`"SL.library <- c("SL.glm","SL.step","SL.glm.interaction","SL.gam","SL.bayesglm")"' _newline ///
			`"n <- nrow(data)"' _newline ///
			`"dt <- data"' _newline ///
			`"dt\$foldid <- NULL"' _newline ///
			`"dt\$QAW <- NULL"' _newline ///
			`"dt\$Q1W <- NULL"' _newline ///
			`"dt\$Q0W <- NULL"' _newline ///
			`"dt\$rowid <- NULL"' _newline ///
			`"nvar <- dim(dt)[[2]]"' _newline ///
			`"Y <- data[,1]"' _newline ///
			`"A <- data[,2]"' _newline ///
			`"X <- data[,2:nvar]"' _newline ///
			`"W <- as.data.frame(data[,3:nvar])"' _newline ///
			`"g <- suppressWarnings(SuperLearner(Y = data[,2], X = W, SL.library = SL.library, family = "binomial", method = "method.NNLS"))"' _newline ///
			`"ps <- g[[4]]"' _newline ///
			`"ps[ps<0.025] <- 0.025"' _newline ///
			`"ps[ps>0.975] <- 0.975"' _newline ///
			`"data <- cbind(data,ps,Y,A)"' _newline ///
			`"write.dta(data, "data2.dta")"'
	qui: file close rcode

	* Specify the operating system
	if "`c(os)'" == "MacOSX" {
		qui shell "/usr/local/bin/r" CMD BATCH SLS.R
	}
	else{
	// Write batch file to find R.exe path and R version
		set more off
		qui: file close _all
		qui: file open bat using setup.bat, write replace
		qui: file write bat ///
		`"@echo off"' _newline ///
		`"SET PATHROOT=C:\Program Files\R\"' _newline ///
		`"echo Locating path of R..."' _newline ///
		`"echo."' _newline ///
		`"if not exist "%PATHROOT%" goto:NO_R"' _newline ///
		`":NO_R"' _newline ///
		`"echo R is not installed in your system."' _newline ///
		`"echo."' _newline ///
		`"echo Download it from https://cran.r-project.org/bin/windows/base/"' _newline ///
		`"echo Install it and re-run this script"' _newline ///
		`":DONE"' _newline ///
		`"echo."' _newline ///
		`"pause"'
		qui: file close bat
		
		// 		`"for /f "delims=" %%r in (' dir /b "%PATHROOT%R*" ') do ("' _newline ///
		// 				`"echo Found %%r"' _newline ///
		// 				`"echo ! "%PATHROOT%%%r\bin\x64\R.exe" CMD BATCH SLS.R > runr.do"' _newline ///
		// 				`"echo All set!"' _newline ///
		// 				`"goto:DONE"' _newline ///
		// 		`")"' _newline ///
		
		// 		//Run batch
		// 		! setup.bat
		// 		//Run R
		// 		do runr.do
		
		* Check you have 'rscript' command installed.
		qui: net install rscript, from("https://raw.githubusercontent.com/reifjulian/rscript/master") replace
		
		* Run SLS.R in Windows - 'rscript' command will search for common folders to find R.exe
		qui: rscript using SLS.R
	}
		
		
// Read Revised Data Back to Stata
	qui clear
	quietly: use "data2.dta", clear
	*import delimited "data2.csv", case(preserve) clear 
	qui cap drop X__000000
	tempvar logQAW logQ1W logQ0W HAW H1W H0W eps1 eps2 eps ATE ICrr ICor

// Q to logit scale
	gen `logQAW' = log(QAW / (1 - QAW))
	gen `logQ1W' = log(Q1W / (1 - Q1W))
	gen `logQ0W' = log(Q0W / (1 - Q0W))

// Clever covariate HAW
	gen  `HAW' = (A / ps) - ((1 - A) / (1 - ps))
	gen  `H1W' = A / ps
	gen  `H0W' = (1 - A) / (1 - ps)

// Estimation of the substitution parameter (Epsilon)
	qui glm Y `H1W' `H0W', fam(binomial) offset(`logQAW') robust noconstant
	mat a= e(b)
	gen `eps1' = a[1,1]
	gen `eps2' = a[1,2]

	qui glm Y `HAW', fam(binomial) offset(`logQAW') robust noconstant
	mat a= e(b)
	gen `eps' = a[1,1]
	
	
// Targeted ATE, update from Q̅^0 (A,W) to Q̅^1 (A,W)
	gen double Qa0star = exp(`H0W'*`eps' + `logQ0W')/(1 + exp(`H0W'*`eps' + `logQ0W'))
	gen double Qa1star = exp(`H1W'*`eps' + `logQ1W')/(1 + exp(`H1W'*`eps' + `logQ1W'))

	gen double Q0star = exp(`H0W'*`eps2' + `logQ0W')/(1 + exp(`H0W'*`eps2' + `logQ0W'))
	gen double Q1star = exp(`H1W'*`eps1' + `logQ1W')/(1 + exp(`H1W'*`eps1' + `logQ1W'))

	// gen double cin = ($b - $a)

	gen double POM1 = cond($flag == 1, Qa1star, (Qa1star * ($b - $a )) + $a , .)
	gen double POM0 = cond($flag == 1, Qa0star, (Qa0star * ($b - $a )) + $a , .)

	di as text " "
	sum POM1 POM0 ps
	di as text " "

// Estimating the updated targeted ATE binary outcome
	gen double ATE = cond($flag == 1, (Qa1star - Qa0star), ((Qa1star * ($b - $a )) + $a ) - ((Qa0star * ($b - $a )) + $a ) , .)
	qui sum ATE
	return scalar ATEtmle = r(mean)

// Relative risk
	qui sum Q1star
	local Q1 = r(mean)
	qui sum Q0star
	local Q0 = r(mean)

// Relative risk and Odds ratio
	local RRtmle = `Q1'/`Q0'
	local logRRtmle = log(`Q1') - log(`Q0')
	local ORtmle = (`Q1' * (1 - `Q0')) / ((1 - `Q1') * `Q0')


// Statistical inference (Efficient Influence Curve)
	gen d1 = cond($flag == 1,(A * (Y - Q1star) / ps) + Q1star - `Q1',(A * (Y - Qa1star) / ps) + Qa1star - `Q1' ,.)
	gen d0 = cond($flag == 1,(1 - A) * (Y - Q0star) / (1 - ps) + Q0star - `Q0',(1 - A) * (Y - Qa0star) / (1 - ps) + Qa0star - `Q0' ,.)
	gen IC = cond($flag == 1,(d1 - d0), ((d1 * ($b - $a )) + $a ) - ((d0 * ($b - $a )) + $a ) , .)
	qui sum IC
	return scalar ATE_SE_tmle = sqrt(r(Var)/r(N))

// Statistical inference ATE
	return scalar ATE_pvalue =  2 * (normalden(abs(return(ATEtmle) / (return(ATE_SE_tmle)))))
	return scalar ATE_LCIa   =  return(ATEtmle) - 1.96 * return(ATE_SE_tmle)
	return scalar ATE_UCIa   =  return(ATEtmle) + 1.96 * return(ATE_SE_tmle)

// Statistical inference RR
	gen `ICrr' = (1/`Q1' * d1) + ((1/`Q0') * d0)
	qui sum `ICrr'
	local varICrr = r(Var)/r(N)
	
	local LCIrr =  exp(`logRRtmle' - 1.96 * sqrt(`varICrr'))
	local UCIrr =  exp(`logRRtmle' + 1.96 * sqrt(`varICrr'))
	
// Statistical inference OR (We do not use the log of the OR)
	local ORtmle = (`Q1' * (1 - `Q0')) / ((1 - `Q1') * `Q0')
	gen `ICor' = (1/ (`Q1' *(1 - `Q1')) * d1) - (1/ (`Q0' *(1 - `Q0')) * d0) // Partial derivatives (use wolfram to check: d/dx log(x/(1-x))+ log(y/(1-y)))
	qui sum `ICor'
	local varICor = r(Var)/r(N)
	local LCIOr = exp(log(`ORtmle') - 1.96 * sqrt(`varICor'))
	local UCIOr = exp(log(`ORtmle') + 1.96 * sqrt(`varICor'))
	*di sqrt(`varICor')
	
// Display Results of ATE
	return scalar CRR = `RRtmle'
	return scalar SE_log_CRR  = sqrt(`varICrr')
	return scalar MOR = `ORtmle'
	return scalar SE_log_MOR  = sqrt(`varICor')

	if $flag==1 {
	disp as text "{hline 63}"
	di "         {c |}" "    ATE         SE     P-value           95% CI"
	disp as text "{hline 63}"
	disp as text "TMLE:    {c |}" %7.4f as result return(ATEtmle) "    " %7.4f as result return(ATE_SE_tmle) "     " %7.4f as result return(ATE_pvalue) as text "     (" %7.4f as result return(ATE_LCIa) as text ","  %7.4f as result return(ATE_UCIa) as text " )"
	disp as text "{hline 63}"
	disp as text " "
	}
	else if $flag!=1{
	disp as text "{hline 63}"
	di "         {c |}" "    ATE         SE     P-value           95% CI"
	disp as text "{hline 63}"
	disp as text "TMLE:    {c |}" %7.1f as result return(ATEtmle) "    " %7.1f as result return(ATE_SE_tmle) "     " %7.4f as result return(ATE_pvalue) as text "     (" %7.1f as result return(ATE_LCIa) as text ","  %7.1f as result return(ATE_UCIa) as text " )"
	disp as text "{hline 63}"
	disp as text " "
	}

	local rrbin ""CRR: "%4.2f `RRtmle'  "; 95%CI:("%3.2f `LCIrr' ", "%3.2f `UCIrr' ")""
	local orbin ""MOR: "%4.2f `ORtmle'  "; 95%CI:("%3.2f `LCIOr' ", "%3.2f `UCIOr' ")""

	* Display the results of CRR and MOR
	disp as text "{hline 51}"
	di "                           Estimate          95% CI"
	disp as text "{hline 51}"
	disp as text "Causal Risk Ratio:      " "{c |}      "  %04.2f as result `RRtmle' as text "     (" %03.2f as result `LCIrr' as text ","  %03.2f as result `UCIrr' as text ")"
	disp as text "Marginal Odds Ratio:    " "{c |}      "  %04.2f as result `ORtmle' as text "     (" %03.2f as result `LCIOr' as text "," %03.2f as result `UCIOr' as text ")"
	disp as text "{hline 51}"

	
// Display covariate balance table

	if $covbalancetable == 0 {
		
			* Create macros from the varlist
				tokenize $variablelist
				local outcome `1'
				macro shift
				local exposure `1'
				macro shift
				local varlist `*'

			* Create the IPW weights
				*di "Capture dropping _ipw 1"
				qui capture drop _ipw
				qui gen _ipw = .
				qui replace _ipw = (`exposure'==1) / ps if `exposure'==1
				qui replace _ipw = (`exposure'==0) / (1- ps) if `exposure'==0
	}


 if $covbalancetable == 1 {
 	
		* Create macros from the varlist
			tokenize $variablelist
			local outcome `1'
			macro shift
			local exposure `1'
			macro shift
			local varlist `*'

		* Create the IPW weights
			capture drop _ipw
			qui gen _ipw = .
			qui replace _ipw = (`exposure'==1) / ps if `exposure'==1
			qui replace _ipw = (`exposure'==0) / (1- ps) if `exposure'==0

		* Layout of the results
			di as text "{hline 67}"
			di as text "                 Standardised Differences            Variance ratio"
			di as text "                          Raw    Weighted           Raw    Weighted"
			di as text "{hline 67}"


		* Calculate the covariate balance
			foreach var in `varlist' {
					di as text "`var'"

					* Raw SMD
					qui summarize `var' if `exposure'==1
					local m1 = r(mean)
					local v1 = r(Var)
					qui summarize `var' if `exposure'==0
					local m0 = r(mean)
					local v0 = r(Var)
					* Calculate the Standardised "mean difference"
					local rSMD = (`m1' - `m0') / sqrt( (`v1' + `v0') /2 )

					* Weighted SMD
					qui summarize `var' [iw=_ipw] if `exposure'==1
					local m1 = r(mean)
					local v1 = r(Var)
					qui summarize `var' [iw=_ipw] if `exposure'==0
					local m0 = r(mean)
					local v0 = r(Var)
					* Calculate the Standardised "mean difference"
					local wSMD = (`m1' - `m0') / sqrt( (`v1' + `v0') /2 )

					* Raw VR
					qui sum `var' if `exposure' ==1
					local v1 = r(Var)
					qui sum `var' if `exposure' ==0
					local v0 = r(Var)
					* Calculate the variance ratio
					local rVR = `v1' / `v0'

					* Weighted VR
					qui sum `var' [iw=_ipw] if `exposure' ==1
					local v1 = r(Var)
					qui sum `var' [iw=_ipw] if `exposure' ==0
					local v0 = r(Var)
					* Calculate the variance ratio
					local wVR `v1' / `v0'

					* Display the values
					di as text "                    " %9.7g as result `rSMD' as text "   " %9.7g as result `wSMD' as text "     " %9.7g as result `rVR' as text "   " %9.7g as result `wVR'
			}
			di as text "{hline 67}"
	
	}
	

// Create density plot to check positivity assumption
	if $covbalancetable == 1 {
	*preserve
	*sort A
	*qui by A: summarize ps
	qui kdensity ps if A==1, generate(x1pointsa d1A) nograph n(10000)
	qui kdensity ps if A==0, generate(x0pointsa d0A) nograph n(10000)
	label variable d1A "A = 1"
	label variable d0A "A = 0"
	twoway (line d0A x0pointsa) || ///
			(line d1A x1pointsa), ///
			xtitle("Propensity score") ///
			ytitle("Density") ///
			graphregion(color(white)) bgcolor(white) plotregion(fcolor(white)) ///
			legend(order(2 "Treated" 1 "Not treated") position(11) cols(1) ring(0) ///
				region(style(none)))
	*restore
	}
	
	
// Drop the variables if the elements option is not specified
		if $keepvars == 0 {
			drop d1 d0
			drop QAW Q1W Q0W
			drop Q1star Qa1star Q0star Qa0star
			drop ATE IC
			rename Y _Y
			rename A _A
			drop _Y _A
			drop _ipw
			drop POM1 POM0 ps
			*drop cin
			drop foldid rowid
			capture drop x1pointsa x0pointsa d1A d0A
		}

// Rename and label the variables if the elements option *is* specified
		if $keepvars == 1 {
			
			capture drop x1pointsa x0pointsa d1A d0A
			
			* Label the variables
				lab var d1 "Parameter for the influence curve"
				lab var d0 "Parameter for the influence curve"
				lab var QAW "Initial prediction of the outcome"
				lab var Q1W "Initial prediction of the outcome for A = 1"
				lab var Q0W "Initial prediction of the outcome for A = 0"
				lab var Q1star "Update of the initial prediction for A = 1"
				lab var Qa1star "Update of the initial prediction for A = 1"
				lab var Q0star "Update of the initial prediction for A = 0"
				lab var Qa0star "Update of the initial prediction for A = 0"
				lab var A "Exposure/Treatment"
				lab var Y "Outcome"
				lab var ATE "Average Treatment Effect"
				lab var IC "Influence Curve"
				lab var Q1star "Update of initial plug-in estimate for A=1"
				lab var Q0star "Update of initial plug-in estimate for A=0"
				lab var POM1 "Potential Outcome Y(1)"
				lab var POM0 "Potential Otucome Y(0)"
				lab var ps "Propensity Score"
				*lab var cin "Range of Y"
				lab var _ipw "Inverse probability of treatment weights"
			* Rename the elements variables
				foreach var of varlist d1 d0 QAW Q1W Q0W Q1star Qa1star Q0star Qa0star ATE IC  Y A  POM1 POM0 ps {
					rename `var' _`var'
				}
				keep rowid _d1 _d0 _QAW _Q1W _Q0W _Q1star _Qa1star _Q0star _Qa0star _ATE _IC _Y _A _POM1 _POM0 _ps _ipw
				qui save "elementsdata.dta", replace
		}


		
restore

	// Drop the temporary variable for Y
	capture drop ytempvar

	// Clean up
	quietly: rm SLS.R
	// quietly: rm SLS.Rout
	quietly: rm data2.dta
	quietly: rm cvdata.csv
	quietly: rm cvdata.dta
	quietly: rm fulldata.csv
	// quietly: rm .RData
	quietly: memory clean

* Note the completion of the program
	*di "Success"
	



	// Merge in the elements
		if $keepvars == 0 {
			qui drop foldid rowid
		}

	// Merge in the elements
		if $keepvars == 1 {
			qui merge 1:1 rowid using "elementsdata.dta"
		}


end






program tmlebal, rclass
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
        `"fulldata <- read.csv("fulldata.csv", sep=",")"' _newline ///
        `"attach(data)"' _newline ///
        `"SL.library <- c("SL.glm","SL.step","SL.glm.interaction")"' _newline ///
        `"n <- nrow(data)"' _newline ///
        `"nvar <- dim(data)[[2]]"' _newline ///
        `"Y <- data[,1]"' _newline ///
        `"A <- data[,2]"' _newline ///
        `"X <- data[,2:nvar]"' _newline ///
        `"W <- as.data.frame(data[,3:nvar])"' _newline ///
        `"X1 <- X0 <- X"' _newline ///
        `"X1[,1] <- 1"' _newline ///
        `"X0[,1] <- 0"' _newline ///
        `"newdata <- rbind(X,X1,X0)"' _newline ///
        `"Q <- try(SuperLearner(Y = data[,1] ,X = X, SL.library=SL.library, family = "binomial", newX=newdata, method ="method.NNLS"), silent=TRUE)"' _newline ///
        `"Q <- as.data.frame(Q[[4]])"' _newline ///
        `"QAW <- Q[1:n,]"' _newline ///
        `"Q1W <- Q[((n+1):(2*n)),]"' _newline ///
        `"Q0W <- Q[((2*n+1):(3*n)),]"' _newline ///
        `"g <- suppressWarnings(SuperLearner(Y = data[,2], X = W, SL.library = SL.library, family = "binomial", method = "method.NNLS"))"' _newline ///
        `"ps <- g[[4]]"' _newline ///
        `"ps[ps<0.025] <- 0.025"' _newline ///
        `"ps[ps>0.975] <- 0.975"' _newline ///
        `"data <- cbind(fulldata,QAW,Q1W,Q0W,ps,Y,A)"' _newline ///
        `"write.dta(data, "data2.dta")"'
qui: file close rcode


	if "`c(os)'" == "MacOSX" {
		qui shell "/usr/local/bin/r" CMD BATCH SLS.R
	}
	else{
	// Write batch file to find R.exe path and R version
		set more off
		qui: file close _all
		qui: file open bat using setup.bat, write replace
		qui: file write bat ///
		`"@echo off"' _newline ///
		`"SET PATHROOT=C:\Program Files\R\"' _newline ///
		`"echo Locating path of R..."' _newline ///
		`"echo."' _newline ///
		`"if not exist "%PATHROOT%" goto:NO_R"' _newline ///
		`":NO_R"' _newline ///
		`"echo R is not installed in your system."' _newline ///
		`"echo."' _newline ///
		`"echo Download it from https://cran.r-project.org/bin/windows/base/"' _newline ///
		`"echo Install it and re-run this script"' _newline ///
		`":DONE"' _newline ///
		`"echo."' _newline ///
		`"pause"'
		qui: file close bat
		
		// 		`"for /f "delims=" %%r in (' dir /b "%PATHROOT%R*" ') do ("' _newline ///
		// 				`"echo Found %%r"' _newline ///
		// 				`"echo ! "%PATHROOT%%%r\bin\x64\R.exe" CMD BATCH SLS.R > runr.do"' _newline ///
		// 				`"echo All set!"' _newline ///
		// 				`"goto:DONE"' _newline ///
		// 		`")"' _newline ///
		
		// 		//Run batch
		// 		! setup.bat
		// 		//Run R
		// 		do runr.do
		
		* Check you have 'rscript' command installed.
		qui: net install rscript, from("https://raw.githubusercontent.com/reifjulian/rscript/master") replace
		
		* Run SLS.R in Windows - 'rscript' command will search for common folders to find R.exe
		qui: rscript using SLS.R
	}

	
// Read Revised Data Back to Stata
	qui clear
	quietly: use "data2.dta", clear
	qui cap drop X__000000
	tempvar logQAW logQ1W logQ0W HAW H1W H0W eps1 eps2 eps ATE ICrr ICor

// Create density plot to check positivity assumption
	preserve
	sort A
	qui by A: summarize ps
	qui kdensity ps if A==1, generate(x1pointsa d1A) nograph n(10000)
	qui kdensity ps if A==0, generate(x0pointsa d0A) nograph n(10000)
	label variable d1A "A = 1"
	label variable d0A "A = 0"
	twoway (line d0A x0pointsa) || ///
			(line d1A x1pointsa), ///
			xtitle("Propensity score") ///
			ytitle("Density") ///
			graphregion(color(white)) bgcolor(white) plotregion(fcolor(white)) ///
			legend(order(2 "Treated" 1 "Not treated") position(11) cols(1) ring(0) ///
				region(style(none)))
	restore

// Q to logit scale
	gen `logQAW' = log(QAW / (1 - QAW))
	gen `logQ1W' = log(Q1W / (1 - Q1W))
	gen `logQ0W' = log(Q0W / (1 - Q0W))

// Clever covariate HAW
	gen  `HAW' = (A / ps) - ((1 - A) / (1 - ps))
	gen  `H1W' = A / ps
	gen  `H0W' = (1 - A) / (1 - ps)

// Estimation of the substitution parameter (Epsilon)
	qui glm Y `H1W' `H0W', fam(binomial) offset(`logQAW') robust noconstant
	mat a= e(b)
	gen `eps1' = a[1,1]
	gen `eps2' = a[1,2]

	qui glm Y `HAW', fam(binomial) offset(`logQAW') robust noconstant
	mat a= e(b)
	gen `eps' = a[1,1]


// Targeted ATE, update from Q̅^0 (A,W) to Q̅^1 (A,W)
	gen double Qa0star = exp(`H0W'*`eps' + `logQ0W')/(1 + exp(`H0W'*`eps' + `logQ0W'))
	gen double Qa1star = exp(`H1W'*`eps' + `logQ1W')/(1 + exp(`H1W'*`eps' + `logQ1W'))

	gen double Q0star = exp(`H0W'*`eps2' + `logQ0W')/(1 + exp(`H0W'*`eps2' + `logQ0W'))
	gen double Q1star = exp(`H1W'*`eps1' + `logQ1W')/(1 + exp(`H1W'*`eps1' + `logQ1W'))

	// gen double cin = ($b - $a)

	gen double POM1 = cond($flag == 1, Qa1star, (Qa1star * ($b - $a )) + $a , .)
	gen double POM0 = cond($flag == 1, Qa0star, (Qa0star * ($b - $a )) + $a , .)

	sum POM1 POM0 ps
	di as text " "

// Estimating the updated targeted ATE binary outcome
	gen double ATE = cond($flag == 1, (Qa1star - Qa0star), ((Qa1star * ($b - $a )) + $a ) - ((Qa0star * ($b - $a )) + $a ) , .)
	qui sum ATE
	return scalar ATEtmle = r(mean)

// Relative risk
	qui sum Q1star
	local Q1 = r(mean)
	qui sum Q0star
	local Q0 = r(mean)

// Relative risk and Odds ratio
	local RRtmle = `Q1'/`Q0'
	local logRRtmle = log(`Q1') - log(`Q0')
	local ORtmle = (`Q1' * (1 - `Q0')) / ((1 - `Q1') * `Q0')

// Statistical inference (Efficient Influence Curve)
	gen d1 = cond($flag == 1,(A * (Y - Q1star) / ps) + Q1star - `Q1',(A * (Y - Qa1star) / ps) + Qa1star - `Q1' ,.)
	gen d0 = cond($flag == 1,(1 - A) * (Y - Q0star) / (1 - ps) + Q0star - `Q0',(1 - A) * (Y - Qa0star) / (1 - ps) + Qa0star - `Q0' ,.)
	gen IC = cond($flag == 1,(d1 - d0), ((d1 * ($b - $a )) + $a ) - ((d0 * ($b - $a )) + $a ) , .)
	qui sum IC
	return scalar ATE_SE_tmle = sqrt(r(Var)/r(N))

// Statistical inference ATE
	return scalar ATE_pvalue =  2 * (normalden(abs(return(ATEtmle) / (return(ATE_SE_tmle)))))
	return scalar ATE_LCIa   =  return(ATEtmle) - 1.96 * return(ATE_SE_tmle)
	return scalar ATE_UCIa   =  return(ATEtmle) + 1.96 * return(ATE_SE_tmle)

// Statistical inference RR
	gen `ICrr' = (1/`Q1' * d1) + ((1/`Q0') * d0)
	qui sum `ICrr'
	local varICrr = r(Var)/r(N)
	local LCIrr =  exp(`logRRtmle' - 1.96 * sqrt(`varICrr'))
	local UCIrr =  exp(`logRRtmle' + 1.96 * sqrt(`varICrr'))

// Statistical inference OR
	// 	gen `ICor' = ((1 - `Q0') / `Q0' / (1 - `Q1')^2) * d1 - (`Q1' / (1 - `Q1') / `Q0'^2) * d0
	// 	qui sum `ICor'
	// 	local varICor = r(Var)/r(N)
	// 	local LCIOr =  `ORtmle' - 1.96 * sqrt(`varICor')
	// 	local UCIOr =  `ORtmle' + 1.96 * sqrt(`varICor')
	
	// Statistical inference OR (We do not use the log of the OR)
	*local ORtmle = (`Q1'  / (1 - `Q1')) / (`Q0' / ((1 - `Q0')))
	gen `ICor' = (1/ (`Q1' *(1 - `Q1')) * d1) - (1/ (`Q0' *(1 - `Q0')) * d0) // Partial derivatives (use wolfram to check: d/dx log(x/(1-x))+ log(y/(1-y)))
	qui sum `ICor'
	local varICor = r(Var)/r(N)
	local LCIOr = exp(log(`ORtmle') - 1.96 * sqrt(`varICor'))
	local UCIOr = exp(log(`ORtmle') + 1.96 * sqrt(`varICor'))
	
	
	// Display Results of ATE
		return scalar CRR = `RRtmle'
		return scalar SE_log_CRR  = sqrt(`varICrr')
		return scalar MOR = `ORtmle'
		return scalar SE_log_MOR  = sqrt(`varICor')

		if $flag==1 {
		disp as text "{hline 63}"
		di "         {c |}" "    ATE         SE     P-value           95% CI"
		disp as text "{hline 63}"
		disp as text "TMLE:    {c |}" %7.4f as result return(ATEtmle) "    " %7.4f as result return(ATE_SE_tmle) "     " %7.4f as result return(ATE_pvalue) as text "     (" %7.4f as result return(ATE_LCIa) as text ","  %7.4f as result return(ATE_UCIa) as text " )"
		disp as text "{hline 63}"
		disp as text " "
		}
		else if $flag!=1{
		disp as text "{hline 63}"
		di "         {c |}" "    ATE         SE     P-value           95% CI"
		disp as text "{hline 63}"
		disp as text "TMLE:    {c |}" %7.1f as result return(ATEtmle) "    " %7.1f as result return(ATE_SE_tmle) "     " %7.4f as result return(ATE_pvalue) as text "     (" %7.1f as result return(ATE_LCIa) as text ","  %7.1f as result return(ATE_UCIa) as text " )"
		disp as text "{hline 63}"
		disp as text " "
		}

		local rrbin ""CRR: "%4.2f `RRtmle'  "; 95%CI:("%3.2f `LCIrr' ", "%3.2f `UCIrr' ")""
		local orbin ""MOR: "%4.2f `ORtmle'  "; 95%CI:("%3.2f `LCIOr' ", "%3.2f `UCIOr' ")""

		* Display the results of CRR and MOR
		disp as text "{hline 51}"
		di "                           Estimate          95% CI"
		disp as text "{hline 51}"
		disp as text "Causal Risk Ratio:      " "{c |}      "  %4.2f as result `RRtmle' as text "     (" %3.2f as result `LCIrr' as text ","  %3.2f as result `UCIrr' as text ")"
		disp as text "Marginal Odds Ratio:    " "{c |}      "  %4.2f as result `ORtmle' as text "     (" %3.2f as result `LCIOr' as text "," %3.2f as result `UCIOr' as text ")"
		disp as text "{hline 51}"


		
// Display covariate balance table
		* Create macros from the varlist
			tokenize $variablelist
			local outcome `1'
			macro shift
			local exposure `1'
			macro shift
			local varlist `*'

		* Create the IPW weights
			capture drop _ipw
			qui gen _ipw = .
			qui replace _ipw = (`exposure'==1) / ps if `exposure'==1
			qui replace _ipw = (`exposure'==0) / (1- ps) if `exposure'==0

		* Layout of the results
			di as text "{hline 67}"
			di as text "                 Standardised Differences            Variance ratio"
			di as text "                          Raw    Weighted           Raw    Weighted"
			di as text "{hline 67}"

			
		* Calculate the covariate balance
			foreach var in `varlist' {
					di as text "`var'"

					* Raw SMD
					qui summarize `var' if `exposure'==1
					local m1 = r(mean)
					local v1 = r(Var)
					qui summarize `var' if `exposure'==0
					local m0 = r(mean)
					local v0 = r(Var)
					* Calculate the Standardised "mean difference"
					local rSMD = (`m1' - `m0') / sqrt( (`v1' + `v0') /2 )

					* Weighted SMD
					qui summarize `var' [iw=_ipw] if `exposure'==1
					local m1 = r(mean)
					local v1 = r(Var)
					qui summarize `var' [iw=_ipw] if `exposure'==0
					local m0 = r(mean)
					local v0 = r(Var)
					* Calculate the Standardised "mean difference"
					local wSMD = (`m1' - `m0') / sqrt( (`v1' + `v0') /2 )

					* Raw VR
					qui sum `var' if `exposure' ==1
					local v1 = r(Var)
					qui sum `var' if `exposure' ==0
					local v0 = r(Var)
					* Calculate the variance ratio
					local rVR = `v1' / `v0'

					* Weighted VR
					qui sum `var' [iw=_ipw] if `exposure' ==1
					local v1 = r(Var)
					qui sum `var' [iw=_ipw] if `exposure' ==0
					local v0 = r(Var)
					* Calculate the variance ratio
					local wVR `v1' / `v0'

					* Display the values
					di as text "                    " %9.7g as result `rSMD' as text "   " %9.7g as result `wSMD' as text "     " %9.7g as result `rVR' as text "   " %9.7g as result `wVR'
			}
			di as text "{hline 67}"

	* Drop the variables if the elements option is not specified
		if $keepvars == 0 {
			drop d1 d0
			drop QAW Q1W Q0W
			drop Q1star Qa1star Q0star Qa0star
			drop ATE IC
			drop Y A
			drop POM1 POM0 ps
			*drop cin
			capture drop ytempvar
		}

	* Rename and label the variables if the elements option *is* specified
		if $keepvars == 1 {
			* Label the variables
				lab var d1 "Parameter for the influence curve"
				lab var d0 "Parameter for the influence curve"
				lab var QAW "Initial prediction of the outcome"
				lab var Q1W "Initial prediction of the outcome for A = 1"
				lab var Q0W "Initial prediction of the outcome for A = 0"
				lab var Q1star "Update of the initial prediction for A = 1"
				lab var Qa1star "Update of the initial prediction for A = 1"
				lab var Q0star "Update of the initial prediction for A = 0"
				lab var Qa0star "Update of the initial prediction for A = 0"
				lab var A "Exposure/Treatment"
				lab var Y "Outcome"
				lab var ATE "Average Treatment Effect"
				lab var IC "Influence Curve"
				lab var Q1star "Update of initial plug-in estimate for A=1"
				lab var Q0star "Update of initial plug-in estimate for A=0"
				lab var POM1 "Potential Outcome Y(1)"
				lab var POM0 "Potential Otucome Y(0)"
				lab var ps "Propensity Score"
				*lab var cin "Range of Y"
				capture drop ytempvar
			* Rename the elements variables
				foreach var of varlist d1 d0 QAW Q1W Q0W Q1star Qa1star Q0star Qa0star ATE IC  Y A  POM1 POM0 ps {
					rename `var' _`var'
				}
		}


	// Clean up
	quietly: rm SLS.R
	// quietly: rm SLS.Rout
	quietly: rm data2.dta
	quietly: rm data.csv
	quietly: rm fulldata.csv
	// quietly: rm .RData
	quietly: memory clean

end




program tmlebgam, rclass
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
        `"fulldata <- read.csv("fulldata.csv", sep=",")"' _newline ///
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
        `"Q <- try(SuperLearner(Y = data[,1] ,X = X, SL.library=SL.library, family = "binomial", newX=newdata, method ="method.NNLS"), silent=TRUE)"' _newline ///
        `"Q <- as.data.frame(Q[[4]])"' _newline ///
        `"QAW <- Q[1:n,]"' _newline ///
        `"Q1W <- Q[((n+1):(2*n)),]"' _newline ///
        `"Q0W <- Q[((2*n+1):(3*n)),]"' _newline ///
        `"g <- suppressWarnings(SuperLearner(Y = data[,2], X = W, SL.library = SL.library, family = "binomial", method = "method.NNLS"))"' _newline ///
        `"ps <- g[[4]]"' _newline ///
        `"ps[ps<0.025] <- 0.025"' _newline ///
        `"ps[ps>0.975] <- 0.975"' _newline ///
        `"data <- cbind(fulldata,QAW,Q1W,Q0W,ps,Y,A)"' _newline ///
        `"write.dta(data, "data2.dta")"'
qui: file close rcode


	if "`c(os)'" == "MacOSX" {
	qui shell "/usr/local/bin/r" CMD BATCH SLS.R
	}
	else{
	// Write batch file to find R.exe path and R version
		set more off
		qui: file close _all
		qui: file open bat using setup.bat, write replace
		qui: file write bat ///
		`"@echo off"' _newline ///
		`"SET PATHROOT=C:\Program Files\R\"' _newline ///
		`"echo Locating path of R..."' _newline ///
		`"echo."' _newline ///
		`"if not exist "%PATHROOT%" goto:NO_R"' _newline ///
		`":NO_R"' _newline ///
		`"echo R is not installed in your system."' _newline ///
		`"echo."' _newline ///
		`"echo Download it from https://cran.r-project.org/bin/windows/base/"' _newline ///
		`"echo Install it and re-run this script"' _newline ///
		`":DONE"' _newline ///
		`"echo."' _newline ///
		`"pause"'
		qui: file close bat
		
		// 		`"for /f "delims=" %%r in (' dir /b "%PATHROOT%R*" ') do ("' _newline ///
		// 				`"echo Found %%r"' _newline ///
		// 				`"echo ! "%PATHROOT%%%r\bin\x64\R.exe" CMD BATCH SLS.R > runr.do"' _newline ///
		// 				`"echo All set!"' _newline ///
		// 				`"goto:DONE"' _newline ///
		// 		`")"' _newline ///
		
		// 		//Run batch
		// 		! setup.bat
		// 		//Run R
		// 		do runr.do
		
		* Check you have 'rscript' command installed.
		qui: net install rscript, from("https://raw.githubusercontent.com/reifjulian/rscript/master") replace
		
		* Run SLS.R in Windows - 'rscript' command will search for common folders to find R.exe
		qui: rscript using SLS.R
	}

	
// Read Revised Data Back to Stata
qui clear
quietly: use "data2.dta", clear
qui cap drop X__000000
tempvar logQAW logQ1W logQ0W HAW H1W H0W eps1 eps2 eps ATE ICrr ICor
// Q to logit scale
gen `logQAW' = log(QAW / (1 - QAW))
gen `logQ1W' = log(Q1W / (1 - Q1W))
gen `logQ0W' = log(Q0W / (1 - Q0W))

// Clever covariate HAW
gen  `HAW' = (A / ps) - ((1 - A) / (1 - ps))
gen  `H1W' = A / ps
gen  `H0W' = (1 - A) / (1 - ps)

// Estimation of the substitution parameter (Epsilon)
qui glm Y `H1W' `H0W', fam(binomial) offset(`logQAW') robust noconstant
mat a= e(b)
gen `eps1' = a[1,1]
gen `eps2' = a[1,2]

qui glm Y `HAW', fam(binomial) offset(`logQAW') robust noconstant
mat a= e(b)
gen `eps' = a[1,1]


// Targeted ATE, update from Q̅^0 (A,W) to Q̅^1 (A,W)
gen double Qa0star = exp(`H0W'*`eps' + `logQ0W')/(1 + exp(`H0W'*`eps' + `logQ0W'))
gen double Qa1star = exp(`H1W'*`eps' + `logQ1W')/(1 + exp(`H1W'*`eps' + `logQ1W'))

gen double Q0star = exp(`H0W'*`eps2' + `logQ0W')/(1 + exp(`H0W'*`eps2' + `logQ0W'))
gen double Q1star = exp(`H1W'*`eps1' + `logQ1W')/(1 + exp(`H1W'*`eps1' + `logQ1W'))

// gen double cin = ($b - $a)

gen double POM1 = cond($flag == 1, Qa1star, (Qa1star * ($b - $a )) + $a , .)
gen double POM0 = cond($flag == 1, Qa0star, (Qa0star * ($b - $a )) + $a , .)

summ POM1 POM0 ps
di as text " "

// Estimating the updated targeted ATE binary outcome
gen double ATE = cond($flag == 1, (Qa1star - Qa0star), ((Qa1star * ($b - $a )) + $a ) - ((Qa0star * ($b - $a )) + $a ) , .)
qui sum ATE
return scalar ATEtmle = r(mean)

// Relative risk and Odds ratio
qui sum Q1star
local Q1 = r(mean)
qui sum Q0star
local Q0 = r(mean)

// Relative risk and Odds ratio
	local RRtmle = `Q1'/`Q0'
	local logRRtmle = log(`Q1') - log(`Q0')
	local ORtmle = (`Q1' * (1 - `Q0')) / ((1 - `Q1') * `Q0')
	local logORtmle = log((`Q1' * (1 - `Q0')) / ((1 - `Q1') * `Q0'))

// Statistical inference (Efficient Influence Curve)
	gen d1 = cond($flag == 1,(A * (Y - Q1star) / ps) + Q1star - `Q1',(A * (Y - Qa1star) / ps) + Qa1star - `Q1' ,.)
	gen d0 = cond($flag == 1,(1 - A) * (Y - Q0star) / (1 - ps) + Q0star - `Q0',(1 - A) * (Y - Qa0star) / (1 - ps) + Qa0star - `Q0' ,.)
	gen IC = cond($flag == 1,(d1 - d0), ((d1 * ($b - $a )) + $a ) - ((d0 * ($b - $a )) + $a ) , .)
	qui sum IC
	return scalar ATE_SE_tmle = sqrt(r(Var)/r(N))

// Statistical inference ATE
	return scalar ATE_pvalue =  2 * (normalden(abs(return(ATEtmle) / (return(ATE_SE_tmle)))))
	return scalar ATE_LCIa   =  return(ATEtmle) - 1.96 * return(ATE_SE_tmle)
	return scalar ATE_UCIa   =  return(ATEtmle) + 1.96 * return(ATE_SE_tmle)

// Statistical inference RR
	gen `ICrr' = (1/`Q1' * d1) + ((1/`Q0') * d0)
	qui sum `ICrr'
	local varICrr = r(Var)/r(N)
	local LCIrr =  exp(`logRRtmle' - 1.96 * sqrt(`varICrr'))
	local UCIrr =  exp(`logRRtmle' + 1.96 * sqrt(`varICrr'))
	
// Statistical inference OR (We do not use the log of the OR
    gen `ICor' = 1/ (`Q1' *(1 - `Q1')) * d1 - 1/ (`Q0' *(1 - `Q0')) * d0 // Partial derivatives (use wolfram to check: d/dx log(x/(1-x))+ log(y/(1-y)))
    qui sum `ICor'
    local varICor = r(Var)/r(N)
	local LCIOr = exp(log(`ORtmle') - 1.96 * sqrt(`varICor'))
	local UCIOr = exp(log(`ORtmle') + 1.96 * sqrt(`varICor'))
	
// Display Results of ATE
	return scalar CRR = `RRtmle'
	return scalar SE_log_CRR  = sqrt(`varICrr')
	return scalar MOR = `ORtmle'
	return scalar SE_log_MOR  = sqrt(`varICor')

	if $flag==1 {
	disp as text "{hline 63}"
	di "         {c |}" "    ATE         SE     P-value           95% CI"
	disp as text "{hline 63}"
	disp as text "TMLE:    {c |}" %7.4f as result return(ATEtmle) "    " %7.4f as result return(ATE_SE_tmle) "     " %7.4f as result return(ATE_pvalue) as text "     (" %7.4f as result return(ATE_LCIa) ","  %7.4f as result return(ATE_UCIa) as text " )"
	disp as text "{hline 63}"
	disp as text " "
	}
	else if $flag!=1{
	disp as text "{hline 63}"
	di "         {c |}" "    ATE         SE     P-value           95% CI"
	disp as text "{hline 63}"
	disp as text "TMLE:    {c |}" %7.1f as result return(ATEtmle) "    " %7.1f as result return(ATE_SE_tmle) "     " %7.4f as result return(ATE_pvalue) as text "     (" %7.1f as result return(ATE_LCIa) ","  %7.1f as result return(ATE_UCIa) as text " )"
	disp as text "{hline 63}"
	disp as text " "
	}

	local rrbin ""CRR: "%4.2f `RRtmle'  "; 95%CI:("%3.2f `LCIrr' ", "%3.2f `UCIrr' ")""
	local orbin ""MOR: "%4.2f `ORtmle'  "; 95%CI:("%3.2f `LCIOr' ", "%3.2f `UCIOr' ")""

	* Display the results of CRR and MOR
	disp as text "{hline 51}"
	di "                           Estimate          95% CI"
	disp as text "{hline 51}"
	disp as text "Causal Risk Ratio:      " "{c |}      "  %4.2f as result `RRtmle' "     (" %3.2f as result `LCIrr' as text ","  %3.2f as result `UCIrr' as text ")"
	disp as text "Marginal Odds Ratio:    " "{c |}      "  %4.2f as result `ORtmle' "     (" %3.2f as result `LCIOr' as text "," %3.2f as result `UCIOr' as text ")"
	disp as text "{hline 51}"

// Drop the variables if the elements option is not specified
	if $keepvars == 0 {
		drop d1 d0
		drop QAW Q1W Q0W
		drop Q1star Qa1star Q0star Qa0star
		drop ATE IC
		drop Y A
		drop POM1 POM0 ps
		*drop cin
		capture drop ytempvar
	}

// Rename and label the variables if the elements option *is* specified
	if $keepvars == 1 {
		* Label the variables
			lab var d1 "Parameter for the influence curve"
			lab var d0 "Parameter for the influence curve"
			lab var QAW "Initial prediction of the outcome"
			lab var Q1W "Initial prediction of the outcome for A = 1"
			lab var Q0W "Initial prediction of the outcome for A = 0"
			lab var Q1star "Update of the initial prediction for A = 1"
			lab var Qa1star "Update of the initial prediction for A = 1"
			lab var Q0star "Update of the initial prediction for A = 0"
			lab var Qa0star "Update of the initial prediction for A = 0"
			lab var A "Exposure/Treatment"
			lab var Y "Outcome"
			lab var ATE "Average Treatment Effect"
			lab var IC "Influence Curve"
			lab var Q1star "Update of initial plug-in estimate for A=1"
			lab var Q0star "Update of initial plug-in estimate for A=0"
			lab var POM1 "Potential Outcome Y(1)"
			lab var POM0 "Potential Otucome Y(0)"
			lab var ps "Propensity Score"
			*lab var cin "Range of Y"
			capture drop ytempvar
		* Rename the elements variables
			foreach var of varlist d1 d0 QAW Q1W Q0W Q1star Qa1star Q0star Qa0star ATE IC  Y A  POM1 POM0 ps  {
				rename `var' _`var'
			}
	}

// Clean up
quietly: rm SLS.R
// quietly: rm SLS.Rout
quietly: rm data2.dta
quietly: rm data.csv
quietly: rm fulldata.csv
// quietly: rm .RData
quietly: memory clean
end




program tmlebgambal, rclass
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
        `"fulldata <- read.csv("fulldata.csv", sep=",")"' _newline ///
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
        `"Q <- try(SuperLearner(Y = data[,1] ,X = X, SL.library=SL.library, family = "binomial", newX=newdata, method ="method.NNLS"), silent=TRUE)"' _newline ///
        `"Q <- as.data.frame(Q[[4]])"' _newline ///
        `"QAW <- Q[1:n,]"' _newline ///
        `"Q1W <- Q[((n+1):(2*n)),]"' _newline ///
        `"Q0W <- Q[((2*n+1):(3*n)),]"' _newline ///
        `"g <- suppressWarnings(SuperLearner(Y = data[,2], X = W, SL.library = SL.library, family = "binomial", method = "method.NNLS"))"' _newline ///
        `"ps <- g[[4]]"' _newline ///
        `"ps[ps<0.025] <- 0.025"' _newline ///
        `"ps[ps>0.975] <- 0.975"' _newline ///
        `"data <- cbind(fulldata,QAW,Q1W,Q0W,ps,Y,A)"' _newline ///
        `"write.dta(data, "data2.dta")"'
qui: file close rcode


	if "`c(os)'" == "MacOSX" {
	qui shell "/usr/local/bin/r" CMD BATCH SLS.R
	}
	else{
	// Write batch file to find R.exe path and R version
		set more off
		qui: file close _all
		qui: file open bat using setup.bat, write replace
		qui: file write bat ///
		`"@echo off"' _newline ///
		`"SET PATHROOT=C:\Program Files\R\"' _newline ///
		`"echo Locating path of R..."' _newline ///
		`"echo."' _newline ///
		`"if not exist "%PATHROOT%" goto:NO_R"' _newline ///
		`":NO_R"' _newline ///
		`"echo R is not installed in your system."' _newline ///
		`"echo."' _newline ///
		`"echo Download it from https://cran.r-project.org/bin/windows/base/"' _newline ///
		`"echo Install it and re-run this script"' _newline ///
		`":DONE"' _newline ///
		`"echo."' _newline ///
		`"pause"'
		qui: file close bat
		
		// 		`"for /f "delims=" %%r in (' dir /b "%PATHROOT%R*" ') do ("' _newline ///
		// 				`"echo Found %%r"' _newline ///
		// 				`"echo ! "%PATHROOT%%%r\bin\x64\R.exe" CMD BATCH SLS.R > runr.do"' _newline ///
		// 				`"echo All set!"' _newline ///
		// 				`"goto:DONE"' _newline ///
		// 		`")"' _newline ///
		
		// 		//Run batch
		// 		! setup.bat
		// 		//Run R
		// 		do runr.do
		
		* Check you have 'rscript' command installed.
		qui: net install rscript, from("https://raw.githubusercontent.com/reifjulian/rscript/master") replace
		
		* Run SLS.R in Windows - 'rscript' command will search for common folders to find R.exe
		qui: rscript using SLS.R
	}

	
	
// Read Revised Data Back to Stata
clear
quietly: use "data2.dta", clear
qui cap drop X__000000
tempvar logQAW logQ1W logQ0W HAW H1W H0W eps1 eps2 eps ATE ICrr ICor

// Create density plot to check positivity assumption
preserve
sort A
qui by A: summarize ps
qui kdensity ps if A==1, generate(x1pointsa d1A) nograph n(10000)
qui kdensity ps if A==0, generate(x0pointsa d0A) nograph n(10000)
label variable d1A "A = 1"
label variable d0A "A = 0"
twoway (line d0A x0pointsa) || ///
		(line d1A x1pointsa), ///
		xtitle("Propensity score") ///
		ytitle("Density") ///
		graphregion(color(white)) bgcolor(white) plotregion(fcolor(white)) ///
		legend(order(2 "Treated" 1 "Not treated") position(11) cols(1) ring(0) ///
				region(style(none)))
restore

// Q to logit scale
gen `logQAW' = log(QAW / (1 - QAW))
gen `logQ1W' = log(Q1W / (1 - Q1W))
gen `logQ0W' = log(Q0W / (1 - Q0W))

// Clever covariate HAW
gen  `HAW' = (A / ps) - ((1 - A) / (1 - ps))
gen  `H1W' = A / ps
gen  `H0W' = (1 - A) / (1 - ps)

// Estimation of the substitution parameter (Epsilon)
qui glm Y `H1W' `H0W', fam(binomial) offset(`logQAW') robust noconstant
mat a= e(b)
gen `eps1' = a[1,1]
gen `eps2' = a[1,2]

qui glm Y `HAW', fam(binomial) offset(`logQAW') robust noconstant
mat a= e(b)
gen `eps' = a[1,1]


// Targeted ATE, update from Q̅^0 (A,W) to Q̅^1 (A,W)
gen double Qa0star = exp(`H0W'*`eps' + `logQ0W')/(1 + exp(`H0W'*`eps' + `logQ0W'))
gen double Qa1star = exp(`H1W'*`eps' + `logQ1W')/(1 + exp(`H1W'*`eps' + `logQ1W'))

gen double Q0star = exp(`H0W'*`eps2' + `logQ0W')/(1 + exp(`H0W'*`eps2' + `logQ0W'))
gen double Q1star = exp(`H1W'*`eps1' + `logQ1W')/(1 + exp(`H1W'*`eps1' + `logQ1W'))

//gen double cin = ($b - $a)

gen double POM1 = cond($flag == 1, Qa1star, (Qa1star * ($b - $a )) + $a , .)
gen double POM0 = cond($flag == 1, Qa0star, (Qa0star * ($b - $a )) + $a , .)

summ POM1 POM0 ps
di as text " "

// Estimating the updated targeted ATE binary outcome
gen double ATE = cond($flag == 1, (Qa1star - Qa0star), ((Qa1star * ($b - $a )) + $a ) - ((Qa0star * ($b - $a )) + $a ) , .)
qui sum ATE
return scalar ATEtmle = r(mean)

// Relative risk and Odss ratio
qui sum Q1star
local Q1 = r(mean)
qui sum Q0star
local Q0 = r(mean)

// Relative risk and Odds ratio
	local RRtmle = `Q1'/`Q0'
	local logRRtmle = log(`Q1') - log(`Q0')
	local ORtmle = (`Q1' * (1 - `Q0')) / ((1 - `Q1') * `Q0')
	local logORtmle = log((`Q1' * (1 - `Q0')) / ((1 - `Q1') * `Q0'))

// Statistical inference (Efficient Influence Curve)
	gen d1 = cond($flag == 1,(A * (Y - Q1star) / ps) + Q1star - `Q1',(A * (Y - Qa1star) / ps) + Qa1star - `Q1' ,.)
	gen d0 = cond($flag == 1,(1 - A) * (Y - Q0star) / (1 - ps) + Q0star - `Q0',(1 - A) * (Y - Qa0star) / (1 - ps) + Qa0star - `Q0' ,.)
	gen IC = cond($flag == 1,(d1 - d0), ((d1 * ($b - $a )) + $a ) - ((d0 * ($b - $a )) + $a ) , .)
	qui sum IC
	return scalar ATE_SE_tmle = sqrt(r(Var)/r(N))

// Statistical inference ATE
	return scalar ATE_pvalue =  2 * (normalden(abs(return(ATEtmle) / (return(ATE_SE_tmle)))))
	return scalar ATE_LCIa   =  return(ATEtmle) - 1.96 * return(ATE_SE_tmle)
	return scalar ATE_UCIa   =  return(ATEtmle) + 1.96 * return(ATE_SE_tmle)

// Statistical inference RR
	gen `ICrr' = (1/`Q1' * d1) + ((1/`Q0') * d0)
	qui sum `ICrr'
	local varICrr = r(Var)/r(N)
	local LCIrr =  exp(`logRRtmle' - 1.96 * sqrt(`varICrr'))
	local UCIrr =  exp(`logRRtmle' + 1.96 * sqrt(`varICrr'))
	
// Statistical inference OR (We do not use the log of the OR
    gen `ICor' = 1/ (`Q1' *(1 - `Q1')) * d1 - 1/ (`Q0' *(1 - `Q0')) * d0 // Partial derivatives (use wolfram to check: d/dx log(x/(1-x))+ log(y/(1-y)))
    qui sum `ICor'
    local varICor = r(Var)/r(N)
	local LCIOr = exp(log(`ORtmle') - 1.96 * sqrt(`varICor'))
	local UCIOr = exp(log(`ORtmle') + 1.96 * sqrt(`varICor'))

// Display Results of ATE
	return scalar CRR = `RRtmle'
	return scalar SE_log_CRR  = sqrt(`varICrr')
	return scalar MOR = `ORtmle'
	return scalar SE_log_MOR  = sqrt(`varICor')

	if $flag==1 {
	disp as text "{hline 63}"
	di "         {c |}" "    ATE         SE     P-value           95% CI"
	disp as text "{hline 63}"
	disp as text "TMLE:    {c |}" %7.4f as result return(ATEtmle) "    " %7.4f as result return(ATE_SE_tmle) "     " %7.4f as result return(ATE_pvalue) as text "     (" %7.4f as result return(ATE_LCIa) ","  %7.4f as result return(ATE_UCIa) as text " )"
	disp as text "{hline 63}"
	disp as text " "
	}
	else if $flag!=1{
	disp as text "{hline 63}"
	di "         {c |}" "    ATE         SE     P-value           95% CI"
	disp as text "{hline 63}"
	disp as text "TMLE:    {c |}" %7.1f as result return(ATEtmle) "    " %7.1f as result return(ATE_SE_tmle) "     " %7.4f as result return(ATE_pvalue) as text "     (" %7.1f as result return(ATE_LCIa) ","  %7.1f as result return(ATE_UCIa) as text " )"
	disp as text "{hline 63}"
	disp as text " "
	}

	local rrbin ""CRR: "%4.2f `RRtmle'  "; 95%CI:("%3.2f `LCIrr' ", "%3.2f `UCIrr' ")""
	local orbin ""MOR: "%4.2f `ORtmle'  "; 95%CI:("%3.2f `LCIOr' ", "%3.2f `UCIOr' ")""

	* Display the results of CRR and MOR
	disp as text "{hline 51}"
	di "                           Estimate          95% CI"
	disp as text "{hline 51}"
	disp as text "Causal Risk Ratio:      " "{c |}      "  %4.2f as result `RRtmle' "     (" %3.2f as result `LCIrr' as text ","  %3.2f as result `UCIrr' as text ")"
	disp as text "Marginal Odds Ratio:    " "{c |}      "  %4.2f as result `ORtmle' "     (" %3.2f as result `LCIOr' as text "," %3.2f as result `UCIOr' as text ")"
	disp as text "{hline 51}"


	* Display covariate balance table
		* Create macros from the varlist
			tokenize $variablelist
			local outcome `1'
			macro shift
			local exposure `1'
			macro shift
			local varlist `*'

		* Create the IPW weights
			capture drop _ipw
			qui gen _ipw = .
			qui replace _ipw = (`exposure'==1) / ps if `exposure'==1
			qui replace _ipw = (`exposure'==0) / (1- ps) if `exposure'==0

		* Layout of the results
			di as text "{hline 67}"
			di as text "                 Standardised Differences            Variance ratio"
			di as text "                          Raw    Weighted           Raw    Weighted"
			di as text "{hline 67}"

			/*
			disp as text "{hline 29}"
			di "TMLE: Causal Risk Ratio (CRR)"
			disp as text "CRR:    " "{c |} "  %4.2f as result `RRtmle'
			disp as text "95% CI: " "{c |} " "(" %3.2f as result `LCIrr' as text ","  %3.2f as result `UCIrr' as text ")"
			disp as text "{hline 29}"


			disp as text "{hline 31}"
			di "TMLE: Marginal Odds Ratio (MOR)"
			disp as text "{hline 31}"
			disp as text "MOR:    " "{c |} "  %4.2f as result `ORtmle'
			disp as text "95% CI: " "{c |} " "(" %3.2f as result `LCIOr' as text "," %3.2f as result `UCIOr' as text ")"
			disp as text "{hline 31}"
			*/


		* Calculate the covariate balance
			foreach var in `varlist' {
					di as text "`var'"

					* Raw SMD
					qui summarize `var' if `exposure'==1
					local m1 = r(mean)
					local v1 = r(Var)
					qui summarize `var' if `exposure'==0
					local m0 = r(mean)
					local v0 = r(Var)
					* Calculate the Standardised "mean difference"
					local rSMD = (`m1' - `m0') / sqrt( (`v1' + `v0') /2 )

					* Weighted SMD
					qui summarize `var' [iw=_ipw] if `exposure'==1
					local m1 = r(mean)
					local v1 = r(Var)
					qui summarize `var' [iw=_ipw] if `exposure'==0
					local m0 = r(mean)
					local v0 = r(Var)
					* Calculate the Standardised "mean difference"
					local wSMD = (`m1' - `m0') / sqrt( (`v1' + `v0') /2 )

					* Raw VR
					qui sum `var' if `exposure' ==1
					local v1 = r(Var)
					qui sum `var' if `exposure' ==0
					local v0 = r(Var)
					* Calculate the variance ratio
					local rVR = `v1' / `v0'

					* Weighted VR
					qui sum `var' [iw=_ipw] if `exposure' ==1
					local v1 = r(Var)
					qui sum `var' [iw=_ipw] if `exposure' ==0
					local v0 = r(Var)
					* Calculate the variance ratio
					local wVR `v1' / `v0'

					* Display the values
					di as text "                    " %9.7g as result `rSMD' as text "   " %9.7g as result `wSMD' as text "     " %9.7g as result `rVR' as text "   " %9.7g as result `wVR'
			}
			di as text "{hline 67}"


// Drop the variables if the elements option is not specified
	if $keepvars == 0 {
		drop d1 d0
		drop QAW Q1W Q0W
		drop Q1star Qa1star Q0star Qa0star
		drop ATE IC
		drop Y A
		drop POM1 POM0 ps
		*drop cin
		capture drop ytempvar
	}

// Rename and label the variables if the elements option *is* specified
	if $keepvars == 1 {
		* Label the variables
			lab var d1 "Parameter for the influence curve"
			lab var d0 "Parameter for the influence curve"
			lab var QAW "Initial prediction of the outcome"
			lab var Q1W "Initial prediction of the outcome for A = 1"
			lab var Q0W "Initial prediction of the outcome for A = 0"
			lab var Q1star "Update of the initial prediction for A = 1"
			lab var Qa1star "Update of the initial prediction for A = 1"
			lab var Q0star "Update of the initial prediction for A = 0"
			lab var Qa0star "Update of the initial prediction for A = 0"
			lab var A "Exposure/Treatment"
			lab var Y "Outcome"
			lab var ATE "Average Treatment Effect"
			lab var IC "Influence Curve"
			lab var Q1star "Update of initial plug-in estimate for A=1"
			lab var Q0star "Update of initial plug-in estimate for A=0"
			lab var POM1 "Potential Outcome Y(1)"
			lab var POM0 "Potential Otucome Y(0)"
			lab var ps "Propensity Score"
			*lab var cin "Range of Y"
			capture drop ytempvar
		* Rename the elements variables
			foreach var of varlist d1 d0 QAW Q1W Q0W Q1star Qa1star Q0star Qa0star ATE IC  Y A  POM1 POM0 ps {
				rename `var' _`var'
			}
	}

// Clean up
quietly: rm SLS.R
// quietly: rm SLS.Rout
quietly: rm data2.dta
quietly: rm data.csv
quietly: rm fulldata.csv
// quietly: rm .RData
quietly: memory clean
end




program tmleglsrf, rclass
// Write R Code dependencies: foreign Surperlearner
set more off
qui: file close _all
qui: file open rcode using SLS.R, write replace
qui: file write rcode ///
		`"set.seed(123)"' _newline ///
        `"list.of.packages <- c("foreign","SuperLearner","glmnet","randomForest")"' _newline ///
        `"new.packages <- list.of.packages[!(list.of.packages %in% installed.packages()[,"Package"])]"' _newline ///
        `"if(length(new.packages)) install.packages(new.packages, repos='http://cran.us.r-project.org')"' _newline ///
        `"library(SuperLearner)"' _newline ///
        `"library(foreign)"' _newline ///
        `"data <- read.csv("data.csv", sep=",")"' _newline ///
        `"fulldata <- read.csv("fulldata.csv", sep=",")"' _newline ///
        `"attach(data)"' _newline ///
        `"SL.library <- c("SL.glm","SL.step","SL.glm.interaction","SL.gam","SL.glmnet","SL.randomForest")"' _newline ///
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
        `"Q <- try(SuperLearner(Y = data[,1] ,X = X, SL.library=SL.library, family = "binomial", newX=newdata, method ="method.NNLS"), silent=TRUE)"' _newline ///
        `"Q <- as.data.frame(Q[[4]])"' _newline ///
        `"QAW <- Q[1:n,]"' _newline ///
        `"Q1W <- Q[((n+1):(2*n)),]"' _newline ///
        `"Q0W <- Q[((2*n+1):(3*n)),]"' _newline ///
        `"g <- suppressWarnings(SuperLearner(Y = data[,2], X = W, SL.library = SL.library, family = "binomial", method = "method.NNLS"))"' _newline ///
        `"ps <- g[[4]]"' _newline ///
        `"ps[ps<0.025] <- 0.025"' _newline ///
        `"ps[ps>0.975] <- 0.975"' _newline ///
        `"data <- cbind(fulldata,QAW,Q1W,Q0W,ps,Y,A)"' _newline ///
        `"write.dta(data, "data2.dta")"'
qui: file close rcode



	if "`c(os)'" == "MacOSX" {
	qui shell "/usr/local/bin/r" CMD BATCH SLS.R
	}
	else{
	// Write batch file to find R.exe path and R version
		set more off
		qui: file close _all
		qui: file open bat using setup.bat, write replace
		qui: file write bat ///
		`"@echo off"' _newline ///
		`"SET PATHROOT=C:\Program Files\R\"' _newline ///
		`"echo Locating path of R..."' _newline ///
		`"echo."' _newline ///
		`"if not exist "%PATHROOT%" goto:NO_R"' _newline ///
		`":NO_R"' _newline ///
		`"echo R is not installed in your system."' _newline ///
		`"echo."' _newline ///
		`"echo Download it from https://cran.r-project.org/bin/windows/base/"' _newline ///
		`"echo Install it and re-run this script"' _newline ///
		`":DONE"' _newline ///
		`"echo."' _newline ///
		`"pause"'
		qui: file close bat
		
		// 		`"for /f "delims=" %%r in (' dir /b "%PATHROOT%R*" ') do ("' _newline ///
		// 				`"echo Found %%r"' _newline ///
		// 				`"echo ! "%PATHROOT%%%r\bin\x64\R.exe" CMD BATCH SLS.R > runr.do"' _newline ///
		// 				`"echo All set!"' _newline ///
		// 				`"goto:DONE"' _newline ///
		// 		`")"' _newline ///
		
		// 		//Run batch
		// 		! setup.bat
		// 		//Run R
		// 		do runr.do
		
		* Check you have 'rscript' command installed.
		qui: net install rscript, from("https://raw.githubusercontent.com/reifjulian/rscript/master") replace
		
		* Run SLS.R in Windows - 'rscript' command will search for common folders to find R.exe
		qui: rscript using SLS.R
	}

	
	
// Read Revised Data Back to Stata
clear
quietly: use "data2.dta", clear
qui cap drop X__000000
tempvar logQAW logQ1W logQ0W HAW H1W H0W eps1 eps2 eps ATE ICrr ICor
// Q to logit scale
gen `logQAW' = log(QAW / (1 - QAW))
gen `logQ1W' = log(Q1W / (1 - Q1W))
gen `logQ0W' = log(Q0W / (1 - Q0W))

// Clever covariate HAW
gen  `HAW' = (A / ps) - ((1 - A) / (1 - ps))
gen  `H1W' = A / ps
gen  `H0W' = (1 - A) / (1 - ps)

// Estimation of the substitution parameter (Epsilon)
qui glm Y `H1W' `H0W', fam(binomial) offset(`logQAW') robust noconstant
mat a= e(b)
gen `eps1' = a[1,1]
gen `eps2' = a[1,2]

qui glm Y `HAW', fam(binomial) offset(`logQAW') robust noconstant
mat a= e(b)
gen `eps' = a[1,1]

// Targeted ATE, update from Q̅^0 (A,W) to Q̅^1 (A,W)
gen double Qa0star = exp(`H0W'*`eps' + `logQ0W')/(1 + exp(`H0W'*`eps' + `logQ0W'))
gen double Qa1star = exp(`H1W'*`eps' + `logQ1W')/(1 + exp(`H1W'*`eps' + `logQ1W'))

gen double Q0star = exp(`H0W'*`eps2' + `logQ0W')/(1 + exp(`H0W'*`eps2' + `logQ0W'))
gen double Q1star = exp(`H1W'*`eps1' + `logQ1W')/(1 + exp(`H1W'*`eps1' + `logQ1W'))

// gen double cin = ($b - $a)

gen double POM1 = cond($flag == 1, Qa1star, (Qa1star * ($b - $a )) + $a , .)
gen double POM0 = cond($flag == 1, Qa0star, (Qa0star * ($b - $a )) + $a , .)

summ POM1 POM0 ps
di as text " "

// Estimating the updated targeted ATE binary outcome
gen double ATE = cond($flag == 1, (Qa1star - Qa0star), ((Qa1star * ($b - $a )) + $a ) - ((Qa0star * ($b - $a )) + $a ) , .)
qui sum ATE
return scalar ATEtmle = r(mean)

// Relative risk and Odss ratio
qui sum Q1star
local Q1 = r(mean)
qui sum Q0star
local Q0 = r(mean)

// Relative risk and Odds ratio
	local RRtmle = `Q1'/`Q0'
	local logRRtmle = log(`Q1') - log(`Q0')
	local ORtmle = (`Q1' * (1 - `Q0')) / ((1 - `Q1') * `Q0')
	local logORtmle = log((`Q1' * (1 - `Q0')) / ((1 - `Q1') * `Q0'))

// Statistical inference (Efficient Influence Curve)
	gen d1 = cond($flag == 1,(A * (Y - Q1star) / ps) + Q1star - `Q1',(A * (Y - Qa1star) / ps) + Qa1star - `Q1' ,.)
	gen d0 = cond($flag == 1,(1 - A) * (Y - Q0star) / (1 - ps) + Q0star - `Q0',(1 - A) * (Y - Qa0star) / (1 - ps) + Qa0star - `Q0' ,.)
	gen IC = cond($flag == 1,(d1 - d0), ((d1 * ($b - $a )) + $a ) - ((d0 * ($b - $a )) + $a ) , .)
	qui sum IC
	return scalar ATE_SE_tmle = sqrt(r(Var)/r(N))

// Statistical inference ATE
	return scalar ATE_pvalue =  2 * (normalden(abs(return(ATEtmle) / (return(ATE_SE_tmle)))))
	return scalar ATE_LCIa   =  return(ATEtmle) - 1.96 * return(ATE_SE_tmle)
	return scalar ATE_UCIa   =  return(ATEtmle) + 1.96 * return(ATE_SE_tmle)

// Statistical inference RR
	gen `ICrr' = (1/`Q1' * d1) + ((1/`Q0') * d0)
	qui sum `ICrr'
	local varICrr = r(Var)/r(N)
	local LCIrr =  exp(`logRRtmle' - 1.96 * sqrt(`varICrr'))
	local UCIrr =  exp(`logRRtmle' + 1.96 * sqrt(`varICrr'))
	
// Statistical inference OR (We do not use the log of the OR
    gen `ICor' = 1/ (`Q1' *(1 - `Q1')) * d1 - 1/ (`Q0' *(1 - `Q0')) * d0 // Partial derivatives (use wolfram to check: d/dx log(x/(1-x))+ log(y/(1-y)))
    qui sum `ICor'
    local varICor = r(Var)/r(N)
	local LCIOr = exp(log(`ORtmle') - 1.96 * sqrt(`varICor'))
	local UCIOr = exp(log(`ORtmle') + 1.96 * sqrt(`varICor'))

// Display Results of ATE
	return scalar CRR = `RRtmle'
	return scalar SE_log_CRR  = sqrt(`varICrr')
	return scalar MOR = `ORtmle'
	return scalar SE_log_MOR  = sqrt(`varICor')

	if $flag==1 {
	disp as text "{hline 63}"
	di "         {c |}" "    ATE         SE     P-value           95% CI"
	disp as text "{hline 63}"
	disp as text "TMLE:    {c |}" %7.4f as result return(ATEtmle) "    " %7.4f as result return(ATE_SE_tmle) "     " %7.4f as result return(ATE_pvalue) as text "     (" %7.4f as result return(ATE_LCIa) ","  %7.4f as result return(ATE_UCIa) as text " )"
	disp as text "{hline 63}"
	disp as text " "
	}
	else if $flag!=1{
	disp as text "{hline 63}"
	di "         {c |}" "    ATE         SE     P-value           95% CI"
	disp as text "{hline 63}"
	disp as text "TMLE:    {c |}" %7.1f as result return(ATEtmle) "    " %7.1f as result return(ATE_SE_tmle) "     " %7.4f as result return(ATE_pvalue) as text "     (" %7.1f as result return(ATE_LCIa) ","  %7.1f as result return(ATE_UCIa) as text " )"
	disp as text "{hline 63}"
	disp as text " "
	}

	local rrbin ""CRR: "%4.2f `RRtmle'  "; 95%CI:("%3.2f `LCIrr' ", "%3.2f `UCIrr' ")""
	local orbin ""MOR: "%4.2f `ORtmle'  "; 95%CI:("%3.2f `LCIOr' ", "%3.2f `UCIOr' ")""

	* Display the results of CRR and MOR
	disp as text "{hline 51}"
	di "                           Estimate          95% CI"
	disp as text "{hline 51}"
	disp as text "Causal Risk Ratio:      " "{c |}      "  %4.2f as result `RRtmle' "     (" %3.2f as result `LCIrr' as text ","  %3.2f as result `UCIrr' as text ")"
	disp as text "Marginal Odds Ratio:    " "{c |}      "  %4.2f as result `ORtmle' "     (" %3.2f as result `LCIOr' as text "," %3.2f as result `UCIOr' as text ")"
	disp as text "{hline 51}"

// Drop the variables if the elements option is not specified
	if $keepvars == 0 {
		drop d1 d0
		drop QAW Q1W Q0W
		drop Q1star Qa1star Q0star Qa0star
		drop ATE IC
		drop Y A
		drop POM1 POM0 ps
		*drop cin
		capture drop ytempvar
	}

// Rename and label the variables if the elements option *is* specified
	if $keepvars == 1 {
		* Label the variables
			lab var d1 "Parameter for the influence curve"
			lab var d0 "Parameter for the influence curve"
			lab var QAW "Initial prediction of the outcome"
			lab var Q1W "Initial prediction of the outcome for A = 1"
			lab var Q0W "Initial prediction of the outcome for A = 0"
			lab var Q1star "Update of the initial prediction for A = 1"
			lab var Qa1star "Update of the initial prediction for A = 1"
			lab var Q0star "Update of the initial prediction for A = 0"
			lab var Qa0star "Update of the initial prediction for A = 0"
			lab var A "Exposure/Treatment"
			lab var Y "Outcome"
			lab var ATE "Average Treatment Effect"
			lab var IC "Influence Curve"
			lab var Q1star "Update of initial plug-in estimate for A=1"
			lab var Q0star "Update of initial plug-in estimate for A=0"
			lab var POM1 "Potential Outcome Y(1)"
			lab var POM0 "Potential Otucome Y(0)"
			lab var ps "Propensity Score"
			*lab var cin "Range of Y"
			capture drop ytempvar
		* Rename the elements variables
			foreach var of varlist d1 d0 QAW Q1W Q0W Q1star Qa1star Q0star Qa0star ATE IC  Y A  POM1 POM0 ps {
				rename `var' _`var'
			}
	}

// Clean up
quietly: rm SLS.R
// quietly: rm SLS.Rout
quietly: rm data2.dta
quietly: rm data.csv
quietly: rm fulldata.csv
// quietly: rm .RData
quietly: memory clean
end




program tmleglsrfbal, rclass
// Write R Code dependencies: foreign Surperlearner
set more off
qui: file close _all
qui: file open rcode using SLS.R, write replace
qui: file write rcode ///
		`"set.seed(123)"' _newline ///
        `"list.of.packages <- c("foreign","SuperLearner","glmnet","randomForest")"' _newline ///
        `"new.packages <- list.of.packages[!(list.of.packages %in% installed.packages()[,"Package"])]"' _newline ///
        `"if(length(new.packages)) install.packages(new.packages, repos='http://cran.us.r-project.org')"' _newline ///
        `"library(SuperLearner)"' _newline ///
        `"library(foreign)"' _newline ///
        `"data <- read.csv("data.csv", sep=",")"' _newline ///
        `"fulldata <- read.csv("fulldata.csv", sep=",")"' _newline ///
        `"attach(data)"' _newline ///
        `"SL.library <- c("SL.glm","SL.step","SL.glm.interaction","SL.gam","SL.glmnet","SL.randomForest")"' _newline ///
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
        `"Q <- try(SuperLearner(Y = data[,1] ,X = X, SL.library=SL.library, family = "binomial", newX=newdata, method ="method.NNLS"), silent=TRUE)"' _newline ///
        `"Q <- as.data.frame(Q[[4]])"' _newline ///
        `"QAW <- Q[1:n,]"' _newline ///
        `"Q1W <- Q[((n+1):(2*n)),]"' _newline ///
        `"Q0W <- Q[((2*n+1):(3*n)),]"' _newline ///
        `"g <- suppressWarnings(SuperLearner(Y = data[,2], X = W, SL.library = SL.library, family = "binomial", method = "method.NNLS"))"' _newline ///
        `"ps <- g[[4]]"' _newline ///
        `"ps[ps<0.025] <- 0.025"' _newline ///
        `"ps[ps>0.975] <- 0.975"' _newline ///
        `"data <- cbind(fulldata,QAW,Q1W,Q0W,ps,Y,A)"' _newline ///
        `"write.dta(data, "data2.dta")"'
qui: file close rcode



	if "`c(os)'" == "MacOSX" {
	qui shell "/usr/local/bin/r" CMD BATCH SLS.R
	}
	else{
	// Write batch file to find R.exe path and R version
		set more off
		qui: file close _all
		qui: file open bat using setup.bat, write replace
		qui: file write bat ///
		`"@echo off"' _newline ///
		`"SET PATHROOT=C:\Program Files\R\"' _newline ///
		`"echo Locating path of R..."' _newline ///
		`"echo."' _newline ///
		`"if not exist "%PATHROOT%" goto:NO_R"' _newline ///
		`":NO_R"' _newline ///
		`"echo R is not installed in your system."' _newline ///
		`"echo."' _newline ///
		`"echo Download it from https://cran.r-project.org/bin/windows/base/"' _newline ///
		`"echo Install it and re-run this script"' _newline ///
		`":DONE"' _newline ///
		`"echo."' _newline ///
		`"pause"'
		qui: file close bat
		
		// 		`"for /f "delims=" %%r in (' dir /b "%PATHROOT%R*" ') do ("' _newline ///
		// 				`"echo Found %%r"' _newline ///
		// 				`"echo ! "%PATHROOT%%%r\bin\x64\R.exe" CMD BATCH SLS.R > runr.do"' _newline ///
		// 				`"echo All set!"' _newline ///
		// 				`"goto:DONE"' _newline ///
		// 		`")"' _newline ///
		
		// 		//Run batch
		// 		! setup.bat
		// 		//Run R
		// 		do runr.do
		
		* Check you have 'rscript' command installed.
		qui: net install rscript, from("https://raw.githubusercontent.com/reifjulian/rscript/master") replace
		
		* Run SLS.R in Windows - 'rscript' command will search for common folders to find R.exe
		qui: rscript using SLS.R
	}

	
	
// Read Revised Data Back to Stata
clear
quietly: use "data2.dta", clear
qui cap drop X__000000
tempvar logQAW logQ1W logQ0W HAW H1W H0W eps1 eps2 eps ATE ICrr ICor

// Create density plot to check positivity assumption
preserve
sort A
qui by A: summarize ps
qui kdensity ps if A==1, generate(x1pointsa d1A) nograph n(10000)
qui kdensity ps if A==0, generate(x0pointsa d0A) nograph n(10000)
label variable d1A "A = 1"
label variable d0A "A = 0"
twoway (line d0A x0pointsa) || ///
		(line d1A x1pointsa), ///
		xtitle("Propensity score") ///
		ytitle("Density") ///
		graphregion(color(white)) bgcolor(white) plotregion(fcolor(white)) ///
		legend(order(2 "Treated" 1 "Not treated") position(11) cols(1) ring(0) ///
				region(style(none)))
restore

// Q to logit scale
gen `logQAW' = log(QAW / (1 - QAW))
gen `logQ1W' = log(Q1W / (1 - Q1W))
gen `logQ0W' = log(Q0W / (1 - Q0W))

// Clever covariate HAW
gen  `HAW' = (A / ps) - ((1 - A) / (1 - ps))
gen  `H1W' = A / ps
gen  `H0W' = (1 - A) / (1 - ps)

// Estimation of the substitution parameter (Epsilon)
qui glm Y `H1W' `H0W', fam(binomial) offset(`logQAW') robust noconstant
mat a= e(b)
gen `eps1' = a[1,1]
gen `eps2' = a[1,2]

qui glm Y `HAW', fam(binomial) offset(`logQAW') robust noconstant
mat a= e(b)
gen `eps' = a[1,1]

// Targeted ATE, update from Q̅^0 (A,W) to Q̅^1 (A,W)
gen double Qa0star = exp(`H0W'*`eps' + `logQ0W')/(1 + exp(`H0W'*`eps' + `logQ0W'))
gen double Qa1star = exp(`H1W'*`eps' + `logQ1W')/(1 + exp(`H1W'*`eps' + `logQ1W'))

gen double Q0star = exp(`H0W'*`eps2' + `logQ0W')/(1 + exp(`H0W'*`eps2' + `logQ0W'))
gen double Q1star = exp(`H1W'*`eps1' + `logQ1W')/(1 + exp(`H1W'*`eps1' + `logQ1W'))

// gen double cin = ($b - $a)

	gen double POM1 = cond($flag == 1, Qa1star, (Qa1star * ($b - $a )) + $a , .)
	gen double POM0 = cond($flag == 1, Qa0star, (Qa0star * ($b - $a )) + $a , .)

summ POM1 POM0 ps
di as text " "

// Estimating the updated targeted ATE binary outcome
	gen double ATE = cond($flag == 1, (Qa1star - Qa0star), ((Qa1star * ($b - $a )) + $a ) - ((Qa0star * ($b - $a )) + $a ) , .)
qui sum ATE
return scalar ATEtmle = r(mean)

// Relative risk and Odds ratio
qui sum Q1star
local Q1 = r(mean)
qui sum Q0star
local Q0 = r(mean)

// Relative risk and Odds ratio
	local RRtmle = `Q1'/`Q0'
	local logRRtmle = log(`Q1') - log(`Q0')
	local ORtmle = (`Q1' * (1 - `Q0')) / ((1 - `Q1') * `Q0')
	local logORtmle = log((`Q1' * (1 - `Q0')) / ((1 - `Q1') * `Q0'))

// Statistical inference (Efficient Influence Curve)
	gen d1 = cond($flag == 1,(A * (Y - Q1star) / ps) + Q1star - `Q1',(A * (Y - Qa1star) / ps) + Qa1star - `Q1' ,.)
	gen d0 = cond($flag == 1,(1 - A) * (Y - Q0star) / (1 - ps) + Q0star - `Q0',(1 - A) * (Y - Qa0star) / (1 - ps) + Qa0star - `Q0' ,.)
	gen IC = cond($flag == 1,(d1 - d0), ((d1 * ($b - $a )) + $a ) - ((d0 * ($b - $a )) + $a ) , .)
	qui sum IC
	return scalar ATE_SE_tmle = sqrt(r(Var)/r(N))

// Statistical inference ATE
	return scalar ATE_pvalue =  2 * (normalden(abs(return(ATEtmle) / (return(ATE_SE_tmle)))))
	return scalar ATE_LCIa   =  return(ATEtmle) - 1.96 * return(ATE_SE_tmle)
	return scalar ATE_UCIa   =  return(ATEtmle) + 1.96 * return(ATE_SE_tmle)

// Statistical inference RR
	gen `ICrr' = (1/`Q1' * d1) + ((1/`Q0') * d0)
	qui sum `ICrr'
	local varICrr = r(Var)/r(N)
	local LCIrr =  exp(`logRRtmle' - 1.96 * sqrt(`varICrr'))
	local UCIrr =  exp(`logRRtmle' + 1.96 * sqrt(`varICrr'))
	
// Statistical inference OR (We do not use the log of the OR
    gen `ICor' = 1/ (`Q1' *(1 - `Q1')) * d1 - 1/ (`Q0' *(1 - `Q0')) * d0 // Partial derivatives (use wolfram to check: d/dx log(x/(1-x))+ log(y/(1-y)))
    qui sum `ICor'
    local varICor = r(Var)/r(N)
	local LCIOr = exp(log(`ORtmle') - 1.96 * sqrt(`varICor'))
	local UCIOr = exp(log(`ORtmle') + 1.96 * sqrt(`varICor'))
	
// Display Results of ATE
	return scalar CRR = `RRtmle'
	return scalar SE_log_CRR  = sqrt(`varICrr')
	return scalar MOR = `ORtmle'
	return scalar SE_log_MOR  = sqrt(`varICor')

	if $flag==1 {
	disp as text "{hline 63}"
	di "         {c |}" "    ATE         SE     P-value           95% CI"
	disp as text "{hline 63}"
	disp as text "TMLE:    {c |}" %7.4f as result return(ATEtmle) "    " %7.4f as result return(ATE_SE_tmle) "     " %7.4f as result return(ATE_pvalue) as text "     (" %7.4f as result return(ATE_LCIa) ","  %7.4f as result return(ATE_UCIa) as text " )"
	disp as text "{hline 63}"
	disp as text " "
	}
	else if $flag!=1{
	disp as text "{hline 63}"
	di "         {c |}" "    ATE         SE     P-value           95% CI"
	disp as text "{hline 63}"
	disp as text "TMLE:    {c |}" %7.1f as result return(ATEtmle) "    " %7.1f as result return(ATE_SE_tmle) "     " %7.4f as result return(ATE_pvalue) as text "     (" %7.1f as result return(ATE_LCIa) ","  %7.1f as result return(ATE_UCIa) as text " )"
	disp as text "{hline 63}"
	disp as text " "
	}

	local rrbin ""CRR: "%4.2f `RRtmle'  "; 95%CI:("%3.2f `LCIrr' ", "%3.2f `UCIrr' ")""
	local orbin ""MOR: "%4.2f `ORtmle'  "; 95%CI:("%3.2f `LCIOr' ", "%3.2f `UCIOr' ")""

	* Display the results of CRR and MOR
	disp as text "{hline 51}"
	di "                           Estimate          95% CI"
	disp as text "{hline 51}"
	disp as text "Causal Risk Ratio:      " "{c |}      "  %4.2f as result `RRtmle' "     (" %3.2f as result `LCIrr' as text ","  %3.2f as result `UCIrr' as text ")"
	disp as text "Marginal Odds Ratio:    " "{c |}      "  %4.2f as result `ORtmle' "     (" %3.2f as result `LCIOr' as text "," %3.2f as result `UCIOr' as text ")"
	disp as text "{hline 51}"

// Display covariate balance table
		* Create macros from the varlist
			tokenize $variablelist
			local outcome `1'
			macro shift
			local exposure `1'
			macro shift
			local varlist `*'

		* Create the IPW weights
			capture drop _ipw
			qui gen _ipw = .
			qui replace _ipw = (`exposure'==1) / ps if `exposure'==1
			qui replace _ipw = (`exposure'==0) / (1- ps) if `exposure'==0

		* Layout of the results
			di as text "{hline 67}"
			di as text "                 Standardised Differences            Variance ratio"
			di as text "                          Raw    Weighted           Raw    Weighted"
			di as text "{hline 67}"

		* Calculate the covariate balance
			foreach var in `varlist' {
					di as text "`var'"

					* Raw SMD
					qui summarize `var' if `exposure'==1
					local m1 = r(mean)
					local v1 = r(Var)
					qui summarize `var' if `exposure'==0
					local m0 = r(mean)
					local v0 = r(Var)
					* Calculate the Standardised "mean difference"
					local rSMD = (`m1' - `m0') / sqrt( (`v1' + `v0') /2 )

					* Weighted SMD
					qui summarize `var' [iw=_ipw] if `exposure'==1
					local m1 = r(mean)
					local v1 = r(Var)
					qui summarize `var' [iw=_ipw] if `exposure'==0
					local m0 = r(mean)
					local v0 = r(Var)
					* Calculate the Standardised "mean difference"
					local wSMD = (`m1' - `m0') / sqrt( (`v1' + `v0') /2 )

					* Raw VR
					qui sum `var' if `exposure' ==1
					local v1 = r(Var)
					qui sum `var' if `exposure' ==0
					local v0 = r(Var)
					* Calculate the variance ratio
					local rVR = `v1' / `v0'

					* Weighted VR
					qui sum `var' [iw=_ipw] if `exposure' ==1
					local v1 = r(Var)
					qui sum `var' [iw=_ipw] if `exposure' ==0
					local v0 = r(Var)
					* Calculate the variance ratio
					local wVR `v1' / `v0'

					* Display the values
					di as text "                    " %9.7g as result `rSMD' as text "   " %9.7g as result `wSMD' as text "     " %9.7g as result `rVR' as text "   " %9.7g as result `wVR'
			}
			di as text "{hline 67}"

// Drop the variables if the elements option is not specified
	if $keepvars == 0 {
		drop d1 d0
		drop QAW Q1W Q0W
		drop Q1star Qa1star Q0star Qa0star
		drop ATE IC
		drop Y A
		drop POM1 POM0 ps
		*drop cin
		capture drop ytempvar
	}


// Rename and label the variables if the elements option *is* specified
	if $keepvars == 1 {
		* Label the variables
			lab var d1 "Parameter for the influence curve"
			lab var d0 "Parameter for the influence curve"
			lab var QAW "Initial prediction of the outcome"
			lab var Q1W "Initial prediction of the outcome for A = 1"
			lab var Q0W "Initial prediction of the outcome for A = 0"
			lab var Q1star "Update of the initial prediction for A = 1"
			lab var Qa1star "Update of the initial prediction for A = 1"
			lab var Q0star "Update of the initial prediction for A = 0"
			lab var Qa0star "Update of the initial prediction for A = 0"
			lab var A "Exposure/Treatment"
			lab var Y "Outcome"
			lab var ATE "Average Treatment Effect"
			lab var IC "Influence Curve"
			lab var Q1star "Update of initial plug-in estimate for A=1"
			lab var Q0star "Update of initial plug-in estimate for A=0"
			lab var POM1 "Potential Outcome Y(1)"
			lab var POM0 "Potential Otucome Y(0)"
			lab var ps "Propensity Score"
			*lab var cin "Range of Y"
			capture drop ytempvar
		* Rename the elements variables
			foreach var of varlist d1 d0 QAW Q1W Q0W Q1star Qa1star Q0star Qa0star ATE IC  Y A  POM1 POM0 ps {
				rename `var' _`var'
			}
	}

// Clean up
quietly: rm SLS.R
// quietly: rm SLS.Rout
quietly: rm data2.dta
quietly: rm data.csv
quietly: rm fulldata.csv
// quietly: rm .RData
quietly: memory clean
end

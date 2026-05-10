*! version 4.0.3  09.May.2026
*! ELTMLE: Stata module for Ensemble Learning Targeted Maximum Likelihood Estimation
*! by Miguel Angel Luque-Fernandez [cre,aut]
*! and Matthew J. Smith [aut]
*! and Camille Maringe [aut]
*! Bug reports:
*! miguel-angel.luque at lshtm.ac.uk
*! matt.smith at lshtm.ac.uk
*! camille.maringe at lshtm.ac.


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
** December 2024
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
* The **seed** option is no longer necessary but remains in the syntax for reproducibility.
* Added a progress bar for cross-validation.

* February 2025
* Updated estimation of the substitution parameter (Epsilon) as MLE weight

* June 2025
* Updated cross-validated procedure for the TMLE estimation

* April 2026
* Updated density plot to check positivity assumption

* April 2026
* Refactored: extracted 10 shared helper programs (_eltmle_*) to eliminate
* copy-paste duplication across estimator variants. Behaviour unchanged.

* May 2026 (v4.0.1)
/*
  Counterfactual Q1star/Q0star: previously H1W = A/ps and H0W = (1-A)/(1-ps)
  were used in the targeting update, so Q1star for A=0 obs equalled Q1W
  (untargeted initial estimate) and Q0star for A=1 obs equalled Q0W.
  Fixed to use counterfactual clever covariates 1/ps and 1/(1-ps) for ALL
  observations, so every unit receives a properly targeted counterfactual Q.
  This aligns POM1, POM0, ATE, and the IC with the canonical TMLE EIF.
  Applies to both two-epsilon (Q1star/Q0star) and single-epsilon (Qa1star/Qa0star) 
  branches.

  Display format (continuous, May 2026): %7.4f caused scientific notation
  for values >= 100 (e.g., bweight ATE = -230g displayed as -2.3e+02).
  Changed to %8.1f for ATE/SE/CI in continuous branch; %6.4f for p-value.
  Results now match help file: -230.0, 24.5, (-278.0, -182.0).

  P-value: normalden() returns the PDF, not a tail probability.
  Fixed to 2*normal(-abs(z)), the correct two-sided Wald p-value.
*/

* May 2026 (v4.0.2)
/*
  FOUR bugs fixed: a binary/continuous flag mis-detection, a data-plumbing
  collision (most impactful for users whose outcome/treatment variables are
  literally named "Y"/"A"), and two targeting specification deviations from
  R tmle::tmle.

  Bug -1 (continuous outcome in [0,1] silently mis-classified as binary):
    The non-CV main entry used `global flag = cond(`yvar'<=1, 1, 0)`, which
    evaluates Y[1] -- so any continuous outcome with values bounded in [0,1]
    (e.g. simulated E[Y|A,W] probabilities, ratios, proportions, or anything
    from a [0,1]-bounded scale) was wrongly flagged binary.  Consequence:
    the back-transform `(b - a)` was skipped, and the displayed ATE was on
    the bounded [0,1] scale rather than the original outcome scale,
    producing apparent bias when compared to teffects (which always reports
    on the original scale).  Fix: use `tabulate; flag = cond(r(r)==2, 1, 0)`,
    matching the CV branch.

  Bug 0 (data plumbing -- only manifests when outcome named "Y" or treatment
         named "A"):
    The R-side script wrote data2.dta as cbind(fulldata, ..., Y, A), where
    fulldata.csv carried the user's ORIGINAL (unscaled) outcome and
    treatment, and the explicit Y, A came from data.csv (SCALED).  When the
    user's outcome is literally named "Y" (or treatment "A"), this produced
    a .dta file with two columns named "Y" -- foreign::write.dta does not
    rename duplicates, and Stata reads both with the same name.  Subsequent
    references to `Y' in _eltmle_tmle_estimate resolved to the FIRST (i.e.
    unscaled) column, while QAW/Q1W/Q0W were predictions on the SCALED
    Y* in [0,1].  The targeting GLM therefore tried to fit unscaled Y to
    a scaled offset, producing biased epsilon and biased Q*(a,W).  Effect
    on sim.do (Y in [0.10, 0.99]): ATE bias of ~+0.04 with 0% CI coverage.
    Fix: drop the user's outcome/treatment columns from fulldata before
    cbind in both _eltmle_write_r_noncv and _eltmle_write_r_cv.

  Bug 1 (catastrophic, continuous outcome path):
    The previous version fit a single-epsilon submodel
        glm Y [pweight=HAW], offset(logQAW)            (intercept-only)
    whose implied fluctuation is a CONSTANT epsilon
        Q*(A,W) = expit(logit(Q(A,W)) + epsilon)
    but then formed the counterfactual updates as
        Qa1star = expit(logit(Q1W) + (1/ps)    * epsilon)
        Qa0star = expit(logit(Q0W) + (1/(1-ps))* epsilon)
    -- a *varying* fluctuation that does NOT correspond to the submodel that
    was fit.  Score equation of the EIF was therefore not solved, the
    counterfactual Q's were biased, and the IC variance estimator severely
    under-estimated the true sampling variance.  

  Bug 2 (mild, binary outcome path):
    The previous version fit
        glm Y H1W H0W [pweight=HAW], fam(binomial) offset(logQAW), nocons
    which combines IPW *both* in the weight AND in the clever covariate.
    The score actually solved is sum (Y-mu)/ps^2 = 0 (not the canonical
    sum (Y-mu)/ps = 0 of the EIF).  In MC bias was small for this DGP, but
    the form is non-canonical and not what R tmle does.

  Fix:
    Use a single canonical two-epsilon targeting GLM for BOTH binary and
    continuous outcomes (continuous Y was already mapped to Y* in [0,1]
    before the SuperLearner fit, so the binomial fluctuation is on the
    bounded scale).  No [pweight=HAW] -- weights are 1.

        glm Y H1W H0W, fam(binomial) offset(logQAW) noconstant
        Q1star = invlogit(logQ1W + eps1 * (1/ps))
        Q0star = invlogit(logQ0W + eps2 * (1/(1-ps)))

    Q1star/Q0star are now used uniformly for ATE, POM, IC, CRR, MOR; the
    redundant Qa1star/Qa0star path has been removed.  IC formula matches
    R tmle's calcParameters exactly.  Sign of the IC for log-RR was also
    corrected (was `+`, should be `-` per R tmle line 1623-1624).
*/

* May 2026 (v4.0.3)
/*
  Bug: non-CV estimators (tmle, tmlebal, tmlebgam, tmlebgambal, tmleglsrf,
  tmleglsrfbal) clobbered the user's in-memory data.  The subprograms did
  `clear; use "data2.dta"`, then `_eltmle_label_rename_noncv` dropped Y, A,
  and the prediction columns -- leaving the user with a partially gutted
  copy of `data2.dta` (covariates from fulldata, but no Y, no A) instead of
  the original dataset.

  Cleanup fix (same release): the CV variants (cvtmle, cvtmlebgam,
  cvtmleglsrf) left the user's dataset polluted when called with the
  `elements` option -- foldid, rowid, and the merge indicator _merge
  were never dropped, and the scratch elementsdata.dta file was never
  deleted.  Consolidated post-restore cleanup into a single helper,
  _eltmle_merge_elements_cv, that performs the merge with `nogenerate`,
  removes elementsdata.dta from disk, and drops foldid/rowid in both
  the keepvars=0 and keepvars=1 paths.  Also added `ytempvar` to the
  capture-drop list at the start of eltmle main so a leftover from a
  prior run cannot interfere with the current call.

  Continuous-CI display fix (same release): the upper-bound 95% CI for
  continuous outcomes was displayed as "...0 )" -- a stray space before
  the closing bracket made it look detached from the number.  Trimmed
  the trailing " )" to ")" so the bracket abuts the value cleanly.
*/



capture program drop eltmle
program define eltmle
	syntax varlist(min=3) [if] [pw] [, tmle tmlebgam tmleglsrf bal elements cvtmle cvtmlebgam cvtmleglsrf cvfolds(int 10) seed(numlist)]
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
				capture drop ytempvar


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

		//marksample touse
		global variablelist `varlist'
		local var  "`varlist'" // if `touse'

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
		// Detect binary outcome by counting distinct values, NOT by checking
		// `yvar'<=1.  The latter is buggy: when the user's outcome is a
		// continuous probability in [0,1] (e.g. simulated E[Y|A,W]) every Y
		// value is <=1 and the flag is wrongly set to 1, causing the back-
		// transform `(b-a)` to be skipped.  The CV branch already uses
		// tabulate; align the non-CV branch with it. (v4.0.2)
		qui tabulate `yvar'
		global flag = cond(r(r)==2, 1, 0)
		qui sum `yvar'
		global b = `r(max)'
		global a = `r(min)'
		// Save the user's original outcome so we can restore it after exporting
		// the rescaled copy to data.csv.  Without this, repeated calls to
		// `eltmle ... tmle' on the same dataset (or a tmle call followed by a
		// cvtmle call, as in sim.do) silently re-rescale an already-rescaled
		// Y, putting subsequent estimates on a different scale than the truth.
		tempvar __eltmle_orig_y
		qui gen double `__eltmle_orig_y' = `yvar'
		qui replace `yvar' = (`yvar' - `r(min)') / (`r(max)' - `r(min)') // if `yvar'>1
		local dir `c(pwd)'
		qui cd "`dir'"
		tempfile data
		qui save "`data'.dta", replace
		qui export delimited `var' using "data.csv", nolabel replace
		qui replace `yvar' = `__eltmle_orig_y'
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

		if "`cvtmle'" == "cvtmle" & "`cvfolds'" == "" {
			cvtmle `varlist'
		}
		else if "`cvtmle'" == "cvtmle" & "`cvfolds'" != "" {
			cvtmle `varlist'
		}
		else if "`cvtmleglsrf'" == "cvtmleglsrf" & "`cvfolds'" == "" {
			cvtmleglsrf `varlist'
		}
		else if "`cvtmleglsrf'" == "cvtmleglsrf" & "`cvfolds'" != "" {
			cvtmleglsrf `varlist'
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


// ====================================================
// INTERNAL HELPER PROGRAMS
// Not intended to be called directly by users.
// ====================================================

// ---------------------------------------------------------------------------
// _eltmle_run_r: OS-specific R execution
// Assumes SLS.R has already been written to disk.
// ---------------------------------------------------------------------------
program _eltmle_run_r
	if "`c(os)'" == "MacOSX" {
		local wdir = c(pwd)
		// Try common R locations on macOS
		local r_found = 0
		foreach rpath in "/usr/local/bin/Rscript" "/opt/homebrew/bin/Rscript" "/usr/bin/Rscript" "/Library/Frameworks/R.framework/Resources/bin/Rscript" {
			quietly: capture confirm file "`rpath'"
			if _rc == 0 & `r_found' == 0 {
				local r_found = 1
				local rscript_path "`rpath'"
			}
		}
		if `r_found' == 0 {
			di as error "Cannot find Rscript on this system."
			di as error "Tried: /usr/local/bin/Rscript, /opt/homebrew/bin/Rscript, /usr/bin/Rscript"
			di as error "/Library/Frameworks/R.framework/Resources/bin/Rscript"
			di as error "Please install R from https://cran.r-project.org/"
			exit 601
		}
		qui shell "`rscript_path'" SLS.R > SLS.Rout 2>&1
		quietly: capture confirm file "data2.dta"
		if _rc != 0 {
			di as error "R script failed to produce data2.dta."
			di as error "Working directory: `wdir'"
			di as error "--- R error log (SLS.Rout) ---"
			type SLS.Rout
			di as error "--- end of R log ---"
			quietly: capture rm SLS.Rout
			exit 601
		}
		quietly: capture rm SLS.Rout
	}
	else {
		// Write batch file (informational only; rscript handles actual execution)
		set more off
		qui: file close _all
		qui: file open bat using setup.bat, write replace
		qui: file write bat `"@echo off"' _newline
		qui: file write bat `"SET PATHROOT=C:\Program Files\R\"' _newline
		qui: file write bat `"echo Locating path of R..."' _newline
		qui: file write bat `"echo."' _newline
		qui: file write bat `"if not exist "%PATHROOT%" goto:NO_R"' _newline
		qui: file write bat `":NO_R"' _newline
		qui: file write bat `"echo R is not installed in your system."' _newline
		qui: file write bat `"echo."' _newline
		qui: file write bat `"echo Download it from https://cran.r-project.org/bin/windows/base/"' _newline
		qui: file write bat `"echo Install it and re-run this script"' _newline
		qui: file write bat `":DONE"' _newline
		qui: file write bat `"echo."' _newline
		qui: file write bat `"pause"'
		qui: file close bat
		// Install rscript if needed and run SLS.R
		qui: net install rscript, from("https://raw.githubusercontent.com/reifjulian/rscript/master") replace
		qui: rscript using SLS.R
		quietly: capture rm SLS.Rout
		quietly: capture confirm file "data2.dta"
		if _rc != 0 {
			di as error "R script failed to produce data2.dta."
			di as error "Check that R and required packages (SuperLearner, foreign, gam, arm) are installed."
			exit 601
		}
	}
end


// ---------------------------------------------------------------------------
// _eltmle_write_r_noncv: write SLS.R for non-CV estimators
// Arg: learner = base | bgam | glsrf
// ---------------------------------------------------------------------------
program _eltmle_write_r_noncv
    args learner
    set more off
    qui: file close _all
    qui: file open rcode using SLS.R, write replace
    qui: file write rcode `"set.seed(123)"' _newline
    local wdir = c(pwd)
    local wdir = subinstr("`wdir'","\","/",.)
    qui: file write rcode `"setwd("`wdir'")"' _newline
    // Learner-specific package list
    if "`learner'" == "base" {
        qui: file write rcode `"list.of.packages <- c("foreign","SuperLearner")"' _newline
    }
    else if "`learner'" == "bgam" {
        qui: file write rcode `"list.of.packages <- c("foreign","SuperLearner","gam","arm")"' _newline
    }
    else if "`learner'" == "glsrf" {
        qui: file write rcode `"list.of.packages <- c("foreign","SuperLearner","glmnet","randomForest")"' _newline
    }
    qui: file write rcode `"new.packages <- list.of.packages[!(list.of.packages %in% installed.packages()[,"Package"])]"' _newline
    qui: file write rcode `"userlib <- Sys.getenv('R_LIBS_USER')"' _newline
    qui: file write rcode `"if (!dir.exists(userlib)) dir.create(userlib, recursive = TRUE)"' _newline
    qui: file write rcode `".libPaths(c(userlib, .libPaths()))"' _newline
    qui: file write rcode `"if(length(new.packages)) install.packages(new.packages, repos='https://cran.r-project.org', lib=userlib)"' _newline
    qui: file write rcode `"library(SuperLearner)"' _newline
    qui: file write rcode `"library(foreign)"' _newline
    // Learner-specific extra libraries
    if "`learner'" == "bgam" {
        qui: file write rcode `"library(gam)"' _newline
        qui: file write rcode `"library(arm)"' _newline
    }
    else if "`learner'" == "glsrf" {
        qui: file write rcode `"library(glmnet)"' _newline
        qui: file write rcode `"library(randomForest)"' _newline
    }
    qui: file write rcode `"data <- read.csv("data.csv", sep=",")"' _newline
    qui: file write rcode `"fulldata <- read.csv("fulldata.csv", sep=",")"' _newline
    qui: file write rcode `"attach(data)"' _newline
    // Learner-specific SL library
    if "`learner'" == "base" {
        qui: file write rcode `"SL.library <- c("SL.glm","SL.step","SL.glm.interaction")"' _newline
    }
    else if "`learner'" == "bgam" {
        qui: file write rcode `"SL.library <- c("SL.glm","SL.step","SL.glm.interaction","SL.gam","SL.bayesglm")"' _newline
    }
    else if "`learner'" == "glsrf" {
        qui: file write rcode `"SL.library <- c("SL.glm","SL.step","SL.glm.interaction","SL.gam","SL.glmnet","SL.randomForest")"' _newline
    }
    // Common R code
    qui: file write rcode `"n <- nrow(data)"' _newline
    qui: file write rcode `"nvar <- dim(data)[[2]]"' _newline
    qui: file write rcode `"Y <- data[,1]"' _newline
    qui: file write rcode `"A <- data[,2]"' _newline
    qui: file write rcode `"X <- data[,2:nvar]"' _newline
    qui: file write rcode `"W <- as.data.frame(data[,3:nvar])"' _newline
    qui: file write rcode `"X1 <- X0 <- X"' _newline
    qui: file write rcode `"X1[,1] <- 1"' _newline
    qui: file write rcode `"X0[,1] <- 0"' _newline
    qui: file write rcode `"newdata <- rbind(X,X1,X0)"' _newline
    // SMART STRATIFICATION 
    qui: file write rcode `"is_bin_Y <- length(unique(data[,1])) == 2"' _newline
    qui: file write rcode `"Q <- SuperLearner(Y = as.numeric(data[,1]), X = X, SL.library=SL.library, family = "binomial", newX=newdata, method ="method.NNLS", cvControl = list(V = 10, stratifyCV = is_bin_Y))"' _newline
    qui: file write rcode `"Q <- as.data.frame(Q[[4]])"' _newline
    qui: file write rcode `"QAW <- Q[1:n,]"' _newline
    qui: file write rcode `"Q1W <- Q[((n+1):(2*n)),]"' _newline
    qui: file write rcode `"Q0W <- Q[((2*n+1):(3*n)),]"' _newline
    qui: file write rcode `"g <- SuperLearner(Y = as.numeric(data[,2]), X = W, SL.library = SL.library, family = "binomial", method = "method.NNLS", cvControl = list(V = 10, stratifyCV = TRUE))"' _newline
    qui: file write rcode `"ps <- g[[4]]"' _newline
    qui: file write rcode `"ps[ps<0.025] <- 0.025"' _newline
    qui: file write rcode `"ps[ps>0.975] <- 0.975"' _newline
    qui: file write rcode `"yname <- colnames(data)[1]; aname <- colnames(data)[2]"' _newline
    qui: file write rcode `"fulldata <- fulldata[, !(colnames(fulldata) %in% c(yname, aname)), drop = FALSE]"' _newline
    qui: file write rcode `"data <- cbind(fulldata,QAW,Q1W,Q0W,ps,Y,A)"' _newline
    qui: file write rcode `"write.dta(data, "data2.dta")"' _newline
    qui: file close rcode
end


// ---------------------------------------------------------------------------
// _eltmle_write_r_cv: write SLS.R for CV estimators (called once per fold)
// Arg: learner = base | bgam | glsrf
// ---------------------------------------------------------------------------
program _eltmle_write_r_cv
    args learner
    qui: set more off
    qui: file close _all
    qui: file open rcode using SLS.R, write replace
    qui: file write rcode `"set.seed(123)"' _newline
    local wdir = c(pwd)
    local wdir = subinstr("`wdir'","\","/",.)
    qui: file write rcode `"setwd("`wdir'")"' _newline
    // Learner-specific package list (all CV variants include dplyr)
    if "`learner'" == "base" {
        qui: file write rcode `"list.of.packages <- c("foreign","SuperLearner","dplyr")"' _newline
    }
    else if "`learner'" == "bgam" {
        qui: file write rcode `"list.of.packages <- c("foreign","SuperLearner","dplyr","gam","arm")"' _newline
    }
    else if "`learner'" == "glsrf" {
        qui: file write rcode `"list.of.packages <- c("foreign","SuperLearner","dplyr","glmnet","randomForest")"' _newline
    }
    qui: file write rcode `"new.packages <- list.of.packages[!(list.of.packages %in% installed.packages()[,"Package"])]"' _newline
    qui: file write rcode `"userlib <- Sys.getenv('R_LIBS_USER')"' _newline
    qui: file write rcode `"if (!dir.exists(userlib)) dir.create(userlib, recursive = TRUE)"' _newline
    qui: file write rcode `".libPaths(c(userlib, .libPaths()))"' _newline
    qui: file write rcode `"if(length(new.packages)) install.packages(new.packages, repos='https://cran.r-project.org', lib=userlib)"' _newline
    qui: file write rcode `"library(SuperLearner)"' _newline
    qui: file write rcode `"library(foreign)"' _newline
    qui: file write rcode `"library(dplyr)"' _newline
    // Learner-specific extra libraries
    if "`learner'" == "bgam" {
        qui: file write rcode `"library(gam)"' _newline
        qui: file write rcode `"library(arm)"' _newline
    }
    else if "`learner'" == "glsrf" {
        qui: file write rcode `"library(glmnet)"' _newline
        qui: file write rcode `"library(randomForest)"' _newline
    }
    qui: file write rcode `"data <- read.csv("cvdata.csv", sep=",")"' _newline
    qui: file write rcode `"attach(data)"' _newline
    // Learner-specific SL library
    if "`learner'" == "base" {
        qui: file write rcode `"SL.library <- c("SL.glm","SL.step","SL.glm.interaction")"' _newline
    }
    else if "`learner'" == "bgam" {
        qui: file write rcode `"SL.library <- c("SL.glm","SL.step","SL.glm.interaction","SL.gam","SL.bayesglm")"' _newline
    }
    else if "`learner'" == "glsrf" {
        qui: file write rcode `"SL.library <- c("SL.glm","SL.step","SL.glm.interaction","SL.gam","SL.glmnet","SL.randomForest")"' _newline
    }
    // Common CV R code (cross-validated fold logic)
    qui: file write rcode `"nvar <- dim(data)[[2]]"' _newline
    qui: file write rcode `"f <- foldnumber[1]"' _newline
    qui: file write rcode `"data <- data %>% dplyr::select(-c("foldnumber"))"' _newline
    qui: file write rcode `"tdata <- data[foldid!=f,]"' _newline
    qui: file write rcode `"vdataAa <- data[foldid==f,]"' _newline
    qui: file write rcode `"vdataA1 <- vdataAa; vdataA1[,2] <- 1"' _newline
    qui: file write rcode `"vdataA0 <- vdataAa; vdataA0[,2] <- 0"' _newline
    qui: file write rcode `"Vn <- nrow(vdataAa)"' _newline
    qui: file write rcode `"Tn <- nrow(tdata)"' _newline
    qui: file write rcode `"Y <- vdataAa[,1]"' _newline
    qui: file write rcode `"A <- vdataAa[,2]"' _newline
    qui: file write rcode `"X <- tdata %>% dplyr::select(-c(1, "foldid", "rowid"))"' _newline
    qui: file write rcode `"W <- tdata %>% dplyr::select(-c(1, 2, "foldid", "rowid"))"' _newline
    
    // SMART STRATIFICATION 
    qui: file write rcode `"is_bin_Y <- length(unique(tdata[,1])) == 2"' _newline
    qui: file write rcode `"tset <- SuperLearner(Y = as.numeric(tdata[,1]), X = X, SL.library = SL.library, family = "binomial", method = "method.NNLS", cvControl = list(V = 10, stratifyCV = is_bin_Y))"' _newline
    qui: file write rcode `"g <- SuperLearner(Y = as.numeric(tdata[,2]), X = W, SL.library = SL.library, family = "binomial", method = "method.NNLS", cvControl = list(V = 10, stratifyCV = TRUE))"' _newline
    
    qui: file write rcode `"QAW <- predict(tset, newdata=vdataAa)[[1]]"' _newline
    qui: file write rcode `"Q1W <- predict(tset, newdata=vdataA1)[[1]]"' _newline
    qui: file write rcode `"Q0W <- predict(tset, newdata=vdataA0)[[1]]"' _newline
    qui: file write rcode `"QAWt <- predict(tset)[[1]]"' _newline
    qui: file write rcode `"ps <- predict(g, newdata=vdataAa)[[1]]"' _newline
    qui: file write rcode `"ps[ps<0.025] <- 0.025"' _newline
    qui: file write rcode `"ps[ps>0.975] <- 0.975"' _newline
    
    // Drop user's outcome/treatment columns from vdataAa to avoid duplicate
    // "Y" and "A" columns in data2.dta
    qui: file write rcode `"yname <- colnames(vdataAa)[1]; aname <- colnames(vdataAa)[2]"' _newline
    qui: file write rcode `"vdataAa <- vdataAa[, !(colnames(vdataAa) %in% c(yname, aname)), drop = FALSE]"' _newline
    qui: file write rcode `"data2 <- cbind(vdataAa, QAW, Q1W, Q0W, ps, Y, A)"' _newline
    qui: file write rcode `"write.dta(data2, "data2.dta")"' _newline
    qui: file close rcode
end

// ---------------------------------------------------------------------------
// _eltmle_tmle_estimate: shared TMLE computation and display
// Assumes data2.dta is already loaded. Uses globals: $flag, $a, $b.
// Creates permanent dataset vars: d1, d0, ATE, IC, Q1star, Q0star,
//   Qa1star, Qa0star, POM1, POM0.
// Returns: ATEtmle, ATE_SE_tmle, ATE_pvalue, ATE_LCIa, ATE_UCIa,
//          CRR, SE_log_CRR, MOR, SE_log_MOR.
// Caller must call "return add" immediately after to propagate scalars.
// ---------------------------------------------------------------------------
program _eltmle_tmle_estimate, rclass
	tempvar logQAW logQ1W logQ0W H1W H0W ICrr ICor

	// Q to logit scale
	gen `logQAW' = log(QAW / (1 - QAW))
	gen `logQ1W' = log(Q1W / (1 - Q1W))
	gen `logQ0W' = log(Q0W / (1 - Q0W))

	// Canonical EIF clever covariates
	gen `H1W' = A / ps
	gen `H0W' = (1 - A) / (1 - ps)

	// ----------------------------------------------------------------------
	// Canonical two-epsilon TMLE targeting step (matches R tmle::tmle).
	//   Submodel:  logit(Q*(A,W)) = logit(Q(A,W)) + eps1*H1W + eps2*H0W
	//   Score eqs: sum (A/ps)(Y - Q*) = 0  AND  sum ((1-A)/(1-ps))(Y - Q*) = 0
	// Same form for binary AND continuous outcomes (continuous Y was already
	// mapped to Y* in [0,1] before the SuperLearner fit, so this binomial
	// fluctuation is on the bounded scale).  No [pweight=HAW] -- weights = 1.
	// See changelog v4.0.2 for the previous, broken specification.
	// ----------------------------------------------------------------------
	qui glm Y `H1W' `H0W', fam(binomial) offset(`logQAW') noconstant
	mat a = e(b)
	local eps1 = a[1,1]   // coefficient on H1W (treated)
	local eps2 = a[1,2]   // coefficient on H0W (controls)

	// Counterfactual updates: 1/ps and 1/(1-ps) applied to ALL observations
	gen double Q1star = invlogit(`logQ1W' + (1/ps)     * `eps1')
	gen double Q0star = invlogit(`logQ0W' + (1/(1-ps)) * `eps2')

	qui sum Q1star
	local Q1 = r(mean)
	qui sum Q0star
	local Q0 = r(mean)

	// Potential outcomes (back-transform to original scale for continuous Y)
	gen double POM1 = cond($flag == 1, Q1star, (Q1star * ($b - $a)) + $a, .)
	gen double POM0 = cond($flag == 1, Q0star, (Q0star * ($b - $a)) + $a, .)

	di as text " "
	sum POM1 POM0 ps
	di as text " "

	// ATE
	gen double ATE = cond($flag == 1, (Q1star - Q0star), ///
		(Q1star - Q0star) * ($b - $a), .)
	qui sum ATE
	return scalar ATEtmle = r(mean)

	local RRtmle    = `Q1' / `Q0'
	local logRRtmle = log(`Q1') - log(`Q0')
	local ORtmle    = (`Q1' * (1 - `Q0')) / ((1 - `Q1') * `Q0')

	// Efficient Influence Curve for the ATE (canonical form):
	//   IC = (A/ps)*(Y - Q1*) + Q1* - mu1
	//        - { ((1-A)/(1-ps))*(Y - Q0*) + Q0* - mu0 }
	// computed on the [0,1] scale; for continuous Y multiply by (b-a) to
	// back-transform to the original outcome scale.
	gen d1 = (A * (Y - Q1star) / ps)             + Q1star - `Q1'
	gen d0 = ((1 - A) * (Y - Q0star) / (1 - ps)) + Q0star - `Q0'
	gen IC = cond($flag == 1, (d1 - d0), (d1 - d0) * ($b - $a), .)
	qui sum IC
	return scalar ATE_SE_tmle = sqrt(r(Var)/r(N))

	// Statistical inference ATE
	return scalar ATE_pvalue = 2 * normal(-abs(return(ATEtmle) / return(ATE_SE_tmle)))
	return scalar ATE_LCIa   = return(ATEtmle) - 1.96 * return(ATE_SE_tmle)
	return scalar ATE_UCIa   = return(ATEtmle) + 1.96 * return(ATE_SE_tmle)

	// Statistical inference for log-RR (sign matches R tmle::calcParameters)
	gen `ICrr' = (1/`Q1' * d1) - ((1/`Q0') * d0)
	qui sum `ICrr'
	local varICrr = r(Var)/r(N)
	local LCIrr   = exp(`logRRtmle' - 1.96 * sqrt(`varICrr'))
	local UCIrr   = exp(`logRRtmle' + 1.96 * sqrt(`varICrr'))

	// Statistical inference OR
	local ORtmle = (`Q1' / (1 - `Q1')) / (`Q0' / (1 - `Q0'))
	gen `ICor' = (1/(`Q1' * (1 - `Q1')) * d1) - (1/(`Q0' * (1 - `Q0')) * d0)
	qui sum `ICor'
	local varICor = r(Var)/r(N)
	local LCIOr = exp(log(`ORtmle') - 1.96 * sqrt(`varICor'))
	local UCIOr = exp(log(`ORtmle') + 1.96 * sqrt(`varICor'))

	// Return scalars
	return scalar CRR        = `RRtmle'
	return scalar SE_log_CRR = sqrt(`varICrr')
	return scalar MOR        = `ORtmle'
	return scalar SE_log_MOR = sqrt(`varICor')

	// Display ATE table
	if $flag == 1 {
		disp as text "{hline 63}"
		di "         {c |}" "    ATE         SE     P-value           95% CI"
		disp as text "{hline 63}"
		disp as text "TMLE:    {c |}" %7.4f as result return(ATEtmle) "    " ///
			%7.4f as result return(ATE_SE_tmle) "     " ///
			%7.4f as result return(ATE_pvalue) as text "     (" ///
			%7.4f as result return(ATE_LCIa) as text "," ///
			%7.4f as result return(ATE_UCIa) as text " )"
		disp as text "{hline 63}"
		disp as text " "
	}
	else {
		// Continuous outcomes: choose %g-style formatting that adapts to
		// magnitude (bweight-style hundreds and sim.do-style fractions).
		disp as text "{hline 71}"
		di "         {c |}" "      ATE         SE     P-value             95% CI"
		disp as text "{hline 71}"
		// %9.0g right-justifies the lower CI bound inside the trailing
		// width of the format spec, so leaving an extra leading space
		// after "(" or before ")" makes the bracket appear visually
		// detached.  Trim the trailing " )" -> ")" so the upper bound
		// abuts the closing bracket cleanly for continuous outcomes.
		disp as text "TMLE:    {c |}" %9.0g as result return(ATEtmle) "  " ///
			%9.0g as result return(ATE_SE_tmle) "    " ///
			%6.4f as result return(ATE_pvalue) as text "     (" ///
			%9.0g as result return(ATE_LCIa) as text "," ///
			%9.0g as result return(ATE_UCIa) as text ")"
		disp as text "{hline 71}"
		disp as text " "
	}

	// Display CRR and MOR table
	disp as text "{hline 51}"
	di "                           Estimate          95% CI"
	disp as text "{hline 51}"
	disp as text "Causal Risk Ratio:      " "{c |}      " ///
		%04.2f as result `RRtmle' as text "     (" ///
		%03.2f as result `LCIrr'  as text ","  ///
		%03.2f as result `UCIrr'  as text ")"
	disp as text "Marginal Odds Ratio:    " "{c |}      " ///
		%04.2f as result `ORtmle' as text "     (" ///
		%03.2f as result `LCIOr'  as text "," ///
		%03.2f as result `UCIOr'  as text ")"
	disp as text "{hline 51}"
end


// ---------------------------------------------------------------------------
// _eltmle_positivity_plot: density plot to check the positivity assumption
// Uses preserve/restore so the calling program's data is unaffected.
// ---------------------------------------------------------------------------
program _eltmle_positivity_plot
	preserve
	sort A
	qui by A: summarize ps
	qui kdensity ps if A==1, generate(x1pointsa d1A) nograph n(10000)
	qui kdensity ps if A==0, generate(x0pointsa d0A) nograph n(10000)
	label variable d1A "A = 1"
	label variable d0A "A = 0"
	set scheme sj
	twoway (line d0A x0pointsa) || ///
		(line d1A x1pointsa), ///
		xtitle("Propensity score") ///
		ytitle("Density") ///
		graphregion(color(white)) bgcolor(white) plotregion(fcolor(white)) ///
		legend(order(2 "Treated" 1 "Not treated") position(1) cols(1) ring(0) ///
			region(style(none)))
	restore
end


// ---------------------------------------------------------------------------
// _eltmle_covbal_table: standardised differences and variance ratio table
// Creates _ipw. Requires: ps in dataset, $variablelist global set.
// ---------------------------------------------------------------------------
program _eltmle_covbal_table
	// data2.dta is what is in memory here, and inside data2.dta the
	// outcome and exposure columns are always named "Y" and "A" (the
	// R-side script drops the user's original outcome/exposure columns
	// from fulldata before cbind, then adds explicit Y and A; see
	// _eltmle_write_r_noncv and changelog v4.0.2 Bug 0).  So we cannot
	// reference the user's original exposure name (e.g. "mbsmoke" with
	// cattaneo2.dta) here -- it does not exist in data2.dta.  Use "A".
	// $variablelist is parsed only to recover the covariate names for
	// the per-row table; outcome/exposure tokens are discarded.
	tokenize $variablelist
	macro shift                  // drop outcome
	macro shift                  // drop exposure
	local covars `*'

	// Inverse probability weights
	capture drop _ipw
	qui gen _ipw = .
	qui replace _ipw = (A==1) / ps       if A==1
	qui replace _ipw = (A==0) / (1 - ps) if A==0

	// Display table header
	di as text "{hline 67}"
	di as text "                 Standardised Differences            Variance ratio"
	di as text "                          Raw    Weighted           Raw    Weighted"
	di as text "{hline 67}"

	// Per-covariate statistics
	foreach var in `covars' {
		di as text "`var'"

		// Raw SMD
		qui summarize `var' if A==1
		local m1 = r(mean)
		local v1 = r(Var)
		qui summarize `var' if A==0
		local m0 = r(mean)
		local v0 = r(Var)
		local rSMD = (`m1' - `m0') / sqrt((`v1' + `v0') / 2)

		// Weighted SMD
		qui summarize `var' [iw=_ipw] if A==1
		local m1 = r(mean)
		local v1 = r(Var)
		qui summarize `var' [iw=_ipw] if A==0
		local m0 = r(mean)
		local v0 = r(Var)
		local wSMD = (`m1' - `m0') / sqrt((`v1' + `v0') / 2)

		// Raw VR
		qui sum `var' if A==1
		local v1 = r(Var)
		qui sum `var' if A==0
		local v0 = r(Var)
		local rVR = `v1' / `v0'

		// Weighted VR
		qui sum `var' [iw=_ipw] if A==1
		local v1 = r(Var)
		qui sum `var' [iw=_ipw] if A==0
		local v0 = r(Var)
		local wVR = `v1' / `v0'

		di as text "                    " ///
			%9.7g as result `rSMD' as text "   " ///
			%9.7g as result `wSMD' as text "     " ///
			%9.7g as result `rVR'  as text "   " ///
			%9.7g as result `wVR'
	}
	di as text "{hline 67}"
end


// ---------------------------------------------------------------------------
// _eltmle_label_rename_noncv: label / drop elements for non-CV programs
// Controlled by global $keepvars (0 = drop, 1 = label and prefix with _).
// ---------------------------------------------------------------------------
program _eltmle_label_rename_noncv
	if $keepvars == 0 {
		drop d1 d0 QAW Q1W Q0W Q1star Q0star ATE IC Y A POM1 POM0 ps
		capture drop ytempvar
	}
	if $keepvars == 1 {
		lab var d1      "Parameter for the influence curve"
		lab var d0      "Parameter for the influence curve"
		lab var QAW     "Initial prediction of the outcome"
		lab var Q1W     "Initial prediction of the outcome for A = 1"
		lab var Q0W     "Initial prediction of the outcome for A = 0"
		lab var Q1star  "Targeted update of Q for A = 1"
		lab var Q0star  "Targeted update of Q for A = 0"
		lab var A       "Exposure/Treatment"
		lab var Y       "Outcome"
		lab var ATE     "Average Treatment Effect"
		lab var IC      "Influence Curve"
		lab var POM1    "Potential Outcome Y(1)"
		lab var POM0    "Potential Outcome Y(0)"
		lab var ps      "Propensity Score"
		capture drop ytempvar
		foreach var of varlist d1 d0 QAW Q1W Q0W Q1star Q0star ATE IC Y A POM1 POM0 ps {
			rename `var' _`var'
		}
	}
end


// ---------------------------------------------------------------------------
// _eltmle_label_rename_cv: label / drop elements for CV programs
// Controlled by global $keepvars. CV variant also handles foldid, rowid,
// _ipw; with keepvars==1 saves elementsdata.dta for later merge.
// ---------------------------------------------------------------------------
program _eltmle_label_rename_cv
	if $keepvars == 0 {
		drop d1 d0 QAW Q1W Q0W Q1star Q0star ATE IC
		rename Y _Y
		rename A _A
		drop _Y _A
		capture drop _ipw
		drop POM1 POM0 ps
		capture drop cin
		drop foldid rowid
		capture drop x1pointsa x0pointsa d1A d0A
		capture drop ytempvar
	}
	if $keepvars == 1 {
		capture drop x1pointsa x0pointsa d1A d0A
		lab var d1      "Parameter for the influence curve"
		lab var d0      "Parameter for the influence curve"
		lab var QAW     "Initial prediction of the outcome"
		lab var Q1W     "Initial prediction of the outcome for A = 1"
		lab var Q0W     "Initial prediction of the outcome for A = 0"
		lab var Q1star  "Targeted update of Q for A = 1"
		lab var Q0star  "Targeted update of Q for A = 0"
		lab var A       "Exposure/Treatment"
		lab var Y       "Outcome"
		lab var ATE     "Average Treatment Effect"
		lab var IC      "Influence Curve"
		lab var POM1    "Potential Outcome Y(1)"
		lab var POM0    "Potential Outcome Y(0)"
		lab var ps      "Propensity Score"
		lab var _ipw    "Inverse probability of treatment weights"
		capture drop ytempvar
		foreach var of varlist d1 d0 QAW Q1W Q0W Q1star Q0star ATE IC Y A POM1 POM0 ps {
			rename `var' _`var'
		}
		keep rowid _d1 _d0 _QAW _Q1W _Q0W _Q1star _Q0star ///
			_ATE _IC _Y _A _POM1 _POM0 _ps _ipw
		qui save "elementsdata.dta", replace
	}
end


// ---------------------------------------------------------------------------
// _eltmle_cleanup_noncv: remove temp files for non-CV estimators
// ---------------------------------------------------------------------------
program _eltmle_cleanup_noncv
	quietly: rm SLS.R
	quietly: rm data2.dta
	quietly: rm data.csv
	quietly: rm fulldata.csv
	quietly: memory clean
end


// ---------------------------------------------------------------------------
// _eltmle_cleanup_cv: remove temp files for CV estimators
// ---------------------------------------------------------------------------
program _eltmle_cleanup_cv
	quietly: rm SLS.R
	quietly: rm data2.dta
	quietly: rm cvdata.csv
	quietly: rm cvdata.dta
	quietly: rm fulldata.csv
	forvalues i = 1(1)$folds {
		quietly: rm "folddata`i'.dta"
	}
	quietly: memory clean
end


// ====================================================
// NON-CV ESTIMATOR PROGRAMS
// Each collapses to: write R → run R → load data →
//   (optional bal helpers) → estimate → return → label/cleanup
// ====================================================

// ---------------------------------------------------------------------------
// _eltmle_make_ipw: generate inverse-probability weights using the exposure
// column "A" in data2.dta.  We CANNOT reference the user's original exposure
// name here (e.g. "mbsmoke" with cattaneo2.dta) -- R's cbind in
// _eltmle_write_r_{cv,noncv} drops the user's original outcome/exposure
// columns from fulldata before adding explicit Y and A columns, so inside
// the subprogram the only exposure column is named "A".  Requires `ps' in
// memory.
// ---------------------------------------------------------------------------
program _eltmle_make_ipw
	capture drop _ipw
	qui gen _ipw = .
	qui replace _ipw = (A==1) / ps       if A==1
	qui replace _ipw = (A==0) / (1 - ps) if A==0
end


// ---------------------------------------------------------------------------
// _eltmle_merge_elements_cv: post-restore cleanup for CV estimators.
// When keepvars=1, merges the saved prediction elements onto the user's
// original dataset by rowid (created in eltmle main, line 281), drops the
// _merge indicator, and removes the scratch elementsdata.dta file.  In all
// cases, drops foldid and rowid -- they are bookkeeping vars created in
// eltmle main that the user should never see in the returned dataset.
// ---------------------------------------------------------------------------
program _eltmle_merge_elements_cv
	if $keepvars == 1 {
		qui merge 1:1 rowid using "elementsdata.dta", nogenerate
		capture rm "elementsdata.dta"
	}
	qui drop foldid rowid
end


// ---------------------------------------------------------------------------
// _eltmle_save_elements_noncv: when keepvars=1, persist the renamed prediction
// vars (with row index) to elementsdata.dta so they can be merged back into
// the user's original dataset after restore.  Called *inside* the subprogram's
// preserve block, after _eltmle_label_rename_noncv has applied the _ prefix.
// ---------------------------------------------------------------------------
program _eltmle_save_elements_noncv
	if $keepvars == 1 {
		qui gen __eltmle_rowid = _n
		qui keep __eltmle_rowid _d1 _d0 _QAW _Q1W _Q0W _Q1star _Q0star ///
			_ATE _IC _Y _A _POM1 _POM0 _ps
		qui save "elementsdata.dta", replace
	}
end

// ---------------------------------------------------------------------------
// _eltmle_merge_elements_noncv: after restore, merge the saved prediction
// elements onto the user's original dataset by row index, then delete the
// scratch file.  No-op when keepvars=0.
// ---------------------------------------------------------------------------
program _eltmle_merge_elements_noncv
	if $keepvars == 1 {
		qui gen __eltmle_rowid = _n
		qui merge 1:1 __eltmle_rowid using "elementsdata.dta", nogenerate
		qui drop __eltmle_rowid
		capture rm "elementsdata.dta"
	}
end


program tmle, rclass
preserve
	_eltmle_write_r_noncv base
	_eltmle_run_r
	qui clear
	quietly: use "data2.dta", clear
	qui cap drop X__000000
	_eltmle_tmle_estimate
	return add
	_eltmle_label_rename_noncv
	_eltmle_save_elements_noncv
restore
	_eltmle_cleanup_noncv
	_eltmle_merge_elements_noncv
end

program tmlebal, rclass
preserve
	_eltmle_write_r_noncv base
	_eltmle_run_r
	qui clear
	quietly: use "data2.dta", clear
	qui cap drop X__000000
	_eltmle_positivity_plot
	_eltmle_tmle_estimate
	return add
	_eltmle_covbal_table
	_eltmle_label_rename_noncv
	_eltmle_save_elements_noncv
restore
	_eltmle_cleanup_noncv
	_eltmle_merge_elements_noncv
end

program tmlebgam, rclass
preserve
	_eltmle_write_r_noncv bgam
	_eltmle_run_r
	qui clear
	quietly: use "data2.dta", clear
	qui cap drop X__000000
	_eltmle_tmle_estimate
	return add
	_eltmle_label_rename_noncv
	_eltmle_save_elements_noncv
restore
	_eltmle_cleanup_noncv
	_eltmle_merge_elements_noncv
end

program tmlebgambal, rclass
preserve
	_eltmle_write_r_noncv bgam
	_eltmle_run_r
	qui clear
	quietly: use "data2.dta", clear
	qui cap drop X__000000
	_eltmle_positivity_plot
	_eltmle_tmle_estimate
	return add
	_eltmle_covbal_table
	_eltmle_label_rename_noncv
	_eltmle_save_elements_noncv
restore
	_eltmle_cleanup_noncv
	_eltmle_merge_elements_noncv
end

program tmleglsrf, rclass
preserve
	_eltmle_write_r_noncv glsrf
	_eltmle_run_r
	qui clear
	quietly: use "data2.dta", clear
	qui cap drop X__000000
	_eltmle_tmle_estimate
	return add
	_eltmle_label_rename_noncv
	_eltmle_save_elements_noncv
restore
	_eltmle_cleanup_noncv
	_eltmle_merge_elements_noncv
end

program tmleglsrfbal, rclass
preserve
	_eltmle_write_r_noncv glsrf
	_eltmle_run_r
	qui clear
	quietly: use "data2.dta", clear
	qui cap drop X__000000
	_eltmle_positivity_plot
	_eltmle_tmle_estimate
	return add
	_eltmle_covbal_table
	_eltmle_label_rename_noncv
	_eltmle_save_elements_noncv
restore
	_eltmle_cleanup_noncv
	_eltmle_merge_elements_noncv
end


// ====================================================
// CROSS-VALIDATION ESTIMATOR PROGRAMS
// Structure: preserve → fold loop (write R, run R, save folddata) →
//   combine folds → TMLE estimate → covbal → label → restore →
//   cleanup → merge elements
// ====================================================

program cvtmle, rclass
preserve
	forval fold = 1/$folds {

		if `fold' == 1 _dots 0, title(Performing $folds folds)

		qui use "cvdata.dta", clear
		capture drop foldnumber
		qui gen foldnumber = `fold'
		qui: export delimited $varusedforcv foldnumber using "cvdata.csv", nolabel replace

		_eltmle_write_r_cv base
		_eltmle_run_r

		quietly: use "data2.dta", clear
		qui save "folddata`fold'.dta", replace

		_dots `fold' 0
	}

	// Combine fold predictions
	qui: use "folddata1.dta", clear
	forval fold = 2/$folds {
		qui: append using "folddata`fold'"
	}
	qui save "data2.dta", replace

	// Load combined data and run TMLE estimation
	qui clear
	quietly: use "data2.dta", clear
	qui cap drop X__000000
	_eltmle_tmle_estimate
	return add

	// IPW weights (always); covbal table and positivity plot (conditional).
	// _eltmle_make_ipw uses the column "A" in data2.dta -- see helper docs
	// for why we cannot use the user's original exposure name here.
	if $covbalancetable == 0 {
		_eltmle_make_ipw
	}
	if $covbalancetable == 1 {
		_eltmle_covbal_table
		// Positivity plot — inline (no nested preserve; already inside one)
		qui kdensity ps if A==1, generate(x1pointsa d1A) nograph n(10000)
		qui kdensity ps if A==0, generate(x0pointsa d0A) nograph n(10000)
		label variable d1A "A = 1"
		label variable d0A "A = 0"
		set scheme sj
		twoway (line d0A x0pointsa) || ///
			(line d1A x1pointsa), ///
			xtitle("Propensity score") ///
			ytitle("Density") ///
			graphregion(color(white)) bgcolor(white) plotregion(fcolor(white)) ///
			legend(order(2 "Treated" 1 "Not treated") position(1) cols(1) ring(0) ///
				region(style(none)))
	}

	_eltmle_label_rename_cv

restore

	capture drop ytempvar
	_eltmle_cleanup_cv
	_eltmle_merge_elements_cv
end

program cvtmlebgam, rclass
preserve
	forval fold = 1/$folds {

		if `fold' == 1 _dots 0, title(Performing $folds folds)

		qui use "cvdata.dta", clear
		capture drop foldnumber
		qui gen foldnumber = `fold'
		qui: export delimited $varusedforcv foldnumber using "cvdata.csv", nolabel replace

		_eltmle_write_r_cv bgam
		_eltmle_run_r

		quietly: use "data2.dta", clear
		qui save "folddata`fold'.dta", replace

		_dots `fold' 0
	}

	// Combine fold predictions
	qui: use "folddata1.dta", clear
	forval fold = 2/$folds {
		qui: append using "folddata`fold'"
	}
	qui save "data2.dta", replace

	// Load combined data and run TMLE estimation
	qui clear
	quietly: use "data2.dta", clear
	qui cap drop X__000000
	_eltmle_tmle_estimate
	return add

	// IPW weights (always); covbal table and positivity plot (conditional).
	// _eltmle_make_ipw uses the column "A" in data2.dta -- see helper docs
	// for why we cannot use the user's original exposure name here.
	if $covbalancetable == 0 {
		_eltmle_make_ipw
	}
	if $covbalancetable == 1 {
		_eltmle_covbal_table
		qui kdensity ps if A==1, generate(x1pointsa d1A) nograph n(10000)
		qui kdensity ps if A==0, generate(x0pointsa d0A) nograph n(10000)
		label variable d1A "A = 1"
		label variable d0A "A = 0"
		set scheme sj
		twoway (line d0A x0pointsa) || ///
			(line d1A x1pointsa), ///
			xtitle("Propensity score") ///
			ytitle("Density") ///
			graphregion(color(white)) bgcolor(white) plotregion(fcolor(white)) ///
			legend(order(2 "Treated" 1 "Not treated") position(1) cols(1) ring(0) ///
				region(style(none)))
	}

	_eltmle_label_rename_cv

restore

	capture drop ytempvar
	_eltmle_cleanup_cv
	_eltmle_merge_elements_cv
end

program cvtmleglsrf, rclass
preserve
	forval fold = 1/$folds {

		if `fold' == 1 _dots 0, title(Performing $folds folds)

		qui use "cvdata.dta", clear
		capture drop foldnumber
		qui gen foldnumber = `fold'
		qui: export delimited $varusedforcv foldnumber using "cvdata.csv", nolabel replace

		_eltmle_write_r_cv glsrf
		_eltmle_run_r

		quietly: use "data2.dta", clear
		qui save "folddata`fold'.dta", replace

		_dots `fold' 0
	}

	// Combine fold predictions
	qui: use "folddata1.dta", clear
	forval fold = 2/$folds {
		qui: append using "folddata`fold'"
	}
	qui save "data2.dta", replace

	// Load combined data and run TMLE estimation
	qui clear
	quietly: use "data2.dta", clear
	qui cap drop X__000000
	_eltmle_tmle_estimate
	return add

	// IPW weights (always); covbal table and positivity plot (conditional).
	// _eltmle_make_ipw uses the column "A" in data2.dta -- see helper docs
	// for why we cannot use the user's original exposure name here.
	if $covbalancetable == 0 {
		_eltmle_make_ipw
	}
	if $covbalancetable == 1 {
		_eltmle_covbal_table
		qui kdensity ps if A==1, generate(x1pointsa d1A) nograph n(10000)
		qui kdensity ps if A==0, generate(x0pointsa d0A) nograph n(10000)
		label variable d1A "A = 1"
		label variable d0A "A = 0"
		set scheme sj
		twoway (line d0A x0pointsa) || ///
			(line d1A x1pointsa), ///
			xtitle("Propensity score") ///
			ytitle("Density") ///
			graphregion(color(white)) bgcolor(white) plotregion(fcolor(white)) ///
			legend(order(2 "Treated" 1 "Not treated") position(1) cols(1) ring(0) ///
				region(style(none)))
	}

	_eltmle_label_rename_cv

restore

	capture drop ytempvar
	_eltmle_cleanup_cv
	_eltmle_merge_elements_cv
end

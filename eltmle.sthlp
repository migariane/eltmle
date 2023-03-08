{smcl}
{right: version 2.2.8  01.03.2023}
{...}

{title:Title}

{p2colset 5 18 20 2}{...}
{p2col:{hi:eltmle}{hline 1}}Ensemble Learning Targeted Maximum Likelihood Estimation{p2colreset}{...}

{title:Syntax}

{p 4 4 2}
{cmd:eltmle} {hi: Y} {hi: X} {hi: Z} [{cmd:,} {it:tmle} {it:tmlebgam} {it:tmleglsrf} {it:bal} {it:elements}]

{p 4 4 2}
where:
{p_end}

{p 4 4 2}
{hi:Y}: Outcome: numeric binary or continuous variable
{p_end}
{p 4 4 2}
{hi:X}: Treatment or exposure: numeric binary variable
{p_end}
{p 4 4 2}
{hi:Z}: Covariates: vector of continuous and categorical variables
{p_end}

{title:Description}

{p 4 4 2 120}
{hi: Modern Epidemiology} has been able to identify significant limitations of classic epidemiological methods, like outcome regression analysis, when estimating causal quantities such as the average treatment effect (ATE)
for observational data. For example, using classical regression models to estimate the ATE requires making the assumption that the effect measure is constant across levels of confounders included in the model,
 i.e. that there is no effect modification. Other methods do not require this assumption, including g-methods (e.g. the {hi:g-formula}) and targeted maximum likelihood estimation ({hi:TMLE}).
{p_end}

{p 4 4 2 120}
The average treatment effect ({hi:ATE}) or risk difference is the most commonly used causal parameter. Many estimators of the ATE but no all rely on parametric modeling assumptions. Therefore, the correct model specification is crucial
to obtain unbiased estimates of the true ATE.
{p_end}

{p 4 4 2 120}
TMLE is a semiparametric, efficient substitution estimator allowing for data-adaptive estimation while obtaining valid statistical inference based on the targeted minimum loss-based estimation. TMLE has the advantage of
being doubly robust. Moreover, TMLE allows inclusion of {hi:machine learning} algorithms to minimise the risk of model misspecification, a problem that persists for competing estimators. Evidence shows that TMLE typically
provides the {hi: least unbiased} estimates of the ATE compared with other double robust estimators.
{p_end}

{p 4 4 2 120}
The following link provides access to a TMLE tutorial: {browse "http://migariane.github.io/TMLE.nb.html":TMLE_tutorial}.
{p_end}

{p 4 4 2 120}
{hi:eltmle} is a Stata program implementing the targeted maximum likelihood estimation for the ATE for a binary or continuous outcome and binary treatment. {hi:eltmle} includes the use of a super-learner called from the {hi:SuperLearner}
package v.2.0-21 (Polley E., et al. 2011). The Super-Learner uses V-fold cross-validation (10-fold by default) to assess the performance of prediction regarding the potential outcomes and the propensity score as weighted
averages of a set of machine learning algorithms. We used the default SuperLearner algorithms implemented in the base installation of the {hi:tmle-R} package v.1.2.0-5 (Susan G. and Van der Laan M., 2017),
which included the following: i) stepwise selection, ii) generalized linear modeling (GLM), iii) a GLM variant that includes second order polynomials and two-by-two interactions of the main terms
included in the model. Additionally, {hi:eltmle} users will have the option to include Bayes Generalized Linear Models and Generalized Additive Models as additional Super-Learner algorithms. Future implementations will offer
more advanced machine learning algorithms.
{p_end}

{title:Options}

{p 4 4 2 120}
{hi:tmle}: this is the default option. If no-option is specified eltmle by default implements the
TMLE algorithm plus the super-Learner ensemble learning for the main three machine learning algorithms described above.
{p_end}

{p 4 4 2 120}
{hi:tmlebgam}: this option may be specified or unspecified. When specified, it does include in addition to the above default
implementation, the Bayes Generalized Linear Models and Generalized Additive Models as Super-Learner algorithms for the tmle estimator.
This option might be suitable for non-linear treatment effects.
{p_end}

{p 4 4 2 120}
{hi:tmleglsrf}: this option may be specified or unspecified. When specified, it does include in addition to the three main learning algorithms
described above, the Lasso (glmnet R package), Random Forest (randomForest R package) and the Generalized Additive Models as Super-Learner algorithms for the tmle estimator.
This option might be suitable for heterogeneous treatment effects.
{p_end}

{p 4 4 2 120}
{hi:bal}: this option may be specified or unspecified. When specified, it provides two additional features. Firstly, a visual diagnostic check of the positivity assumption
based on the estimation of kernel density plots for the propensity score by levels of the treatment. Secondly, a table displaying the differences in distributions of each of 
the covariates Z between treatment groups: Standardised mean differences and variance ratios are reported, for both the raw and weighted covariate values. Note that perfect 
covariate balance between treatment groups is indicated by Standardised Mean Differences of 0 and variance ratios of 1. Both are calculated using formulas from Austin (2009) 
{it: Balance Diagnostics for Comparing the Distribution of Baseline Covariates Between Treatment Groups in Propensity-Score Matched Samples}.
{p_end}

{p 4 4 2 120}
{hi:elements}: this option may be specified or unspecified. When specified, the data set will retain the variables used for each step in TMLE such as the initial predictions, propensity score, updated predictions, etc.
{p_end}


{title:Results}

{p 4 4 2 120}
In addition to the ATE, the ATE's standard error and p-value, the marginal odds ratio (MOR), and the causal risk ratio (CRR),
including their respective type Wald 95%CIs, {hi:eltmle} output provides a descriptive summary for the potential outcomes (POM)
and the propensity score (ps):
{hi: POM1}: Potential outcome among the treated
{hi: POM0}: Potential outcome among the non-treated
{hi: ps}: Propensity score


{title:Example}

We provide the following examples:
1) TMLE for:
  a) Binary outcome
  b) Continuous outcome
2) Advanced machine-learning techniques:
  a) tmleglsrf
  b) tmlebgam
3) Retaining potential outcomes using the {hi:elements} option
4) Covariate balance tables

***********************************************************
* eltmle Y X Z [if] [,tmle tmlebgam tmleglsrf bal elements]
***********************************************************

.clear
.use http://www.stata-press.com/data/r14/cattaneo2.dta
.describe
.gen lbw = cond(bweight<2500,1,0.)
.lab var lbw "Low birthweight, <2500 g"
.save "your path/cattaneo2.dta", replace
.cd "your path"


**********************************
* 1) a) TMLE with a binary outcome
**********************************

.eltmle lbw mbsmoke mage medu prenatal mmarried, tmle


 Variable |        Obs        Mean    Std. dev.       Min        Max
----------+---------------------------------------------------------
     POM1 |      4,642    .1025904    .0404756   .0197383   .3806479
     POM0 |      4,642    .0515831    .0253505   .0212925   .1713066
       ps |      4,642    .1861267    .1111685   .0356226     .85422

---------------------------------------------------------------
         |    ATE         SE     P-value           95% CI
---------------------------------------------------------------
TMLE:    | 0.0510     0.0122      0.0001     ( 0.0270, 0.0750 )
---------------------------------------------------------------

---------------------------------------------------
                           Estimate          95% CI
---------------------------------------------------
Causal Risk Ratio:      |      2.00     (1.53,2.60)
Marginal Odds Ratio:    |      2.11     (1.50,2.72)
---------------------------------------------------



**************************************
* 1) b) TMLE with a continuous outcome
**************************************

.eltmle bweight mbsmoke mage medu prenatal mmarried, tmle


 Variable |        Obs        Mean    Std. dev.       Min        Max
----------+---------------------------------------------------------
     POM1 |      4,642    2832.914    74.79904   2574.164    2961.55
     POM0 |      4,642     3062.81    89.92703   2863.544   3169.349
       ps |      4,642    .1861267    .1111685   .0356226     .85422

---------------------------------------------------------------
         |    ATE         SE     P-value           95% CI
---------------------------------------------------------------
TMLE:    | -229.9       24.6      0.0000     ( -278.2, -181.6 )
---------------------------------------------------------------

---------------------------------------------------
                           Estimate          95% CI
---------------------------------------------------
Causal Risk Ratio:      |      0.93     (0.91,0.94)
Marginal Odds Ratio:    |      0.83     (0.80,0.87)
---------------------------------------------------



*****************
* 2) a) tmleglsrf
*****************

.eltmle lbw mbsmoke mage medu prenatal mmarried, tmleglsrf


 Variable |        Obs        Mean    Std. dev.       Min        Max
----------+---------------------------------------------------------
     POM1 |      4,642    .0804235    .0392915   .0201732   .4046224
     POM0 |      4,642    .0399219    .0217449   .0150639   .2473742
       ps |      4,642    .1536965    .1122575       .025    .635824

---------------------------------------------------------------
         |    ATE         SE     P-value           95% CI
---------------------------------------------------------------
TMLE:    | 0.0405     0.0146      0.0171     ( 0.0119, 0.0691 )
---------------------------------------------------------------

---------------------------------------------------
                           Estimate          95% CI
---------------------------------------------------
Causal Risk Ratio:      |      1.68     (1.16,2.44)
Marginal Odds Ratio:    |      1.74     (1.05,2.44)
---------------------------------------------------


*****************
* 2) b) tmlebgam
*****************

.eltmle lbw mbsmoke mage medu prenatal mmarried, tmlebgam


 Variable |        Obs        Mean    Std. dev.       Min        Max
----------+---------------------------------------------------------
     POM1 |      4,642    .1024584     .040525   .0200955   .3841472
     POM0 |      4,642     .051588    .0253702    .020537   .1600263
       ps |      4,642    .1861267    .1105893   .0313728   .5872202

---------------------------------------------------------------
         |    ATE         SE     P-value           95% CI
---------------------------------------------------------------
TMLE:    | 0.0509     0.0125      0.0002     ( 0.0263, 0.0754 )
---------------------------------------------------------------

---------------------------------------------------
                           Estimate          95% CI
---------------------------------------------------
Causal Risk Ratio:      |      2.00     (1.53,2.62)
Marginal Odds Ratio:    |      2.11     (1.49,2.74)
---------------------------------------------------


*************
* 3) Elements
*************

.eltmle lbw mbsmoke mage medu prenatal mmarried, elements


*****************************
* 4) Covariates balance table
*****************************

.eltmle lbw mbsmoke mage medu prenatal mmarried, bal


 Variable |        Obs        Mean    Std. dev.       Min        Max
----------+---------------------------------------------------------
     POM1 |      4,642    .1025904    .0404756   .0197383   .3806479
     POM0 |      4,642    .0515831    .0253505   .0212925   .1713066
       ps |      4,642    .1861267    .1111685   .0356226     .85422

---------------------------------------------------------------
         |    ATE         SE     P-value           95% CI
---------------------------------------------------------------
TMLE:    | 0.0510     0.0122      0.0001     ( 0.0270, 0.0750 )
---------------------------------------------------------------

---------------------------------------------------
                           Estimate          95% CI
---------------------------------------------------
Causal Risk Ratio:      |      2.00     (1.53,2.60)
Marginal Odds Ratio:    |      2.11     (1.50,2.72)
---------------------------------------------------

-------------------------------------------------------------------
             Standardised Differences            Variance ratio
                      Raw    Weighted           Raw    Weighted
-------------------------------------------------------------------
mage
                 -.300179   -.0211604      .8818025    .8687959
medu
                -.5474357   -.0947559      .7315846    .5692243
prenatal
                 .2339922     .010213      1.774373    1.108046
mmarried
                -.5953009   -.0172262      1.335944    1.014663
-------------------------------------------------------------------



**********************************************************************************************

{title:Remarks}

{p 4 4 2 120}
Remember 1: Y must be a binary or continuous variable; X must be numeric binary
variable coded (0, 1) and, Z a vector of numeric covariates.
{p_end}

{p 4 4 2 120}
Remember 2: You must change your working directory to the location of the Stata.dta file.
{p_end}

{p 4 4 2 120}
Remember 3: eltmle automatically implements a complete case analysis (i.e., listwise delection for missing data). However, if you would like to impute your missing values 
before running eltmle, see the example here below using the sys auto data:
{p_end}

		{title:Example of imputation using predictive mean matching}
		.clear
		.sysuse auto
		.describe
		.misstable summarize
		.mi set wide
		.mi register imputed rep78
    // Impute using predictive mean matching: only one dataset and knn(# nearest neighbors) to draw from
		.mi impute pmm rep78, add(1) knn(5)
		.describe
		.drop _mi_miss rep78
		.rename _1_rep78 rep78
		.eltmle price foreign rep78 weight length, tmle

{p 4 4 2 120}
Remember 4: Mac users must have installed R software on their personal computer as
eltmle calls R to implement the Super Learner. The R executable file must be located at
the following path: {hi:"/usr/local/bin/r"}.
{p_end}

{p 4 4 2 120}
Remember 5: Windows users must have installed R software on their personal computer
as eltmle calls R to implement the Super Learner. The R executable file must be located
at the following path: {hi:"C:\Program Files\R\R-X.X.X\bin\x64\R.exe"} (where X stands for the number of the version).
{p_end}

{p 4 4 2 120}
Remember 6: Windows users must have only one version of R software installed on their personal computer
at the following path: {hi:"C:\Program Files\R\R-X.X.X\bin\x64\R.exe"}. In case more than one different version
is located in the above highlighted path users would like to keep the latest.
{p_end}


{title:Stored results}

eltmle stores the following in {hi: r()}:

Scalars

	{hi: r(SE_log_MOR)}          Standard error marginal odds ratio
               {hi: r(MOR)} 		Marginal odds ratio
        {hi: r(SE_log_CRR)} 		Standard error causal risk ratio
        {hi:        r(CRR)} 		Causal risk ratio
        {hi:   r(ATE_UCIa)} 		Risk difference upper 95%CI
        {hi:   r(ATE_LCIa)} 		Risk difference lower 95%CI
        {hi: r(ATE_pvalue)} 		Risk difference pvalue
        {hi: (ATE_SE_tmle)} 		Standard error Risk difference
        {hi:    r(ATEtmle)} 		Risk difference

{title:Version in development: updates}

{browse "https://github.com/migariane/eltmle/tree/master": https://github.com/migariane/eltmle/tree/master}

{title:References}

{p 4 4 2 120}
Miguel Angel Luque‐Fernandez, M Schomaker, B Rachet, M Schnitzer (2018). Targeted maximum likelihood estimation for a binary treatment: A tutorial.
Statistics in medicine. {browse "https://onlinelibrary.wiley.com/doi/abs/10.1002/sim.7628":link}.
{p_end}

{p 4 4 2 120}
Miguel Angel Luque-Fernandez (2017). Targeted Maximum Likelihood Estimation for a
Binary Outcome: Tutorial and Guided Implementation {browse "http://migariane.github.io/TMLE.nb.html":link}.
{p_end}

{p 4 4 2 120}
Matthew James Smith, Mohammad A. Mansournia, Camille Maringe, Clemence Leyrat, Aurelien Belot, Bernard Rachet, Paul Zivich, Stephen R. Cole, Miguel Angel Luque Fernandez (2021). Tutorial: Introduction to computational causal inference for applied researchers and epidemiologists {browse "https://github.com/migariane/TutorialCausalInferenceEstimators":link}.
{p_end}

{p 4 4 2 120}
StataCorp. 2015. Stata Statistical Software: Release 14. College Station, TX: StataCorp LP.
{p_end}

{p 4 4 2 120}
Gruber S, Laan M van der. (2011). Tmle: An R package for targeted maximum likelihood
estimation. UC Berkeley Division of Biostatistics Working Paper Series.
{p_end}

{p 4 4 2 120}
Laan M van der, Rose S. (2011). Targeted learning: Causal inference for observational
and experimental data. Springer Series in Statistics.626p.
{p_end}

{p 4 4 2 120}
Van der Laan MJ, Polley EC, Hubbard AE. (2007). Super learner. Statistical applications
in genetics and molecular biology 6.
{p_end}

{title:Author and developers}

{phang}Miguel Angel Luque-Fernandez [Author and developer] {p_end}
{phang}Department of Statistics and Operations Research (Biostatistics group) University of Granada,
Granada, Spain.{p_end}
{phang}Department of Epidemiology and Population Health, ICON group, London School of Hygiene and Tropical Medicine. London, UK.{p_end}
{phang}E-mail: {browse "mailto:miguel-angel.luque@lshtm.ac.uk":miguel-angel.luque@lshtm.ac.uk}{p_end}

{phang}Camille Maringe [Developer] {p_end}
{phang}Department of Epidemiology and Population Health, ICON group, London School of Hygiene and Tropical Medicine. London, UK.{p_end}
{phang}E-mail: {browse "mailto:camille.maringe@lshtm.ac.uk":camille.maringe@lshtm.ac.uk}{p_end}

{phang}Matthew J. Smith [Developer] {p_end}
{phang}Department of Epidemiology and Population Health, ICON group, London School of Hygiene and Tropical Medicine. London, UK.{p_end}
{phang}E-mail: {browse "mailto:matthew1.smith@lshtm.ac.uk":matthew1.smith@lshtm.ac.uk}{p_end}

{title:Also see}

{psee}
Online:  {help teffects}
{p_end}

{smcl}
{right: version 2.2.2 September 17th, 2018}
{...}

{title:Title}

{p2colset 5 18 20 2}{...}
{p2col:{hi:eltmle}{hline 1}}Ensemble Learning Targeted Maximum Likelihood Estimation{p2colreset}{...}

{title:Syntax}

{p 4 4 2}
{cmd:eltmle} {hi: Y} {hi: X} {hi: Z} [{cmd:,} {it:tmle} {it:tmlebgam} {it:tmleglsrf}]

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
{hi:Z}: Covariates: vector of numeric and categorical variables
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
{p_end}

{p 4 4 2 120}
{hi:tmleglsrf}: this option may be specified or unspecified. When specified, it does include in addition to the three main learning algorithms 
described above, the Lasso (glmnet R package), Random Forest (randomForest R package) and the Generalized Additive Models as Super-Learner algorithms for the tmle estimator.
{p_end}

{title:Example}

**********************************************
* eltmle Y X Z [if] [,tmle tmlebgam tmleglsrf] 
**********************************************

.clear
.use http://www.stata-press.com/data/r14/cattaneo2.dta
.describe
.gen lbw = cond(bweight<2500,1,0.)
.lab var lbw "Low birthweight, <2500 g"
.save "your path/cattaneo2.dta", replace
.cd "your path"

******************
// Binary outcome 
******************

******************************************************
.eltmle lbw mbsmoke mage medu prenatal mmarried, tmle
******************************************************

    Variable |        Obs        Mean    Std. Dev.       Min        Max
-------------+---------------------------------------------------------
        POM1 |      4,642    .1023406    .0401616   .0201151   .3747893
        POM0 |      4,642     .051377    .0251473   .0208754   .1706158
          PS |      4,642    .1861267     .110755   .0372202   .8494988
		  
________________________________
TMLE: Average Treatment Effect
________________________________

Risk Differences: 0.05; SE: 0.0122; p-value: 0.0001; 95%CI:(0.03, 0.07)

________________________________
TMLE: Causal Relative Risk (CRR)
________________________________

CRR:     1.99; 95%CI:(1.53,2.59)

________________________________
TMLE: Marginal Odds Ratio (MOR)
________________________________

MOR:     2.11; 95%CI:(1.50,2.72)

**********************
// Continuous outcome 
**********************

***********************************************************
.eltmle bweight mbsmoke mage medu prenatal mmarried, tmle
***********************************************************

    Variable |        Obs        Mean    Std. Dev.       Min        Max
-------------+---------------------------------------------------------
        POM1 |      4,642    2832.925    74.85617   2579.981   2958.981
        POM0 |      4,642    3062.132    89.56902   2864.008   3166.092
          PS |      4,642    .1861267     .110755   .0372202   .8494988
		  
_____________________________
TMLE: Additive Causal Effect
____________________________

Risk Differences: -229.21; SE: 0.3413; p-value: 0.0000; 95%CI:(-229.88, -228.54)

________________________________
TMLE: Causal Relative Risk (CRR)
________________________________

CRR: 0.93; 95%CI:(0.91,0.94)

________________________________
TMLE: Marginal Odds Ratio (MOR)
________________________________

MOR: 0.83; 95%CI:(0.80,0.87)

***************************************************************
// Continuous outcome: preserving original dataset and using 
// more advance machine learning techniques
***************************************************************

***************************************************************
.preserve
.eltmle bweight mbsmoke mage medu prenatal mmarried, tmleglsrf
.restore
***************************************************************

    Variable |        Obs        Mean    Std. Dev.       Min        Max
-------------+---------------------------------------------------------
        POM1 |      4,642    2834.385    74.96674    2582.38   2967.601
        POM0 |      4,642    3063.005     89.5585   2867.587    3167.28
          PS |      4,642     .154641    .1111959       .025   .6236143
		  
________________________________
TMLE: Additive Causal Effect
________________________________

Risk Differences:   -228.62; SE: 0.4193; p-value: 0.0000; 95%CI:(-229.44,-227.80)

________________________________
TMLE: Causal Relative Risk (CRR)
________________________________

CRR:     0.93; 95%CI:(0.91,0.94)

________________________________
TMLE: Marginal Odds Ratio (MOR)
________________________________

MOR:     0.83; 95%CI:(0.80,0.87)

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
Remember 3: Mac users must have installed R software on their personal computer as
eltmle calls R to implement the Super Learner. The R executable file must be located at 
the following path: {hi:"/usr/local/bin/r"}.
{p_end}

{p 4 4 2 120}
Remember 4: Windows users must have installed R software on their personal computer
as eltmle calls R to implement the Super Learner. The R executable file must be located 
at the following path: {hi:"C:\Program Files\R\R-X.X.X\bin\x64\R.exe"} (where X stands for the number of the version).
{p_end}

{p 4 4 2 120}
Remember 5: Windows users must have only one version of R software installed on their personal computer
at the following path: {hi:"C:\Program Files\R\R-X.X.X\bin\x64\R.exe"}. In case more than one different version 
is located in the above highlighted path users would like to keep the latest.
{p_end}

{p 4 4 2 120}
Remember 6: In case you want to preserve the original dataset you can use the 
preserve restore Stata functionality in combination to the Stata {hi: eltmle} command as shown in the previous example.
{p_end}

{title:References}

{p 4 4 2 120}
Miguel Angel Luque‐Fernandez, M Schomaker, B Rachet, M Schnitzer (2018). Targeted maximum likelihood estimation for a binary treatment: A tutorial.
Statistics in medicine. {browse "https://onlinelibrary.wiley.com/doi/abs/10.1002/sim.7628":Download here}.
{p_end}

{p 4 4 2 120}
Miguel Angel Luque-Fernandez (2017). Targeted Maximum Likelihood Estimation for a
Binary Outcome: Tutorial and Guided Implementation {browse "http://migariane.github.io/TMLE.nb.html":Download here}.
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

{title:Author and developer}

{phang}Miguel Angel Luque-Fernandez{p_end}
{phang}Biomedical Research Institute of Granada, Noncommunicable Disease and Cancer Epidemiolgy Group. University of Granada,
Granada, Spain.{p_end}
{phang}Department of Epidemiology and Population Health, London School of Hygiene and Tropical Medicine. 
London, UK.{p_end}
{phang}E-mail: {browse "mailto:miguel-angel.luque@lshtm.ac.uk":miguel-angel.luque@lshtm.ac.uk}{p_end}  

{title:Also see}

{psee}
Online:  {helpb teffects}
{p_end}

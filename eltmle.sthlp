{smcl}
{right:*!version 1.7 05.APRIL.2017}{...}

{phang}
{cmd:help eltmle}
{hline}

{title:Title}

{p2colset 5 18 20 2}{...}
{p2col:{hi:eltmle}{hline 1}}Ensemble Learning Targeted Maximum Likelihood Estimation{p2colreset}{...}

{title:Syntax}

{p 4 4 2}
{cmd:eltmle} {hi: Y} {hi: X} {hi: Z} [{cmd:,} {it:slapiw} {it:slaipwbgam} {it:tmle} {it:tmlebgam}]

{p 4 4 2}
where:
{p_end}

{p 4 4 2}
{hi:Y}: Outcome: numeric binary variable
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
{hi:eltmle} is a Stata program implementing the targeted maximum likelihood estimation for the ATE for a binary outcome and binary treatment. {hi:eltmle} includes the use of a super-learner called from the {hi:SuperLearner}
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
{hi:slaipw}: this option may be specified or unspecified. When specified, it does estimate the augmented
inverse probability weighting algorithm plus the Super Learner ensemble learning for the main three machine 
learning algorithms described above.
{p_end}

{p 4 4 2 120}
{hi:slaipwbgam}: this option may be specified or unspecified. When specified, it does include in addition to the above default
implementation, the Bayes Generalized Linear Models and Generalized Additive Models as Super-Learner algorithms for the slaipw estimator.
{p_end}

{title:Example}

*******************************************************
* eltmle Y X Z [if] [,slaipw slaipwbgam tmle tmlebgam] 
*******************************************************

.clear
.use http://www.stata-press.com/data/r14/cattaneo2.dta
.describe
.gen lbw = cond(bweight<2500,1,0.)
.lab var lbw "Low birthweight, <2500 g"
.save "your path/cattaneo2.dta", replace
.cd "your path"

*******************************************************
.eltmle lbw mbsmoke mage medu prenatal mmarried, tmle
*******************************************************

    Variable |        Obs        Mean    Std. Dev.       Min        Max
-------------+---------------------------------------------------------
      Q1star |      4,642    .1055429    .0403702   .0205745   .3789994
      Q0star |      4,642    .0509794    .0249112   .0207113   .1664136
          ps |      4,642    .1861267     .110755   .0372202   .8494988


TMLE: Average Treatment Effect

ATE:  0.0546; SE:0.0122; p-value: 0.0000; 95%CI:(0.030700,0.078427)

TMLE: Relative Risk

RR:   2.0703; 95%CI:(2.0132,2.1290)

*********************************************************
.eltmle lbw mbsmoke mage medu prenatal mmarried, slaipw
*********************************************************

    Variable |        Obs        Mean    Std. Dev.       Min        Max
-------------+---------------------------------------------------------
        aQ1W |      4,642    .0959911    1.717125  -2.481395     20.267
        aQ0W |      4,642    .0516773    .3212793  -3.679752   1.911056
          ps |      4,642    .1861267     .110755   .0372202   .8494988

AIPW ensemble learning: Average Treatment Effect

ATE:  0.0542; SE:0.0122; p-value: 0.0000; 95%CI:(0.030299,0.078026)

AIPW ensemble learning: Relative Risk

RR:   1.8575; 95%CI:(1.8020,1.9147)
************************************************************************

{title:Remarks} 

{p 4 4 2 120}
Remember 1: Y must be a binary numeric variable coded (0,1); X must be numeric binary
variable and, Z a vector of covariates. 
{p_end}

{p 4 4 2 120}
Remember 2: You must change your working directory to the location of the Stata.dta file.
{p_end}

{p 4 4 2 120}
Remember 3: Mac users install {hi:meltmle.ado} file. 
{p_end}

{p 4 4 2 120}
Remember 4: Windows users intall {hi:weltmle.ado} file.
{p_end}

{p 4 4 2 120}
Remember 5: Mac users must have installed R software in their personal computer as
eltmle calls R to implement the Super Learner. The R executable file must be located at 
the following path: {hi:"/usr/local/bin/r"}.
{p_end}

{p 4 4 2 120}
Remember 6: Windows users must have installed R software in their personal computer
as eltmle calls R to implement the Super Learner. The R executable file must be located 
at the following path: {hi:"C:\Program Files\R\R-3.1.2\bin\x64\R.exe"} and the version must be v.3.0.0 or later.
{p_end}

{p 4 4 2 120}
Remember 7: Windows users must have only one version of R software installed in their personal computer
at the following path: {hi:"C:\Program Files\R\R-3.1.2\bin\x64\R.exe"}. In case more than one different version 
are located in the above highlighted path users might want to keep the latest.
{p_end}

{p 4 4 2 120}
Remember 8: Check the SLS.R file for problems with R.
{p_end}

{title:References}

{p 4 4 2 120}
Luque-Fernandez, Miguel Angel. (2017). Targeted Maximum Likelihood Estimation for a
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
{phang}Department of Epidemiology and Population Health, London School of Hygiene and Tropical Medicine. 
London, UK.{p_end}
{phang}E-mail: {browse "mailto:miguel-angel.luque@lshtm.ac.uk":miguel-angel.luque@lshtm.ac.uk}{p_end}  

{title:Also see}

{psee}
Online:  {helpb teffects}
{p_end}

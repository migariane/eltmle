{smcl}
{* 30.MARCH.2017}{...}
{cmd:help eltmle}
{hline}

{title:Title}

{p2colset 5 18 20 2}{...}
{p2col:{hi:eltmle}{hline 1}}Ensemble Learning Targeted Maximum Likelihood Estimation{p2colreset}{...}

{title:Syntax}

{p 8 17 2}
{cmd:eltmle}
{hi: Y} 
{hi: X}
{hi: Z}
,
[slapiw slaipwrf slaipwbgam tmle tmlerf tmlebgam]

where:

{hi:Y}: Outcome: Numeric binary variable
{hi:X}: Treatment or exposure: Numeric binary variable
{hi:Z}: Covariates: Vector of numeric variables 

{title:Description}

{pstd} 

{hi: Modern Epidemiology} has been able to identify significant limitations of classic epidemiological methods when the focus 
is to explain the main effect of a risk factor on a disease or outcome:   

1. Non-collapsibility of the odds and hazard ratios.  
2. Impact of paradoxical effects due to conditioning on colliders.  
3. Selection bias related to the vague understanding of the effect of time on exposure and outcome.  
4. Effect of time-dependent confounding and mediators, etc.  

Classical epidemiological methods to control for confounding require making the assumption that the effect measure is constant 
across levels of confounders included in the model. Alternatively, James Robins in 1986 showed that using {hi:standardisation}
implemented through the use of the {hi: G-formula} allows obtaining a unconfounded marginal estimation of the causal average 
treatment effect (ATE) under causal assumptions.      

The most commonly used estimator for a binary treatment effect is the risk difference or ATE. The ATE estimation relies on 
parametric modelling assumptions. Therefore, the correct model specification is crucial to obtain unbiased estimates of the 
true ATE. Professor Mark van der Laan and collaborators have developed a double-robust estimation procedure to reduce bias 
against misspecification. The targeted maximum likelihood estimation (TMLE) is a semiparametric, efficient substitution 
estimator. {hi:TMLE} allows for data-adaptive estimation while obtaining valid statistical inference based on the targeted 
minimum loss-based estimation and machine learning algorithmS to minimise the risk of model misspecification.  

Evidence shows that {hi:TMLE] provides the less unbiased estimate of the ATE compared with other double robust estimators.      

The following link provides access to a TMLE tutorial: {browse "http://migariane.github.io/TMLE.nb.html":TMLE_tutorial}.        
    
{hi:eltmle} is a Stata program implementing the targeted maximum likelihood estimation for a binary outcome and treatment including 
the ensemble learning algorithm "Super Learner" called from the {hi: SuperLearner} package v.2.0-21 (Polley E., et al. 2011). 
The Super-Learner uses V-fold cross-validation (10-fold by default) to assess the performance of prediction regarding the potential 
outcomes and the propensity score as weighted averages of a set of machine learning algorithms. We used the default SuperLearner 
algorithms implemented in the base installation of the {hi:tmle-R} package v.1.2.0-5 (Susan G. and Van der Laan M., 2017), which 
included the following: i) stepwise selection, ii) generalized linear modeling (glm), iii) a glm variant that included second order 
polynomials and two-by-two interactions of the main terms included in the model. The first phase of the Super Learner algorithm is 
computationally equivalent to performing model selection via 10-folds cross-validation. The latter phase of the Super Learner 
algorithm (the meta-learning step) is just training another single model (no cross-validation) on the level one data. 
    
{title:Options}


{phang} 
{hi:tmle}: this is the default option. If no-option is specified eltmle by default implements the
TMLE based on the main three machine learning algorithms described before. 
{p_end}

{phang} 
{hi:tmlebgam}: this option may be specified or unspecified. When specified, it does include in addition to the above default
implementation for the SuperLearner call the Bayes Generalized Linear Models and the Generalized Additive Models libraries. 
{p_end}

{phang} 
{hi:slaipw}: this option may be specified or unspecified. When specified, it does estimate the augmented
inverse probability weighting algorithm plus the Super Learner ensemble learning for the main three machine 
learning algorithms described before. 
{p_end}

{phang} 
{hi:slaipwbgam}: this option may be specified or unspecified. When specified, it does include in addition to the above default
implementation for the SuperLearner call the Bayes Generalized Linear Models and the Generalized Additive Models libraries for 
the slaipw estimator. 
{p_end}


{title:Example}

*******************************************************
* eltmle Y X Z [if] [,slaipw slaipwbgam tmle tmlebgam] 
*******************************************************

.clear
.use http://www.stata-press.com/data/r14/cattaneo2.dta
.describe
.gen lbw = cond(bweight<2500,1,0.)
.gen lab var lbw “Low birthweight, <2500 g”
.save "your path/cattaneo2.dta", replace
.cd “your path”

*******************************************************
.eltmle lbw mbsmoke mage medu prenatal mmarried, tmle
*******************************************************

    Variable |        Obs        Mean    Std. Dev.       Min        Max
-------------+---------------------------------------------------------
      Q1star |      4,642    .1055429    .0403702   .0205745   .3789994
      Q0star |      4,642    .0509794    .0249112   .0207113   .1664136
          ps |      4,642    .1861267     .110755   .0372202   .8494988


TMLE: Average Treatment Effect

ATE:   0.0546; SE:0.0122; p-value: 0.0000; 95%CI:(0.030700,0.078427)

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

ATE:   0.0542; SE:0.0122; p-value: 0.0000; 95%CI:(0.030299,0.078026)

AIPW ensemble learning: Relative Risk

RR:   1.8575; 95%CI:(1.8020,1.9147)
************************************************************************

{title:Remarks} 

{pstd} 
Remember 1: Y must be a binary numeric variable coded (0,1); X must be numeric binary
variable and, Z a vector of covariates. 
{p_end}

{pstd} 
Remember 2: You must change your working directory to the location of the Stata.dta file.
{p_end}

{pstd} 
Remember 3: Mac users install {hi:meltmle.ado} file. 
{p_end}

{pstd} 
Remember 4: Windows users intall {hi:weltmle.ado} file.
{p_end}

{pstd} 
Remember 5: Mac users must have installed R software in their personal computer as 
eltmle calls R to implement the Super Learner. The R executable file must be located at 
the following path: {hi:”/usr/local/bin/r”}.
{p_end}

{pstd} 
Remember 6: Windows users must have installed R software in their personal computer 
as eltmle calls R to implement the Super Learner. The R executable file must be located 
at the following path: {hi:”C:\Program Files\R\R-3.1.2\bin\x64\R.exe”} and the version must be v.3.0.0 or later.
{p_end}

{pstd} 
Remember 7: Windows users must have only one version of R software installed in their personal computer  
at the following path: {hi:”C:\Program Files\R\R-3.1.2\bin\x64\R.exe”}. In case more than one different version 
are located in the above highlighted path users might want to keep the latest.
{p_end}

{pstd}
Remember 8: Check the SLS.R file for problems with R.
{p_end}

{title:References}

{phang}
Luque-Fernandez, Miguel Angel. (2017). Targeted Maximum Likelihood Estimation for a 
Binary Outcome: Tutorial and Guided Implementation {browse "http://migariane.github.io/TMLE.nb.html":Download here}.
{p_end}

{phang}
StataCorp. 2015. Stata Statistical Software: Release 14. College Station, TX: StataCorp LP.
{p_end}

{phang}
Gruber S, Laan M van der. (2011). Tmle: An r package for targeted maximum likelihood
estimation. UC Berkeley Division of Biostatistics Working Paper Series.
{p_end}

{phang}
Laan M van der, Rose S. (2011). Targeted learning: Causal inference for observational 
and experimental data. Springer Series in Statistics.626p.
{p_end}

{phang}
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

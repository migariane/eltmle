#eltmle: Ensemble Learning Targeted Maximum Likelihood Estimation (Stata program)  

**Modern Epidemiology** has been able to identify significant limitations of classic epidemiologic methods when the focus is to explain the main effect of a risk factor on a disease or outcome:   

1. Non-collapsibility of the odds and hazard ratios.  
2. Impact of paradoxical effects due to conditioning on colliders.  
3. Selection bias related to the vague understanding of the effect of time on exposure and outcome.  
4. Effect of time-dependent confounding and mediators, etc.  

Classical epidemiologic methods to control for confounding require making the assumption that the effect measure is constant across levels of confounders included in the model. Alternatively, James Robins in 1986 showed that using **standardisation** implemented through the use of the **G-formula** allows obtaining a unconfounded marginal estimation of the causal average treatment effect (ATE) under causal assumptions.    

The most commonly used estimator for a binary treatment effect is the risk difference or ATE. The ATE estimation relies on parametric modelling assumptions. Therefore, the **correct model specification** is crucial to obtain unbiased estimates of the true ATE.  

However, Mark van der Laan and collaborators have developed a double-robust estimation procedure to reduce bias against misspecification. The **targeted maximum likelihood estimation (TMLE)** is a semiparametric, efficient substitution estimator. **TMLE** allows for data-adaptive estimation while obtaining valid statistical inference based on the **targeted minimum loss-based estimation** and **machine learning algorithms** to minimise the risk of model misspecification.  

Evidence shows that **TMLE** provides the less unbiased estimate of the ATE compared with other double robust estimators.  
The following link provides access to a TMLE tutorial:  http://migariane.github.io/TMLE.nb.html   

**eltmle** is a Stata program implementing the targeted maximum likelihood estimation for a binary outcome and treatment including the ensemble learning algorithm "Super Learner" for three machine learning algorithms: "SL.glm" (logistic regression using A and W as covariates), "SL.step" (stepwise model selection using Akaike's Information Criterion for subsets of the full model as well as inclusion of nonlinear covariates based on second order polynomials) and, "SL.glm.interaction" (a logistic regression variant that includes both second-order polynomials and two-way interactions of the terms included in the model). A common task in machine learning is to perform model selection by specifying a number of models with different parameters. 
    The first phase of the Super Learner algorithm is computationally equivalent to performing model selection via 10-folds cross-validation. The latter
    phase of the Super Learner algorithm (the metalearning step) is just training another single model (no cross-validation) on the level one data.  

#Installation note  

NOTE: To install eltmle directly from github you need to use a Stata module for installing Stata packages from GitHub, including previous releases of a package. You can install the latest version of the github command by executing the following code in your Stata sesion:

    net install github, from("https://haghish.github.io/github/")

    then, you can install eltmle simply using the following code in Stata:

    1) For MAC users: 
    
    github install migariane/meltmle
    
    2) For Windows users:

    github install migariane/weltmle
     
    For both, MAC and Windows users, in case you want to uninstall the package type:  
	
    ado unistall eltmle  
 
#Authors  
Miguel Angel Luque-Fernandez    
Email: miguel-angel.luque@lshtm.ac.uk  
Michael Schomaker    
Email: micheal.schomaker@uct.ac.uk  
We would like to thank Karla Diaz-Ordaz and Rhian Daniel for their input and comments to improve the program.  

In case you have updates or changes that you would like to make, please send me a pull request.  
Alternatively, if you have any questions, please e-mail me.     
You can cite this repository as:  
Luque-Fernandez MA, (2017) et al. Ensemble Taregeted Maximum Likelihood Estimation for a Binary Outcome and Treatment. 
GitHub repository, https://github.com/migariane/STATA-TMLE      
Miguel Angel Luque-Fernandez    
E-mail: miguel-angel.luque at lshtm.ac.uk  
Twitter @WATZILEI  

#Copyright
This software is distributed under the GPL-2 license.



#eltmle: Ensemble Learning Targeted Maximum Likelihood Estimation (Implementation for Stata software)  

**Modern Epidemiology** has been able to identify significant limitations of classic epidemiological methods when the focus is to explain the main effect of a risk factor on a disease or outcome:   

1. Non-collapsibility of the odds and hazard ratios.  
2. Impact of paradoxical effects due to conditioning on colliders.  
3. Selection bias related to the vague understanding of the effect of time on exposure and outcome.  
4. Effect of time-dependent confounding and mediators, etc.  

Classical epidemiological methods to control for confounding require making the assumption that the effect measure is constant across levels of confounders included in the model. Alternatively, James Robins in 1986 showed that using **standardisation** implemented through the use of the **G-formula** allows obtaining a unconfounded marginal estimation of the causal average treatment effect (ATE) under causal assumptions.    

The most commonly used estimator for a binary treatment effect is the risk difference or ATE. The ATE estimation relies on parametric modelling assumptions. Therefore, the **correct model specification** is crucial to obtain unbiased estimates of the true ATE. Professor Mark van der Laan and collaborators have developed a double-robust estimation procedure to reduce bias against misspecification. The **targeted maximum likelihood estimation (TMLE)** is a semiparametric, efficient substitution estimator. **TMLE** allows for data-adaptive estimation while obtaining valid statistical inference based on the **targeted minimum loss-based estimation** and **machine learning algorithms** to minimise the risk of model misspecification. Evidence shows that **TMLE** provides the less unbiased estimate of the ATE compared with other double robust estimators.      

The following link provides access to a TMLE tutorial:  http://migariane.github.io/TMLE.nb.html     

**eltmle** is a Stata program implementing the targeted maximum likelihood estimation for a binary outcome and treatment including the ensemble learning algorithm "Super Learner" called from the **SuperLearner** package v.2.0-21 (Polley E., et al. 2011). The Super-Learner uses V-fold cross-validation (10-fold by default) to assess the performance of prediction regarding the potential outcomes and the propensity score as weighted averages of a set of machine learning algorithms. We used the default SuperLearner algorithms implemented in the base installation of the **tmle-R** package v.1.2.0-5 (Susan G. and Van der Laan M., 2017), which included the following: i) stepwise selection, ii) generalized linear modeling (glm), iii) a glm variant that included second order polynomials and two-by-two interactions of the main terms included in the model. The first phase of the Super Learner algorithm is computationally equivalent to performing model selection via 10-folds cross-validation. The latter phase of the Super Learner algorithm (the meta-learning step) is just training another single model (no cross-validation) on the level one data.    

#Installation note    

To install eltmle directly from github you need to use a Stata module for installing Stata packages from GitHub, including previous releases of a package. You can install the latest version of the github command by executing the following code in your Stata session:

    net install github, from("https://haghish.github.io/github/")

    then, you can install eltmle simply using the following code in Stata:

    1) For MAC users (please, read carefully the help file before using eltmle in Stata):  
    
    github install migariane/meltmle  
   
    which eltmle   

    help eltmle   

    2) For Windows users (please, read carefully the help file before using eltmle in Stata):  

    github install migariane/weltmle  
 
    which eltmle 

    help eltmle   
     
    For both, MAC and Windows users, in case you want to uninstall the package type:    
	
    ado unistall eltmle   
     
 
#Author 
[Author and Developer]
Miguel Angel Luque-Fernandez, LSHTM, NCDE, Cancer Survival Group, London, UK    
Email: miguel-angel.luque@lshtm.ac.uk 

#Contributors:
[Intellectual advice and suggestions for improvement]

Michael Schomaker, CIDER, UCT, Cape Twon, South Africa      
Email: michael.schomaker at uct.ac.za    

Bernard Rachet, LSHTM, NCDE, Cancer Survival Group, London, UK  
Email: bernard.rachet at lsthm.ac.uk  

Aurelien Belot, LSHTM, NCDE, Cancer Survival Group, London, UK  
Email: Aurelien.belot at lshtm.ac.uk  

Camille Maringe, LSHTM, NCDE, Cancer Survival Group, London, UK  
Email: camille.maringe at lshtm.ac.uk  

Karla Diaz Ordaz, LSHTM, Department of Medical Statistics, London, UK    
Email: Karla.diaz-ordaz at lshtm.ac.uk    

Rhian Daniel, LSHTM, Department of Medical Statistics, London, UK    
Email: Rhian.daniel at lshtm.ac.uk    

#Acknowledgements  
  
I would like to thank Professors Bianca De Stavola (LSHTM), Simon Cousens (LSHTM), Aurelio Tobias (CSIC), Michel Coleman (LSHTM) for their comments and support and Haghish E. F. (CMBMI, Freiburg, Germany) for his wonderful **Github** and **MarkDoc** Stata packages. 
  
In case you have updates or changes that you would like to make, please send me a pull request.  
Alternatively, if you have any questions, please e-mail me.     
You can cite this repository as:  
Luque-Fernandez MA, (2017). Ensemble Taregeted Maximum Likelihood Estimation for a Binary Outcome and Treatment. 
GitHub repository, https://github.com/migariane/weltmle (Windows users) or https://github.com/migariane/meltmle (MAC users)        
Miguel Angel Luque-Fernandez    
E-mail: miguel-angel.luque at lshtm.ac.uk  
Twitter @WATZILEI  

#Copyright

This software is distributed under the GPL-2 license.

![Figure Link](https://github.com/migariane/meltmle/blob/master/acknowledgement.png)   

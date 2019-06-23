# eltmle: Ensemble Learning Targeted Maximum Likelihood Estimation (Implementation for Stata software)  

**Modern Epidemiology** has been able to identify significant limitations of classic epidemiological methods, like outcome regression analysis, when estimating causal quantities such as the average treatment effect (ATE) or the causal odds ratio, for observational data.       

For example, using classical regression models to estimate the ATE requires making the assumption that the effect measure is constant across levels of confounders included in the model, i.e. that there is no effect modification. Other methods do not require this assumption, including g-methods (e.g. the **g-formula**) and targeted maximum likelihood estimation (**TMLE**).     

The latter estimator has the advantage of being **doubly robust**. Moreover, TMLE allows inclusion of **machine learning** algorithms to minimise the risk of model misspecification, a problem that persists for competing estimators. Evidence shows that TMLE typically provides the **least unbiased** estimates of the ATE compared with other double robust estimators.         

The following link provides access to a TMLE tutorial:  http://migariane.github.io/TMLE.nb.html       

**eltmle** is a Stata program implementing the targeted maximum likelihood estimation for the ATE for a binary or continuous outcome and binary treatment. Future implementations will offer more general settings. **eltmle** includes the use of a "Super Learner" called from the **SuperLearner** package v.2.0-21 (Polley E., et al. 2011). The Super-Learner uses V-fold cross-validation (10-fold by default) to assess the performance of prediction regarding the potential outcomes and the propensity score as weighted averages of a set of machine learning algorithms. We used the default SuperLearner algorithms implemented in the base installation of the **tmle-R** package v.1.2.0-5 (Susan G. and Van der Laan M., 2017), which included the following: i) stepwise selection, ii) generalized linear modeling (glm), iii) a glm variant that included second order polynomials and two-by-two interactions of the main terms included in the model.    

# Installation note    

To install eltmle directly from github you need to use a Stata module for installing Stata packages from GitHub, including previous releases of a package. You can install the latest version of the github command by executing the following code in your Stata session.    
**Note**: you need a Stata version greater or equal than 13.2, otherwise you can install the package manually downloading the files from the Github repository and placing it in your Stata ADO PERSONAL folder:   

    net install github, from("https://haghish.github.io/github/")

    then, you can install eltmle simply using the following code in Stata:

    1) Please, read carefully the help file before using eltmle in Stata):  
    
    github install migariane/eltmle  
   
        .which eltmle   

        .help eltmle   
   
    For both, MAC and Windows users, in case you want to uninstall the package type:    
	
        .ado uninstall eltmle   
     
 
# Author 
[Author and Developer]
Miguel Angel Luque-Fernandez, LSHTM, NCDE, Cancer Survival Group, London, UK    
Email: miguel-angel.luque@lshtm.ac.uk 

# Contributors:
[Intellectual advice and suggestions for improvement]

Michael Schomaker, CIDER, UCT, Cape Town, South Africa      
Email: michael.schomaker at uct.ac.za 

Mireille E. Schnitzer, Faculté de pharmacie,     
Université de Montréal, Montréal, Canada  
Email: mireille.schnitzer at umontreal.ca   
 
# Updates
In case you have updates or changes that you would like to make, please send me a pull request.  
Alternatively, if you have any questions, please e-mail me:   
Miguel Angel Luque-Fernandez    
E-mail: miguel-angel.luque at lshtm.ac.uk  
Twitter @WATZILEI  

# Citation    
You can cite this repository as:  
Luque-Fernandez MA, (2018). Ensemble Targeted Maximum Likelihood Estimation for a Binary Treatment. 
GitHub repository, https://github.com/migariane/eltmle          

[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.2560828.svg)](https://doi.org/10.5281/zenodo.2560828)

# Copyright
This software is distributed under the GPL-2 license.

# Acknowledgments  
Miguel Angel Luque Fernandez is supported by a Miguel Servet I Investigator Award (grant CP17/00206) from the Carlos III National Institute of Health, Madrid, Spain.  
  
I would like to thank Professors Bianca De Stavola (LSHTM), Simon Cousens (LSHTM), Aurelio Tobias (CSIC), Michel Coleman (LSHTM) for their comments and support and Haghish E. F. (CMBMI, Freiburg, Germany) for his wonderful **Github** and **MarkDoc** Stata packages. 

![Figure Link](https://github.com/migariane/eltmle/blob/master/Acknowledgment.png)   

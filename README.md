# Bayesian_estimation_of_Pa_viscoelasticity
Bayesian estimation of Pseudomonas aeruginosa viscoelastic properties based on creep of wild type, rugose, and mucoid variant biofilms

Here, we provide Matlab files to estimate the viscoelastic properties of P. aeruginosa. Please refer to the manuscript entitled "Bayesian estimation of Pseudomonas aeruginosa viscoelastic properties based on creep of wild type, rugose, and mucoid variant biofilms" for more details on the mathematical model and discussions. The Matlab file 'main.m' imports observational data, initiate the model parameters, and calls the forward model and the 'mcmcstat' package. The function 'forward_burger.m' provides the mathematical fowrward problem to calculate biofilm strain response during a creep-recovery test. The function 'likelihood.m' computes likelihood, the probability of the observed data creep-recovery, given the model parameters and the error variance.

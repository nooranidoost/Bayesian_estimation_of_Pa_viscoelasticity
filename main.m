clear all;
clc;

% This code estimates viscoelastic properties of P. aeruginosa biofilms
% during a creep-recovery test base on observational data. We described 
% strain response of Wild Type P. aeruginosa, and isogenic RSCV and mucoid
% variant biofilms using a Burger viscoelastic model. More details about
% this analysis can be found in:
% M. Nooranidoost, et al., 2023. Bayesian estimation of Pseudomonas
% aeruginosa viscoelastic properties based on creep of wild type, rugose,
% and mucoid variant biofilms, Biofilm (Under review)    
%
% This code uses MCMCSTAT package to generate MCMC chains using
% multivariate Gausssian proposal distribution. Please refer to the
% following paper for more detailed information about the MCMC algoirthm.
% H. Haario, et al., 2006. DRAM: Efficient adaptive MCMC, Statistics and
% Computing 16, pp. 339-354.


% import observation data
strain_obs1 = importdata('strain_observation1.txt');
strain_obs2 = importdata('strain_observation1.txt');
strain_obs3 = importdata('strain_observation1.txt');
strain_obs4 = importdata('strain_observation1.txt');
time_obs    = importdata('time_observation');

data.xdata = time_obs;
data.ydata = (strain_obs1 + strain_obs2 + strain_obs3 + strain_obs4);

% set initial guess for viscoelastic parameters
E_k_init   = 50;
eta_k_init = 50;
E_m_init   = 50;
eta_m_init = 50;

theta = [E_k_init eta_k_init E_m_init eta_m_init];

% set model parameters first guess and prior
params = {
    {'E_k',    theta(1), 0, inf}
    {'\eta_k', theta(2), 0, inf}
    {'E_m',    theta(3), 0, inf}
    {'\eta_m', theta(4), 0, inf}
    };

% set mathematical model
ssfun         = @likelihood;
model.ssfun   = ssfun;

sd            = std([strain_obs1 strain_obs2 strain_obs3 strain_obs4]);
model.sigma2  = sd;
model.N       = length(time_obs);

% set MCMC algorithm properties
options.method      = 'dram'; 
options.updatesigma = 1;
options.nsimu       = 1000; 

% generate the inital chain (burn-in)
[results, chain, s2chain] = mcmcrun(model,data,params,options);

options.nsimu = 1000; 

% re-run and generate the main chain using the previous run results as the initial condition 		      
[results, chain, s2chain] = mcmcrun(model,data,params,options,results);

% save statistics 
save    ('results.mat','-struct','results');
dlmwrite('chain.txt'  ,chain  ,'delimiter','\t','precision',5);
dlmwrite('s2chain.txt',s2chain,'delimiter','\t','precision',5);



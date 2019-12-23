%% Local Approximations for Subset Simulation method
%{
---------------------------------------------------------------------------
Created by:
Kenan Šehić (kense@dtu.dk; kenosehic@gmail.com)
Department of Applied Mathematics and Computer Science
Technical University of Denmark
Licence: Copyright (C) 2019 Kenan Šehić DTU Compute, Technical University of Denmark

Cite: Šehić K., Karamehmedović M., Marzouk Y.: Estimation of Failure Probabilities via Local Subset Approximations, TBD
---------------------------------------------------------------------------
Version December 2019
---------------------------------------------------------------------------
Description:
* Space: Standard Normal N(0,1)
* MCMC: Component-wise Metropolis Algorithm
* Regresion: Gaussian Process
* Dimensionality Reduction: Partial Least Square regression (PLS1)
--------------------------------------------------------------------------
%}
%% Test local approximation approach

seed=144; % seed

N = 5000; % initial number of samples at j=0

p0 = 0.1; % conditional probability which defines failure threshold

g = 8; % select function to test local approximation approach

d = 300; % dimension

off = 1;
   
[Pf_SuS,delta_SuS,prob,b,Pf,b_line,Pf_line,samplesU,gruns] = SuS_Local(N, p0, g,d,off,seed);
    
filename='test.mat';
 
save(filename)

%% THE END
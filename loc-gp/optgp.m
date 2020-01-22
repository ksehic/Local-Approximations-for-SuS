function Q = optgp(gpmodel,X,b,Rpls11,Xmean)
%% Local Learning based on Gaussian process for Subset Simulation method
%{
---------------------------------------------------------------------------
Created by:
Kenan Šehić (kense@dtu.dk; kenosehic@gmail.com)
Department of Applied Mathematics and Computer Science
Technical University of Denmark
Licence: Copyright (C) 2019 Kenan Šehić DTU Compute, Technical University of Denmark

Cite: Šehić K., Karamehmedović M.: Estimation of Failure Probabilities via Local Subset Approximations, TBD
---------------------------------------------------------------------------
Version December 2019
---------------------------------------------------------------------------
Description:
*   Find optimal sample that minimize U-function for trained Gaussian process
    and specific failure threshold
---------------------------------------------------------------------------
* Main Ref:
      
    R. Schöbi, B. Sudret, S. Marelli, Rare event estimation using polynomial-chaos kriging, 
    ASCE-ASME Journal of Risk and Uncertainty in Engineering Systems, Part A: Civ. Eng. 3 (2).

%}
%% Procedure

X1 = X - Xmean; % centering

[Z,Zsd]= predict(gpmodel,X1*Rpls11); % predict with trained GP

Q = abs(Z+b)./Zsd; % estimate U-value

return

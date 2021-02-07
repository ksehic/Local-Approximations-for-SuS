function [x,Yadd,run] = findpoint(gpmodel,g,b,R,Rpls11,Xmean,theta)
%% Local Learning based on Gaussian process for Subset Simulation method
%{
---------------------------------------------------------------------------
Created by:
Kenan Šehić (kense@dtu.dk; mrsehickenan@gmail.com)
Department of Applied Mathematics and Computer Science
Technical University of Denmark
Licence: Copyright (C) 2019 Kenan Šehić DTU Compute, Technical University of Denmark

Cite: Šehić K., Karamehmedović M.: Estimation of Failure Probabilities via Local Subset Approximations, TBD
---------------------------------------------------------------------------
Version December 2019
---------------------------------------------------------------------------
Description:

* Model Refinement to find the bext next sample within the ball to improve predictions

%}

%% Procedure

options = optimset('Largescale','off','Display','off');

lb = [];
ub = [];

x0 = theta;

[x,~] = fmincon(@(x) optgp(gpmodel,x,b,Rpls11,Xmean),x0.',[],[],[],[],lb,ub,@(x)cons(x,x0.',R),options); % search for the best sample

x=x.';

Yadd=g(x);

run=1;

return

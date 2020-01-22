function [c,ceq] = cons(theta0,seed,R)
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

* Constraint part for Local Model Refinement to select a sample within the ball

%}
%% Procedure

c1=norm(theta0-seed)-R; % Radius of the ball

c=c1;

ceq=[];

return

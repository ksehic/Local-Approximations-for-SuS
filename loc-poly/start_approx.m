function [lsf1,eplus] = start_approx(theta,seeds,g_limit,skipsample,g)
%% Local Polynomial Approximations for Subset Simulation method
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
* Fit locally polynomial regression at j=0
---------------------------------------------------------------------------
* Main Ref:
      
    P. R. Conrad, Y. M. Marzouk, N. S. Pillai and A. Smith: Accelerating Asymptotically Exact MCMC for Computationally Intensive Models via Local Approximations.
    Journal of the American Statistical Association, Volume 111, 2016 - Issue 516, Pages 1591-1607

%}
%% Procedure
if isempty(skipsample)~=1
    for ii=1:length(skipsample(1,:))

        [~,cc]=find(seeds==skipsample(:,ii));

        seeds(:,cc(1)) = [];

        g_limit(cc(1)) = [];

    end
end
[seeds,idseed,~]=unique(seeds.','row','stable'); % remove the same points

seeds=seeds.';

g_limit=g_limit(idseed);
%% collecting samples
d=length(theta);%d=2; % dimensional of the problem

Ndef=(d+1)*(d+2)/2; % how many points include in Rdef circle

N=ceil(sqrt(d)*Ndef); % how many points include total in R circle

%% Least Square quadratic regression
[WZP,W,Fi,id,~,~] = WZP0(theta,seeds,N,Ndef);

gsort=g_limit(id(1:N));

gsort=gsort.';

rhs = Fi.'*W*gsort;

Z = WZP\rhs;

lsf1=Z(1);

[~,~,seed0,~,~] = RmaxRefine(theta,seeds,g_limit,[]);

eplus=start_errors(theta,seeds,g_limit,seed0,lsf1,g);

return

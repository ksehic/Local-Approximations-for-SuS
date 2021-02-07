function [Rmax,Rdef,seeds0,eval0,seeds1] = RmaxRefine(theta,seeds,g_limit,skipsample)
%% Local Polynomial Approximations for Subset Simulation method
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
* Select Radius and samples within the ball
---------------------------------------------------------------------------
* Main Ref:
      
    P. R. Conrad, Y. M. Marzouk, N. S. Pillai and A. Smith: Accelerating Asymptotically Exact MCMC for Computationally Intensive Models via Local Approximations.
    Journal of the American Statistical Association, Volume 111, 2016 - Issue 516, Pages 1591-1607

%}
%% Procedure

if isempty(skipsample)~=1 % which sample to neglect

    [~,cc]=find(seeds==skipsample);

    seeds(:,cc(1)) = [];

end

[seeds,~,~]=unique(seeds.','row','stable'); % remove the same points

seeds=seeds.';

%% collecting samples
d=length(theta);%d=2; % dimensional of the problem

Ndef=(d+1)*(d+2)/2; % how many points include in Rdef circle

N=ceil(sqrt(d)*Ndef); % how many points include total in R circle

R=zeros(length(seeds),1);

for i=1:length(seeds)
    
    R(i)=norm(seeds(:,i)-theta);   
    
end

[Rsort,id]=sort(R);

Rmax = max(Rsort(1:N));

Rdef = max(Rsort(1:Ndef));

seeds0 = seeds(:,id(1:N));

eval0 = g_limit(1,id(1:N));

seeds1 = seeds(:,id(1:Ndef));

return

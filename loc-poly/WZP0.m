function [WZP,W,Fi,id,Rmax,Rdef] = WZP0(theta,seeds,N,Ndef)
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
* Prepare samples as a part of quadratic polynomial regression within AX=B
---------------------------------------------------------------------------
* Main Ref:
      
    P. R. Conrad, Y. M. Marzouk, N. S. Pillai and A. Smith: Accelerating Asymptotically Exact MCMC for Computationally Intensive Models via Local Approximations.
    Journal of the American Statistical Association, Volume 111, 2016 - Issue 516, Pages 1591-1607

%}
%% Procedure

R=zeros(length(seeds)-1,1);

for i=1:length(seeds)-1
    
    R(i)=norm(seeds(:,i)-theta);   
    
end

[Rsort,id]=sort(R);

Rmax = max(Rsort(1:N));

Rdef = max(Rsort(1:Ndef));

w = zeros(1,N); % weights for samples within the ball

for i=1:N
   
    wcheck = norm(seeds(:,id(i))-theta);
    
    if wcheck <= Rdef
        
        w(i) = 1;
        
    elseif wcheck > Rmax
        
        w(i) = 0;
        
    else
        
        w(i) = (1-((norm(seeds(:,id(i))-theta)-Rdef)/(Rmax-Rdef))^3)^3;
    
    end
    
    
end

W=diag(w);

seeds_scale = (seeds(:,id(1:N))-theta)/Rmax;

Fi=Fin(N,seeds_scale); % multidimensional matrix A for regression

WZP = Fi.'*W*Fi;

return
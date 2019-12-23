function [Z,ZV,ZUP,R,Rpls11,XmeanLoc,gpmodel1] = locapprox(theta,seeds,g_limit,skipsample,g)
%% Local Learning based on Gaussian process for Subset Simulation method
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
* Local approximation pre-processing
---------------------------------------------------------------------------
Input:
* theta         : MCMC proposal
* seeds         : Samples
* g_limit       : Evalutions
* skipsample    : sample that is omitted
* g             : Numerical model
---------------------------------------------------------------------------
Output:
* Z             : Prediction
* ZV            : Standard deviation
* ZUP           : Confidence interval
* R             : Radius of ball
* Rpls11        : Low-dimensional projection matrix
* XmeanLoc      : Mean value for initial samples within ball
* gpmodel1      : trained Gaussian process
---------------------------------------------------------------------------
%}
%% Procedure to prepare data for local approximations

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

gsort = g_limit;

gsort=gsort.';

%% Run local approximations with Gaussian process for MCMC proposal

[Z,ZV,ZUP,R,Rpls11,XmeanLoc,gpmodel1]  = localpls_update(seeds.',gsort,theta.',g);

return

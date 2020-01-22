function eplus=start_errors(thetaplus,seeds,gsort,seeds1,tryrun1,g)
%% Local Polynomial Approximations for Subset Simulation method
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
* Leave-one-out error estimation at j=0
---------------------------------------------------------------------------
* Main Ref:
      
    P. R. Conrad, Y. M. Marzouk, N. S. Pillai and A. Smith: Accelerating Asymptotically Exact MCMC for Computationally Intensive Models via Local Approximations.
    Journal of the American Statistical Association, Volume 111, 2016 - Issue 516, Pages 1591-1607

%}
%% Procedure

        tryrunplus = zeros(length(seeds1),1);
        
        for i=1:length(seeds1) % potential parallel computing
            
           tryrunplus(i) = locapprox(thetaplus,seeds,gsort,seeds1(i),g);
        
        end
        
        eplus = max(abs(tryrun1-tryrunplus));

return

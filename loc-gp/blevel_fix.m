function [samples_sort,geval_sort,gruns,newsamples,newruns] = blevel_fix(samples,geval,g,N,p0,ZUP,newsamples,newruns)
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
* Intermediate failure threshold fix
---------------------------------------------------------------------------
Input:
* samples       : Current samples
* geval         : Current evalutions
* g             : Numerical model
* N             : Initial number of samples
* p0            : Conditional probability of each subset
* ZUP           : Diff for the confidence interval of GP
* newsamples    : New samples
* newruns       : New evalutions
---------------------------------------------------------------------------
Output:
* samples_sort      : Sorted samples
* geval_sort        : Sorted evalutions
* gruns             : Number of evalutions
* newsamples        : New samples
* newruns           : New evalutions
---------------------------------------------------------------------------
* Main Ref:
      
    J. Li, D. Xiu, Evaluation of failure probability via surrogate models, Journalof Computational Physics 229 (23) (2010) 8966–8980.

%}

%% Procedure to correct intermediate failure thresholds

[geval_sort, idx] = sort(geval);
samples_sort= samples(:,idx);   % order the samples by idx

[geval_unique,ia,ic] = unique(geval_sort,'last');

seeds_unique = samples_sort(:,ia);

gruns = 0;

Nuq = length(ia);

b0 = prctile(geval_sort, p0*100); % estimate initial failure threshold

dM = round(1e-1*Nuq); % sample set for fix

dM0 = sum(geval_unique<=b0+0.1*b0); % upper limit which samples check
dM1 = sum(geval_unique<=b0+0.1*b0)-dM; % lower limit which samples check

low_bc = sum(geval_unique<=b0-0.1*b0); % correct failure threshold within 10%

%% Iteration

while true
    
   arM = dM1:dM0; % range which samples check
   
   if dM1 <= low_bc
       
       arM = low_bc:dM0;
       
   end
   
%% Use numerical model   
    for i=1:length(arM) % run numerical model and replace evalutions
        
        j=arM(i);
        
        runn=g(seeds_unique(:,j));
        
        if runn==geval_unique(j)
        
            gruns = gruns + 0;
            
        else
            
            gruns = gruns + 1;
            
        end
    
        geval_unique(j) = runn;
               
        [~,new_c] = find(newsamples==seeds_unique(:,j));
        
        newruns(new_c) = geval_unique(j);
    
    end
    
    gsort_new = geval_unique(ic);
    
%% Update failure threshold
    
    b1 = prctile(gsort_new,p0*100); % estimate updated failure threshold

    if abs(b1 - b0) <= 1e-3 || arM(1)==low_bc % check relative change within failure threshold
        
        break;
        
    end
    
%% Next iteration
    
    dM0 = dM1 - 1;
    
    dM1 = dM1 - dM;
    
    b0=b1;
    
end

geval_sort = geval_unique(ic); % collect samples

samples_sort = seeds_unique(:,ic); % collect evalutions

return
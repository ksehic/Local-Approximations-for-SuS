function [Pf,runs] = lastfix(seeds,geval,b0,N,g)
%% Local Approximations for Subset Simulation method
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
* Final failure threshold fix
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

%% Procedure to correct final failure thresholds

[gsort, idx] = sort(geval);
seedssort= seeds(:,idx);   % order the samples by idx

[geval_unique,ia,ic] = unique(gsort,'last');

seeds_unique = seedssort(:,ia);

runs = 0;

Nuq = length(ia);

Pf = sum((gsort<=b0))/N; % estimate initial failure probability

dM = round(1e-1*Nuq); % define 10% range

dM0 = 1;
dM1 = dM;

if isempty(Nuq)==1 || isnan(Nuq)==1
   
    return
    
end

%% Iteration
while true
      
   arM = dM0:dM1;
   
   if dM1 > Nuq
       
       arM = dM0:Nuq; % define which samples
       
       dM1 = Nuq;
       
   end

%% Use numerical model

    for i=1:length(arM)
        
        j=arM(i);
    
        runn = g(seeds_unique(:,j));
        
        [~,crun]=find(gsort==runn);
        
        if isempty(crun)==1
            
            runs = runs + 0;
            
        else
            
            runs = runs + 1;
            
        end
        
        geval_unique0(i) = g(seeds_unique(:,j)); % run numerical model
    
    end
    
    if isempty(geval_unique0)==1 || isempty(dM1)==1 || isnan(dM1)==1 || sum(isnan(geval_unique0))==1
        
        break
        
    end
    
    [cc0,~] = find(ic==dM0);
    
    [cc1,~] = find(ic==dM1);

%% Correct failure probability
    
    grun = geval_unique0(ic(cc0(1):cc1(end))-ic(cc0(1))+1);
   
    Pf1 = Pf + sum(-(gsort(1,cc0(1):cc1(end))<=b0)+(grun<=b0))/N; % correct failure probability

%% Check improvement

    if Pf1<=0 || sum(isnan(Pf1))==1 || isempty(Pf1)==1 % check improvement
        
        break
        
    end
    
    if abs(Pf1 - Pf) <= 1e-3 || dM1==Nuq  % check improvement
        
        Pf=Pf1;
        
        break;
        
    end
    
%% Next iteration
    
    dM0 = dM1 + 1;
    
    dM1 = dM1 + dM;
    
    Pf=Pf1;
    
    geval_unique0 = [];
    
    grun = [];
    
end


return

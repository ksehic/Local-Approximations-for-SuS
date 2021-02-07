function [Ge,seeds,gseeds,run1,check,epadd] = fulllocal_v2(v,x0,seeds,gseeds,G_LSF,t,j,run1,epadd)
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
* Local Polynomial approximations based on qudratic polynomial fit
* Leave-one-out error estimations
* N=d^0.5*(d+1)*(d+2)/2
* Random and error refinement
* Based on L2, optimal sample is selected within the ball
---------------------------------------------------------------------------
* Main Ref:
      
    P. R. Conrad, Y. M. Marzouk, N. S. Pillai and A. Smith: Accelerating Asymptotically Exact MCMC for Computationally Intensive Models via Local Approximations.
    Journal of the American Statistical Association, Volume 111, 2016 - Issue 516, Pages 1591-1607

%}
%% Procedure

      bs=@(b0,b1,kk,ll) b0*kk.^(-b1*ll^2); % random refinement
      
      [Rmax,Rdef,seed0,~,seed1] = RmaxRefine(v,seeds,gseeds,[]); % radius R and samples at a candidate v
      
      [RmaxEX,RdefEX,seed0EX,~,seed1EX] = RmaxRefine(x0,seeds,gseeds,[]);  % radius R and samples at a candidate x0
      
      %% Approximate at v and x0
      
      [Ge,run] = locapprox(v,seeds,gseeds,[],G_LSF);
          
      run1 = run1  + run;
      
      [GeEX,run] = locapprox(x0,seeds,gseeds,[],G_LSF);
      
      run1 = run1  + run;
      %% Leave-one out error estimations

      [eplus,eminus]=errors(v,x0,seeds,gseeds,seed0EX,GeEX,seed0,Ge,G_LSF);
      
      %% Refinment
      if eplus>0.5 || eplus==epadd
          
              Ge=G_LSF(v);

              check=0;
              
              run1=run1 + 1;
          
      elseif bs(1,0.01,t,j)<rand(1,1) % random refiment
          
              if rand(1,1)<rand(1,1)
              
                  newseed = refinenear(x0,seed0EX,RmaxEX*(1+j)/2); % refine

                  neweval = G_LSF(newseed); % running model for a new point

                  seeds(:,end+1) = newseed; % adding new point to collection

                  gseeds(:,end+1) = neweval; % adding new eval to collection

                  run1=run1 + 1;

                  check=1;
              
              else
                  
                  newseed = refinenear(v,seed0,Rmax*(1+j)/2); % refine

                  neweval = G_LSF(newseed); % running model for a new point

                  seeds(:,end+1) = newseed; % adding new point to collection

                  gseeds(:,end+1) = neweval; % adding new eval to collection

                  run1=run1 + 1;

                  check=1;
              
              end
              %fprintf('Random\n');
                        
      elseif eplus>=eminus && eplus >= 0.1 %
          
              newseed = refinenear(v,seed1,Rdef); % refine

              neweval = G_LSF(newseed); % running model for a new point

              seeds(:,end+1) = newseed; % adding new point to collection

              gseeds(:,end+1) = neweval; % adding new eval to collection

              run1=run1 + 1;
              
              check=1;
              
              %fprintf('Eplus %g \n',eplus);
              
      elseif eplus<eminus && eminus >= 0.1 %
          
              newseed = refinenear(x0,seed1EX,RdefEX); % refine

              neweval = G_LSF(newseed); % running model for a new point

              seeds(:,end+1) = newseed; % adding new point to collection

              gseeds(:,end+1) = neweval; % adding new eval to collection

              run1=run1 + 1;
              
              check=1;
              
              %fprintf('Eminus %g \n',eplus);
          
      else
          
          check=0;
          
      end
      

      epadd=eplus;
      
return


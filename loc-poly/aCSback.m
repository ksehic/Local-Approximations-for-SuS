function [u_jp1,geval,new_lambda,sigma,accrate,run1,err1,seeds,gseeds] = aCSback(N, old_lambda, b, bp, u_j,g_j, G_LSF,totalu,totalg,j,off)
%% Local Approximations for Subset Simulation method
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
* Space: Standard Normal N(0,1)
* MCMC: Component-wise Metropolis Algorithm
* Regresion: Gaussian Process
* Dimensionality Reduction: Partial Least Square regression (PLS1)
---------------------------------------------------------------------------
Input:
* N          : number of samples to be generated
* old_lambda : scaling parameter lambda
* b          : actual intermediate level
* u_j        : seeds used to generate the new samples
* G_LSF      : numerical model
* totalu     : sample data
* totalg     : evalution data
* j          : failure level
* off        : local approximation turn on/off
---------------------------------------------------------------------------
Output:
* u_jp1      : next level samples
* geval      : limit state function evaluations of the new samples
* new_lambda : next scaling parameter lambda
* sigma      : standard deviation of the proposal
* accrate    : acceptance rate
* run1       : number of evalutions
* err1       : average relative error
* seeds      : sample data
* gseeds     : evalution data
* geval_up   : Gaussian 95%CI
---------------------------------------------------------------------------
* Main Ref:
      
    Subset Simulation Method - Felipe Uribe (felipe.uribe@tum.de)
    https://www.bgu.tum.de/era/software/software00/subset-simulation/
    
    Papaioannou I., Betz W., Zwirglmaier K., Straub D.: MCMC algorithms for subset simulation. Probabilistic Engineering Mechanics, 41: 89-103.

    P. R. Conrad, Y. M. Marzouk,N. S. Pillai and A. Smith: Accelerating Asymptotically Exact MCMC for Computationally Intensive Models via Local Approximations.
    Journal of the American Statistical Association, Volume 111, 2016 - Issue 516, Pages 1591-1607

%}

%% MCMC setup: 0. initialize variables Ref: https://www.bgu.tum.de/era/software/software00/subset-simulation/
n  = size(u_j,1);      % number of uncertain parameters (dimension)
Ns = size(u_j,2);      % number of seeds
Na = ceil(100*Ns/N);   % number of chains after which the proposal is adapted

% number of samples per chain
Nchain = ones(Ns,1)*floor(N/Ns);
Nchain(1:mod(N,Ns)) = Nchain(1:mod(N,Ns))+1;

% initialization
u_jp1  = zeros(n,N);                 % generated samples
geval  = zeros(1,N);                 % store lsf evaluations
acc    = zeros(1,N);                 % store acceptance
mu_acc = zeros(1,floor(Ns/Na)+1);    % store mean acceptance until adaptation
hat_a  = zeros(1,floor(Ns/Na));      % average acceptance rate of the chains
lambda = zeros(1,floor(Ns/Na)+1);    % scaling parameter \in (0,1)

%% MCMC setup: 1. compute the standard deviation Ref: https://www.bgu.tum.de/era/software/software00/subset-simulation/
opc = 'a';
switch opc
   case 'a'   % 1a. sigma = ones(n,1);
      sigma_0 = ones(n,1);
      
   case 'b'   % 1b. sigma = sigma_hat; (sample standard deviations)
      mu_hat  = mean(u_j,2);    % sample mean
      var_hat = zeros(n,1);     % sample std
      for i = 1:n   % dimensions
         for k = 1:Ns   % samples
            var_hat(i) = var_hat(i) + (u_j(i,k)-mu_hat(i))^2;
         end
         var_hat(i) = var_hat(i)/(Ns-1);
      end
      sigma_0 = sqrt(var_hat);
   otherwise
      error('Choose a or b');
end

%% MCMC setup: 2. iteration Ref: https://www.bgu.tum.de/era/software/software00/subset-simulation/
star_a    = 0.44;         % optimal acceptance rate
lambda(1) = old_lambda;   % initial scaling parameter \in (0,1)

% a. compute correlation parameter
i         = 1;                                  % index for adaptation of lambda
sigma     = min(lambda(i)*sigma_0,ones(n,1));   % Ref. 1 Eq. 23
rho       = sqrt(1-sigma.^2);                   % Ref. 1 Eq. 24
mu_acc(i) = 0;

%% Local Approximation within MCMC

run1 = 0; % number of extra evalutions
err1 = 0; % relative error estimation

seeds = totalu; % collect data
gseeds = totalg;

bs=@(b0,b1,kk,ll) b0*kk.^(-b1*ll^2); % random refinement function

kkk=1;

%% MCMC steps
for k = 1:Ns
    
   idx          = sum(Nchain(1:k-1))+1;  %((k-1)/pa+1);
   
   u_jp1(:,idx) = u_j(:,k);              % pick a seed at random
   
   geval(idx) = g_j(k);                  % pick corresponding evalution
   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% MCMC proposal
   for t = 1:Nchain(k)-1
      % generate candidate sample
      v = normrnd(rho.*u_jp1(:,idx+t-1), sigma);
      x0 = u_jp1(:,idx+t-1); % previous step
      
%% Use numerical model to check correct estimation. This is only for testing. It should be remove for serious calculations.
      Ge0 = G_LSF(v); % Run the model

%% Local polynomial approximation for generated candidate sample
      if off==1
          
          check=1;
          epadd=100;
          maxite=1;
      
          while check==1

              [Ge,seeds,gseeds,run1,check,epadd] = fulllocal_v2(v,x0,seeds,gseeds,G_LSF,t,j,run1,epadd); % predict for a candidate sample

              maxite=maxite+1;

              if maxite==2

                  %fprintf('Break max iterations\n');

                  break;             

              end
          end

          err = abs((Ge-Ge0)/Ge); % Checking relative error with the real evalution    

          err1(kkk)=err; %
          
          kkk=kkk+1;
      
      else
          
          Ge = Ge0;
          
          err1 = 0;
          
          run1 = 0;
          
      end

%% Check that samples satisfy nestedness

      %Ge=Ge0;
      
      if Ge <= b && not(all(Ge <= bp(1:j-1))) && j > 1 % check nestedness of the (intermediate) failure domains
         
          Ge = G_LSF(v);
          
          fprintf("Nestedness.... %g \n", idx+t-1); % use numerical model if a sample does not satisfy nestedness
           
          run1 = run1 + 1;
          
      end

%% Check if a candidate is in failure or not
      
      if Ge <= b
          
         u_jp1(:,idx+t) = v;    % accept the candidate in failure region
         geval(idx+t)   = Ge;   % store the lsf evaluation
         acc(idx+t)     = 1;    % note the acceptance
         
         seeds(:,end+1) = v; % accepted point using within seeds
            
         gseeds(:,end+1) = Ge; % accepted evaluation using within seeds
      
      else
          
         u_jp1(:,idx+t) = u_jp1(:,idx+t-1);   % reject the candidate and use the same state
         geval(idx+t)   = geval(idx+t-1);     % store the lsf evaluation
         acc(idx+t)     = 0;                  % note the rejection
    
      end
      
      
      
   end
   % average of the accepted samples for each seed Ref: https://www.bgu.tum.de/era/software/software00/subset-simulation/
   mu_acc(i) = mu_acc(i) + min(1, mean(acc(idx+1:idx+Nchain(k)-1)));
   
   if mod(k,Na) == 0
      if (Nchain(k) > 1)   % only if the length of the chain is larger than 1
         % c. evaluate average acceptance rate
         hat_a(i) = mu_acc(i)/Na;   % Ref. 1 Eq. 25
         
         % d. compute new scaling parameter
         zeta        = 1/sqrt(i);   % ensures that the variation of lambda(i) vanishes
         lambda(i+1) = exp(log(lambda(i)) + zeta*(hat_a(i)-star_a));  % Ref. 1 Eq. 26
         
         % update parameters
         sigma = min(lambda(i+1)*sigma_0,ones(n,1));   % Ref. 1 Eq. 23
         rho   = sqrt(1-sigma.^2);                     % Ref. 1 Eq. 24
         
         % update counter
         i = i+1;
      end
   end
end
%% Ref: https://www.bgu.tum.de/era/software/software00/subset-simulation/
% next level lambda
new_lambda = lambda(i);

% compute mean acceptance rate of all chains
if i ~= 1
   accrate = mean(hat_a(1:i-1));
else
   accrate = sum(acc(1:mod(N,Ns)))/mod(N,Ns);   % almost all seeds are in the failure domain
end

err1=mean(err1); % mean value for relative error

return;
%%END

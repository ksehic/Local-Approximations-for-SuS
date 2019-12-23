function [Pf_SuS,delta_SuS,prob,b,Pf,b_line,Pf_line,samplesU,gruns] = SuS_Local(N, p0, ff,n,off,n_seed)
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
* Space: Standard Normal N(0,1)
* MCMC: Component-wise Metropolis Algorithm
* Regresion: Gaussian Process
* Dimensionality Reduction: Partial Least Square regression (PLS1)
---------------------------------------------------------------------------
Input:
* N     : Number of samples per level
* p0    : Conditional probability of each subset
* ff    : Pick your advanture
* n     : Dimension of a problem
* off   : Local approximations TURN ON/OFF
* n_seed : Seed number
---------------------------------------------------------------------------
Output:
* Pf_SuS    : failure probability estimate of subset simulation
* delta_SuS : coefficient of variation estimate of subset simulation 
* b         : intermediate failure levels
* b_line    : limit state function values
* Pf        : intermediate failure probabilities
* Pf_line   : failure probabilities corresponding to b_line
* samplesU  : samples in the Gaussian standard space for each level
* gruns     : number of evalutions for a numerical model
* prob      : intermediate probabilities
---------------------------------------------------------------------------
* Main Ref:
      
    Subset Simulation Method - Felipe Uribe (felipe.uribe@tum.de)
    https://www.bgu.tum.de/era/software/software00/subset-simulation/
    
    Papaioannou I., Betz W., Zwirglmaier K., Straub D.: MCMC algorithms for subset simulation. Probabilistic Engineering Mechanics, 41: 89-103.

    P. R. Conrad, Y. M. Marzouk,N. S. Pillai and A. Smith: Accelerating Asymptotically Exact MCMC for Computationally Intensive Models via Local Approximations.
    Journal of the American Statistical Association, Volume 111, 2016 - Issue 516, Pages 1591-1607

%}

%% Seed
%rng('shuffle'); % random shuffle for sample

rng(n_seed);

%% Pick your advanture

if ff==1
    
    G_LSF = @(x) (x(1,:)-11).^2 - (x(2,:)-6).^2; % Rosenbrock function

elseif ff==2
    
    beta = 30;

    G_LSF = @(x) -((0-0.1*x(1,:).^2).^2+0.5*(x(2,:)-x(1,:).^2).^2) + beta;  % Rosenbrock function
    
elseif ff==3 % Four-Branch Function with Pmc = 4.46 * 1e-03 and N=1e8
    
    g_fun = @(x) min([ 0.1.*(x(:,1)-x(:,2)).^2-(x(:,1)+x(:,2))./sqrt(2)+3,...
               0.1.*(x(:,1)-x(:,2)).^2+(x(:,1)+x(:,2))./sqrt(2)+3,...
               x(:,1)-x(:,2) + 7./sqrt(2),...
               x(:,2)-x(:,1) + 7./sqrt(2) ], [], 2);

    G_LSF = @(x) g_fun(x');
  
elseif ff==4 % linear limit-state function
    
    G_LSF = @(x) -sum(x)/sqrt(n) +4; % Pes = 3.17*1e-5
    
elseif ff==5
    
    G_LSF = @(x) min([3.2 + (x(1,:)+x(2,:))./sqrt(2),0.1*(x(1,:)-x(2,:)).^2-(x(1,:)+x(2,:))./sqrt(2)+2.5]);
    
elseif ff==6 % nonlinear limit-state function
    
    kk=0.2;
    
    G_LSF = @(x) -sum(x)/sqrt(n) +4 - kk/4*(x(1)-x(2))^2; % kk=0.2 Pex=6.41*1e-5 kk=1 Pex=8.99*1e-3.
    
elseif ff==7 % elliptical limit-state function
    
    r = 8; c1 = 1; c2 = 0.5; theta_x = pi/4;
    
    G_LSF = @(x) r^2 - (x(1)*cos(theta_x)+x(2)*sin(theta_x))^2/c1^2 - (x(1)*sin(theta_x)+x(2)*cos(theta_x))^2/c2^2; % r=8 Pex = 6.85*1e-3

elseif ff==8 % Nonlinear oscillator Pf_SuS 3.71 10^(-4)
    
    G_LSF = @(x) bouc_wen(x)+0.3;
    
elseif ff==9 % hypersphere
    
    G_LSF = @(x) 1-norm(x)^2/4.6^2-x(1)/52.6*((1-(norm(x)/4.6)^2)/(1+(norm(x)/4.6)^2));
    
end

%% Initialization of variables and storage

j      = 1;                % initial conditional level
Nc     = N*p0;             % number of markov chains
Ns     = 1/p0;             % number of samples simulated from each Markov chain
lambda = 0.6;              % recommended initial value for lambda
max_it = 150;               % maximum number of iterations
%
geval = zeros(1,N);        % space for the LSF evaluations
gsort = zeros(max_it,N);   % space for the sorted LSF evaluations
delta = zeros(max_it,1);   % space for the coefficient of variation
nF    = zeros(max_it,1);   % space for the number of failure point per level
prob  = zeros(max_it,1);   % space for the failure probability at each level
b     = zeros(max_it,1);   % space for the intermediate leveles

gruns = zeros(max_it,1);   % runs of model
bruns = zeros(max_it,1);   % intermediate runs
%% Initial Latin Hybercube sampling at j=0
plhmc = haltonset(n,'Skip',1e3,'Leap',1e2);
plhmc = scramble(plhmc,'RR2');
X=net(plhmc,N);
% transfer to the standard normal space
mu=0;
sigma=1;
Xg = mu +sqrt(2)*sigma*erfinv(2*X-1);
u_j=Xg.'; %%% input d x 1

%u_j = randn(n,N); % Monte Carlo sampling

%% Run Numerical Model for generated samples

fprintf('Evaluating a numerical model:\t');

% Define when to start local approximate a model
Naprox = 100;
s_runs = 99; % number of evaluations

idx_start = zeros(1,99);

for i = 1:N
   
   if i < Naprox 
        
       geval(i) = G_LSF(u_j(:,i)); % evaluating a numerical model
   
   else 
       
    if off==1

        [Z0,~,ZUP0] = start_approx(u_j(:,1:i-1).',geval(1:i-1).',u_j(:,i).',G_LSF); % local approximate a numerical model

        if abs(diff(ZUP0))/abs(Z0)*100<5 % check uncertainty in prediction

            geval(i) = Z0;

            s_runs=s_runs+0;
            
            idx_start(i)=1;

        else % if prediction is uncertain, use a numerical model

            geval(i) = G_LSF(u_j(:,i)); 

            s_runs=s_runs+1;
            
            idx_start(i)=0;

        end
    
    else

        geval(i) = G_LSF(u_j(:,i)); % evaluating a numerical model when local approx. is turn off
        
    end

   end
   
   if off~=1
       
   geval(i) = G_LSF(u_j(:,i));

   end
   %fprintf("Number of iteration: %g of %g with num evalutions %g \n",i,N,s_runs);

end

fprintf('Done! \n');

%% Fix failure threshold at j=1 
if off==1
    
    fprintf('Start failure fix...\n');

    [u_j,geval,s1_runs] =  start_fix(u_j,geval,G_LSF,p0,idx_start);

    s_runs = s_runs + s1_runs;

    fprintf('Done!\n');

end
%% Collect data
totalu = u_j;
totalg = geval;

%% Subset Simulation Method
while true 
   [gsort(j,:), idx] = sort(geval);     % sort values in ascending order
   
   u_j_sort          = u_j(:,idx);      % order the samples by idx
   samplesU.total{j} = u_j_sort;        % store the ordered samples
   
   b(j) = prctile(geval, p0*100);       % estimate intermediate failure level at j for p0
   
   if j>1 && off==1 && b(j) > 0         % correction intermediate level
      
       [u_j,geval,brun,totalu,totalg] = blevel_fix(u_j,geval,G_LSF,N,p0,ZUP,totalu,totalg);
       
       bruns(j) = brun;
       
       [gsort(j,:), idx] = sort(geval);     % sort values in ascending order

       u_j_sort          = u_j(:,idx);      % order the samples by idx
       samplesU.total{j} = u_j_sort;        % store the ordered samples
       
       b(j) = prctile(geval, p0*100);       % estimate updated intermediate failure level at j for p0

   end
   
   nF(j) = sum(geval <= max(b(j), 0)); % number of failure points in the next failure level
   
   %% Assign conditional probability to the level
   if b(j) <= 0   
       
      b(j)    = 0;
      
      prob(j) = nF(j)/N; % Monte Carlo esimation of failure probability at final level
      
      if off==1 % Fix final failure probability
          
          % Ref:  J. Li, D. Xiu, Evaluation of failure probability via surrogate models, Journalof Computational Physics 229 (23) (2010) 8966–8980.

          [pfupdate,lastruns] = lastfix(u_j_sort,gsort(j,:),b(j),N,G_LSF); %% correct failure point

          fprintf('\n-Threshold intermediate level %g = %g \n', j-1, b(j));

          fprintf('\n- Pf correction at the end %g \n', lastruns);

          prob(j) = pfupdate;

          gruns(j,1) = lastruns;
          
      end

   else
       
      prob(j) = p0;
      
      fprintf('\n-Threshold intermediate level %g = %g \n', j-1, b(j));
   
   end
   
   %% compute coefficient of variation
   % Ref:     Papaioannou I., Betz W., Zwirglmaier K., Straub D.: MCMC algorithms for subset simulation. Probabilistic Engineering Mechanics, 41: 89-103.
   
   if j == 1
      
       delta(j) = sqrt(((1-p0)/(N*p0)));
      
   else
       
      I_Fj     = reshape(geval <= b(j),Ns,Nc);
      p_j      = (1/N)*sum(I_Fj(:));                   
      gamma    = corr_factor(I_Fj,p_j,Ns,Nc);          
      delta(j) = sqrt( ((1-p_j)/(N*p_j))*(1+gamma) );  
      eff = (1+gamma)^(-1);
      
   end
   
   %% Select failure samples for MCMC
   
   samplesU.seeds{j} = u_j_sort(:,1:nF(j));   % store ordered level seeds
  
   geval_0 = gsort(j,1:nF(j)); % store ordered evalutions
   
   idx_rnd   = randperm(nF(j));  % randomize the ordering of the samples (to avoid bias)
   
   rnd_seeds = samplesU.seeds{j}(:,idx_rnd);   % non-ordered seeds
   
   rnd_geval = geval_0(:,idx_rnd);
   
   %% MCMC chains with Local Approximations
   
   if b(j)>0
   
       [u_j, geval, lambda, sigma, acc,run,err1,totalu,totalg,ZUP] = aCSback(N, lambda, b(j), b(1:j), rnd_seeds, rnd_geval,G_LSF,totalu,totalg,j,off);
        
       gruns(j) = run;
       fprintf('\t*aCS lambda = %g \t*aCS sigma = %g \t *aCS accrate = %g \n', lambda, sigma(1), acc);
       fprintf('The total number of extra runs %g which is %g percent \n',run,run/N*100);
       fprintf('With the relative mean error %g \n',err1);
     
   end

   %% next failure level
   j = j+1;   
   
   if b(j-1) <= 0 || j-1 == max_it % check if final failure region is reached
       
      break;
      
   end
   
end

m = j-1;
samplesU.total{j} = u_j;   % store final failure samples (non-ordered)

% delete unnecesary data
if m < max_it
   gsort(m+1:end,:) = [];
   prob(m+1:end)    = [];
   b(m+1:end)       = [];
   delta(m+1:end)   = [];
end

%% probability of failure
% failure probability estimate
Pf_SuS = prod(prob);   % or p0^(m-1)*(Nf(m)/N);

% coeficient of variation estimate
delta_SuS = sqrt(sum(delta.^2)); 

%% Probability of failure estimation
Pf           = zeros(m,1);
Pf(1)        = p0;

Pf_line(1,:) = linspace(p0,1,Nc);
b_line(1,:)  = prctile(gsort(1,:),Pf_line(1,:)*100);

for i = 2:m
    
   Pf(i)        = Pf(i-1)*p0;
   Pf_line(i,:) = Pf_line(i-1,:)*p0;
   b_line(i,:)  = prctile(gsort(i,:),Pf_line(1,:)*100);
   
end

% different failure thresholds and correspoding failure probabilities
Pf_line = sort(Pf_line(:));

b_line  = sort(b_line(:));

%% Couting evalutions

if off==1 

    gruns(end+1) = sum(bruns);

    gruns(end+1) = s_runs;

end

return
%% THE END
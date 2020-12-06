Local Polynomial Approximations for Subset Simulation method
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
* Space: Standard Normal N(0,1)
* MCMC: Component-wise Metropolis Algorithm
* Regresion: Polynomial Quadratic
* Dimensionality Reduction: /
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

---------------------------------------------------------------------------

1. You will need to have MATLAB2019B or similar.
3. Download files to your folder.
4. Open SuS_poly.m
5. Follow comments
6. Run and explore

If you find any mistake or bug or you have comments, please contact me at kense@dtu.dk or kenosehic@gmail.com.

Thank you!

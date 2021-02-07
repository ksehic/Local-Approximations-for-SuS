function [Z1,ZSD1,ZUP,R,Rpls11,XmeanLoc,gpmodel1] = localpls_update(X,Y1,X01,g)
%% Local Learning based on Gaussian process for Subset Simulation method
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
* Train Gaussian process on low-dimensional subspace span by PLS projections
---------------------------------------------------------------------------
Input:
* X             : Samples
* Y1            : Evalutions
* X01           : MCMC proposal
* g             : Numerical model
---------------------------------------------------------------------------
Output:
* Z1            : Prediction
* ZSD1          : Standard deviation
* ZUP           : Confidence interval
* R             : Radius of ball
* Rpls11        : Low-dimensional projection matrix
* XmeanLoc      : Mean value for initial samples within ball
* gpmodel1      : trained Gaussian process
---------------------------------------------------------------------------
* Main Ref:
      
    M. Bouhlel, N. Bartoli, A. Otsmane, J.J. Morlier: Improving  kriging surrogates of high-dimensional design models 
    by partial least squaresdimension  reduction, Structural  and  Multidisciplinary  Optimization  53(2016) 935–952.

    I. Papaioannou, M. Ehre, D. Straub, PLS-based adaptation for efficient PCA representation in high dimensions, 
    Journal of Computational Physics 387(2019) 186–204.

%}

%% Procedure to train Gaussian process on low-dimensional subspace span by PLS projections

epsilon=1e-2; % tolerance

[Rpls1,Tr1] = pls1(X,Y1,epsilon); % estimate global low-dimensional subspace

Xmean=mean(X);

X1 = X01 - Xmean;

T1 = X1*Rpls1;

%% Euclidean distance for samples

ED =@(Q,D) sqrt(sum((Q-D).^2,2)); 

ED1 = ED(T1,Tr1); % estimate distance between samples and MCMC proposal

[EDS1,id1]=sort(ED1);

d1=length(T1(1,:));

Ndef1=d1+1; % Number of design samples for Gaussian process

i=0;

[Rpls11,Tr1] = pls1(X(id1(1:Ndef1+i),:),Y1(id1(1:Ndef1+i)),epsilon); % estimate local low-dimensional subspace

Y11 = Y1(id1(1:Ndef1+i));

%% Train Gaussian process on local low-dimensional subspace

gpmodel1 = fitrgp(Tr1,Y11,'KernelFunction','ardsquaredexponential','BasisFunction','constant','FitMethod','exact','PredictMethod','exact');

%% Make prediction

XmeanLoc = mean(X(id1(1:Ndef1+i),:));

T11 = X01-XmeanLoc;

T11R = T11*Rpls11; % project MCMC proposal to local low-dimensional subspace

[Z1,ZSD1,ZUP] = predict(gpmodel1,T11R);

%% Estimate radius

R = EDS1(Ndef1+i);
       
return


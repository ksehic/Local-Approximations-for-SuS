function [R,Tr] = pls1(X,Y,epsilon)
%% Local Learning based on Gaussian process for Subset Simulation method
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
*   PLS1 - Partial least squares regression
---------------------------------------------------------------------------
* Main Ref:
      
    I. Papaioannou, M. Ehre, D. Straub, PLS-based adaptation for efficient pce650representation in high dimensions, 
    Journal of Computational Physics 387(2019) 186–204.

%}
%% Procedure

X = X - mean(X); % centering

Y = Y - mean(Y); % centering

E = X; F = Y; i=1;

[xr,xc]=size(X);

W = zeros(xc,xr);

T = zeros(xr,xc);

P = zeros(xc,xr);

while true
   
   F1 = norm(F); 
    
   W(:,i) = (E.'*F)/(norm(E.'*F)); %compute weights
   
   W(:,i) = W(:,i)/norm(W(:,i)); % normalization

   T(:,i) = E*W(:,i); % compute score
   
   P(:,i) = (E.'*T(:,i))/(T(:,i).'*T(:,i));% compute load
   
   b = (T(:,i).'*F)/(T(:,i).'*T(:,i)); % compute regressor
   
   E = E - T(:,i)*P(:,i).'; F = F - b*T(:,i); % deflate
   
   i=i+1;
   
   if norm(F)<epsilon || i==length(X(1,:))+1 || (round(norm(F),2)==round(norm(F1),2)) % error check
       
       break;
       
   end
    
end

W(:,i:end) = [];

P(:,i:end) = [];

T(:,i:end) = [];

R = W*(P.'*W)'; % an alternative way to represent the weights

Tr = X*R;

return

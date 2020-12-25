function [IRF] = IRFar1(phi,mu,sigma,T)
%--------------------------------------------------------------------------
% Purpose: Calculate the IRF of an AR(1) for different values of phi
%--------------------------------------------------------------------------
% Inputs:   phi    : 1x1 coefficient of the AR(1) process
%           mu     : 1x1 mean of the white noise error term
%           sigma  : 1x1 variance of the white noise process error term 
%           T      : 1x1 sample size
%--------------------------------------------------------------------------
% Output    py     : Tx1 IRF function
%--------------------------------------------------------------------------


noise    = sigma*randn(1,T)+mu;
e        = noise';
y        = NaN(T,1);
y(1,1)   = e(1,1);
IRF      = NaN(T,1);


for i=1:T
    IRF(i,1)  = phi^(i) 
end 




return 

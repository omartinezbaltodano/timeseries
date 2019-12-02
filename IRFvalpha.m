function [IRF,t] = IRFvalpha(phi1, phi2,alpha,beta,J)
%--------------------------------------------------------------------------
% Purpose :   This function gives the IRF for the AR(3) process
%--------------------------------------------------------------------------
% Inputs  : phi1  : 1x1  first coefficient of the productivity process
%           phi2  : 1x1  second coefficient of the productivity process
%           alpha : 1x1  share of the capital in the income
%           beta  : 1x1  Discount factor
%           j     : 1x1  number of periods for the IRF
%--------------------------------------------------------------------------
% Output  : Vstar : Jx1  IRF for J periods
%--------------------------------------------------------------------------

gamma0 = (1-phi1-phi2)*log(alpha*beta);
gamma1 =  alpha+phi1;
gamma2 = -(phi1*alpha-phi2);
gamma3 =  -alpha*phi2;

F = [gamma1 gamma2  gamma3; 1  0  0; 0  1  0]
e = eig(F)

% valores propios
lambda1 = e(1,1)
lambda2 = e(2,1)
lambda3 = e(3,1)



% Coeficientes C
c1 = lambda1^2/[(lambda1-lambda2)*(lambda1-lambda3)]
c2 = lambda2^2/[(lambda2-lambda1)*(lambda2-lambda3)]
c3 = lambda3^2/[(lambda3-lambda1)*(lambda3-lambda2)]

%IRF
%J   = 40;
IRF = NaN(J,1);
t   = NaN(J,1);
t(1,1) = 1
for i=1:J
    IRF(i,:) = c1*lambda1^i+c2*lambda2^i+c3*lambda3^i
    t(i,:)   = i
end 



return

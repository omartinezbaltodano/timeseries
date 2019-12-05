function [IRF, R] = IRFar2_complex(gamma1, gamma2,J)
%--------------------------------------------------------------------------
% Purpose :   This function gives the IRF for the AR(2) process when the
%             roots are complex
%--------------------------------------------------------------------------
% Inputs  : phi1  : 1x1  first coefficient of process (t-1 component)
%           phi2  : 1x1  second coefficient of process (t-2 component)
%--------------------------------------------------------------------------
% Output  : IRF : Jx1  IRF for J periods
%--------------------------------------------------------------------------
i = sqrt(-1)
F = [gamma1 gamma2; 1  0]
e = eig(F)

% valores propios
lambda1 = e(1,1)
lambda2 = e(2,1)
a       = 0.5*gamma1
b       = 0.5*sqrt(-gamma1^2-4*gamma2)
R       = sqrt(a^2+b^2)  
theta   = acos(a/R)                        %acos: inversa de coseno


%IRF
%J   = 40;
IRF = NaN(J,1);
t   = NaN(J,1);
t(1,1) = 1
for i=1:J
    IRF(i,:) = R^i*cos(theta*i)+ [cos(theta)/sin(theta)]*R^i*sin(theta*i)
    t(i,:)   = i
end 



return

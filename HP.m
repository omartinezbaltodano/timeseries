function Y = HP(X,freq,tipo)
%--------------------------------------------------------------------------
%Proposito: Aplica el filtro Hodrick-Prescott  para calcular el componente 
%           ciclico de una serie de tiempo
%--------------------------------------------------------------------------
% Output  : Y          : TxK, Componente ciclico de la serie
%--------------------------------------------------------------------------
% Input   : X          : TxK  Serie original
%           tipo       : 1x1  Toma valor cero si el componente cíclico está 
%                             determinado por el parámetro lambda (predeterminado)
%                             Toma valor uno si el componente cíclico está 
%                             determinado por la frecuencia de corte
%           freq       : 1x1  En caso que el tipo sea cero: 
%                             valores sugeridos son:
%                             6       para datos anuales
%                             1600    para datos trimestrales
%                             129600  para datos mensuales
%                             En caso que el tipo sea uno: 
%                             lambda = (2*sin(pi/freq))^-4
%                             freq   = pi/arcsin(1/2*lambda^-1/4)
%--------------------------------------------------------------------------

T               = size(X,1);

I               = speye(T);
LT              = spdiags(ones(T-1,1),-1,T,T);
LT              = (I-LT)^2;
Q               = LT(3:end,:)';
SIGMA_R         = Q'*Q;
SIGMA_T         = speye(T-2);
 
lambda          = freq;
if nargin>2 && tipo ==1
    lambda      = (2*sin(pi/freq))^-4;
end
    
g               = Q'*X;
b               = (SIGMA_T+lambda*SIGMA_R)\g;
Y               = lambda*Q*b;
end 
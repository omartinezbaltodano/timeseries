function  Y  = BK(X,lambda_min,lambda_max,o)
%--------------------------------------------------------------------------
%Proposito: Aplica el filtro Baxter-King para calcular el componente 
%           ciclico de una serie de tiempo
%--------------------------------------------------------------------------
% Output  : Y          : (T-2o)xK, Componente ciclico de la serie
%--------------------------------------------------------------------------
% Input   : X          : TxK  Serie original
%           lambda_min : 1x1  Minimo numero de periodos por ciclo
%           lambda_max : 1x1  Maximo numero de periodos por ciclo
%           o          : Orden del Filtro (sugerencia: 12)
%--------------------------------------------------------------------------

lambda_min  = max(2,lambda_min);

omegamin    = 2*pi/lambda_max;
omegamax    = 2*pi/lambda_min;

J           = [1:1:o]';
B           = (sin(omegamax*J)-sin(omegamin*J))./(pi*J);
B           = [(omegamax-omegamin)/pi;B];

%Baxter-King approximation
s           = B(1) + 2*sum(B(2:end));
BB          = B - s/(2*o+1);
BB          = [flipud(BB(2:end));BB];

%filtering
T           = size(X,1);
for t = o+1:1:T-o
    Y(t,:)  = BB'*X(t-o:t+o,:);    
end
Y(1:o,:)=[];
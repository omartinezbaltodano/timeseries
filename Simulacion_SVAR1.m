%--------------------------------------------------------------------------
% Simulación de series y luego estimación del modelo VAR estructural
% Propósito: Este programa está hecho para mostrar a los alumnos de MC3 en 
%            la Universidad de Chile, como simular series (DGP) para 
%            estimar un VAR(1). Luego a partir de estas series simuladas, 
%            se estima el proceso de forma reducida y se calculan las IRF 
%            del DGP y del proceso estimado
%--------------------------------------------------------------------------

clc
clear
close all

%% Creando el DGP

% Largo del periodo
T = 150;
% Numero de variables
k = 2;  

% Parametros de la forma estructural (Si el proceso es triangular)
B       = [2.5, 0; -1.19, 1]; 
% Si quiere simular un proceso donde la matriz B es no triangular, puede
% usar la siguiente B
% B     = [1, -5.12; 2.19, 1]; % 
Phi_1   = [0.6, -0.3; -0.3, 0.6];

% Los parametros de la forma reducida serian:
inv_B   = B^(-1);
A_uno   = inv_B*Phi_1;

% Verificar si el proceso es estacionario
[V, C] = eig(A_uno);
eigenvalues = max(C); 
dominant_root = max(max(abs(C)))
assert( abs(dominant_root) <1) 

% Simulando shocks estructurales
% Fijando la semilla
randn('seed',0)                  
% Shocks ortogonales
varepsilon = mvnrnd(zeros(k,1),eye(k),T)'; 

% Los chocks de forma reducida
u = inv_B*varepsilon; 

% Generando los datos. Con punto de partida cero
Y = zeros(k,T); 
Y(:,1) =  u(:,1);
for t=2:T
    Y(:,t) = A_uno*Y(:,t-1) + u(:,t);
end

%% Calculando las verdaderas IRF a partir del DGP

T_irf    = 10;
IRF_true = zeros(2*k,T_irf);

epsilon1 = 1*std(varepsilon(1,:));
epsilon2 = 1*std(varepsilon(2,:));

v = inv_B*[epsilon1; 0];  
IRF_true(1:2,1) = v;
for t = 2:T_irf
    IRF_true(1:2,t) = A_uno*IRF_true(1:2,t-1);
end

v = inv_B*[0; epsilon2];  
IRF_true(3:4,1) = v;
for t = 2:T_irf
    IRF_true(3:4,t) = A_uno*IRF_true(3:4,t-1);
end

%% Estimando las IRF

X = Y(:,1:T-1)'; 
y = Y(:,2:T)';

% estimacion
A_uno_hat = (X'*X)^(-1)*X'*y;
Sigma_hat = (y-X*A_uno_hat)'*(y-X*A_uno_hat)/T

A_uno_hat = A_uno_hat';
A_uno_hat
A_uno

Sigma_hat
Sigma = inv_B*inv_B'

inv_B_hat = chol(Sigma_hat)'
inv_B_hat
inv_B

% IRFs
IRF_hat = zeros(2*k,T_irf);
T_irf    = 10;
v = inv_B_hat*[epsilon1; 0];  
IRF_hat(1:2,1) = v;
for t = 2:T_irf
    IRF_hat(1:2,t) = A_uno_hat*IRF_hat(1:2,t-1);
end

v = inv_B_hat*[0; epsilon2];  
IRF_hat(3:4,1) = v;
for t = 2:T_irf
    IRF_hat(3:4,t) = A_uno_hat*IRF_hat(3:4,t-1);
end

%% Graficas

figure(1)
 
    subplot(2,2,1)
plot(IRF_true(1,:), 'linewidth', 2), hold on
plot(IRF_hat(1,:), '-* r'), hold on
plot(zeros(T_irf,1), 'k')
title('IRF of Variab. 1, shock 1'), axis([1 T_irf -Inf Inf]), set(gca,'box','off')
legend('Verdadera', 'Estimatada', 'Location', 'NorthEast'), legend('boxoff')

    subplot(2,2,2)
plot(IRF_true(2,:), 'linewidth', 2), hold on
plot(IRF_hat(2,:), '-* r'), hold on
plot(zeros(T_irf,1), 'k')
title('IRF of Variab. 2, shock 1'), axis([1 T_irf -Inf Inf]), set(gca,'box','off')

    subplot(2,2,3)
plot(IRF_true(3,:), 'linewidth', 2), hold on
plot(IRF_hat(3,:), '-* r'), hold on
plot(zeros(T_irf,1), 'k')
title('IRF of Variab. 1, shock 2'), axis([1 T_irf -Inf Inf]), set(gca,'box','off')
    
    subplot(2,2,4)
plot(IRF_true(4,:), 'linewidth', 2), hold on
plot(IRF_hat(4,:), '-* r'), hold on
plot(zeros(T_irf,1), 'k')
title('IRF of Variab. 2, shock 2'), axis([1 T_irf -Inf Inf]), set(gca,'box','off')

set(gcf, 'PaperPositionMode', 'auto');
set(gcf,'Position',[0 0 1000*.7 1200*.7])
%print('-painters', '-dpdf','-r600', 'D:\Dropbox\..');




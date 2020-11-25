function LLN_AR1(n,phi)
%-------------------------------------------------------------------------- 
% Propósito:  Esta función muestra cómo funciona la ley de los grandes 
%             números y el teorema central del límite para un AR(1) 
%                       y_t = phi y_{t-1} + epsilon_t 
%             epsilon_t tiene distribución normal estándar 
%-------------------------------------------------------------------------- 
% Inputs   :  n     : 1x1  tamaño de la muestra  
%             mu    : 1x1  media de la distribución del error 
%             sigma : 1x1  varianza de la distribución del error 
%             phi   : 1x1  coeficiente del AR(1) 
%-------------------------------------------------------------------------- 
% Output   :  graph1 : grafico que muestra como la media muestral converge 
%             graph2 : grafico que muestra la distribución de la media de y          
%-------------------------------------------------------------------------- 

%% Ley de los grandes números
epsilon = randn(n,1);

y       = NaN(n,1);
y(1,1)  = epsilon(1,1);
for t=2:n
    y(t,1) = phi*y(t-1,1)+epsilon(t,1);
end 

ybar = ones(1,n)*y/n;

xbar_iter = NaN(n,1);

for i=1:n
    step1        = ones(1,i)*y(1:i,1)/i;
    xbar_iter(i) = step1;
end
clear step1

figure(1)
plot(xbar_iter), title('Promedios Muestrales: $\frac{1}{n}\sum_{i=1}^n x_i$','Interpreter','latex','FontSize',16), xlabel({'n'}),  legend('Media Muestral','Media Poblacional')

%% Teorema del límite central

y           = NaN(n,n);
approx      = NaN(n,1);
y(1,1:end)  = epsilon(1,1);

for j=1:n
    for t=2:n
        y(t,j) = phi*y(t-1,1)+randn(1,1);
    end 
end 
    
iter = n;
for i=1:iter
    step1     = ones(1,n)*y(:,i)/n; 
    step2     = n^(0.5)*(step1);
    approx(i) = step2;
end

clear step1 step2

figure(3) 
hist(approx,80), title('Densidad de la media muestral estandarizada','FontSize',16) 


end

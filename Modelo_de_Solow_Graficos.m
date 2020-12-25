clc
clear all

%% Parametos
A0      = 2
gammaA  = 0.9
rho     = 0.1
T       = 50
mu      = 0
sigma   = 1
epsilon = normrnd(mu,sigma,T,1)
z       = NaN(T,1)
z(1,1)  = 0
A       = NaN(T,1)

%% Primer Grafico
for i=2:T
    z(i,1) = rho*z(i-1,1)+epsilon(i,1)
end 


for i=1:T
    A(i,1) = A0*(1+gammaA)^(i)*exp(z(i,1))
end


t = NaN(T,1)
for i=1:T
    t(i,1) = i
end 

plot(t,A), xlabel({'t'}), title({'Proceso de Productividad'}), ylabel({'A_t'}) 





%% Grafico 2
clc
clear all

% Parametos
A0      = 2
gammaA  = 0
rho     = 0
T       = 500
mu      = 0
sigma   = 1
epsilon = normrnd(mu,sigma,T,1)
z       = NaN(T,1)
z(1,1)  = 0
A       = NaN(T,1)


for i=2:T
    z(i,1) = rho*z(i-1,1)+epsilon(i,1)
end 


for i=1:T
    A(i,1) = A0*(1+gammaA)^(i)*exp(z(i,1))
end

t = NaN(T,1)
for i=1:T
    t(i,1) = i
end 

plot(t,A), xlabel({'t'}), title({'Proceso de Productividad'}), ylabel({'A_t'}) 




%% IRF del modelo 

clc
clear all
% parametos
delta   =  0.01;
Sigma   =  0.2;
n       =  0.02;
A_0     =  1;
T       =  1000;          
theta   =  0.36;

Omega   = (delta+n)/(1+n);
phi     = (1+theta*n+delta*(1-theta))/(1+n);


IRF = NaN(T,1)
t   = NaN(T,1)

for i=1:T
    IRF(i,1) = phi^(i-1)*Omega
    t(i,1)   = i-1
end 

plot(t,IRF), xlabel({'t'}), title({'IRF'}) 






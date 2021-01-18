
cd('/home/octavio/Documentos/Analisis Espectral/Matlab')

clc; clear all;
close all;
%% Calibracion
omega = [-pi:0.01:pi];
i     = sqrt(-1);

%% Coeficiente del AR(2)
rho_1=[1.25,1.3];
rho_2=[-0.3,-0.4];
 
for l=1:2 
    for j=1:length(omega)
        spectrum(l,j)=1/(2*pi*(1-rho_1(l)*exp(-i*omega(j))-rho_2(l)*exp(-2*i*omega(j)))...
            *(1-rho_1(l)*exp(i*omega(j))-rho_2(l)*exp(2*i*omega(j))));
    end
end


figure(1)
plot(omega,spectrum(1,:),'LineWidth',2)
hold on
plot(omega,spectrum(2,:),'-.','LineWidth',2)
title('Spectral Density (in quarters per cycle)')
xlabel('Quarters per cycle (\omega)');
ylabel('Spectrum (s_x(\omega))');
legend('\phi_1=1.25; \phi_2= - 0.3','\phi_1=1.3; \phi_2= - 0.4')
x0      = 0.15
y0      = 0.2
width   = 9
height  = 5.9
figure('Units','inches','Position',[x0 y0 width height],'PaperPositionMode','auto','PaperSize',[width (height-0.5)]);
plot(omega,spectrum(1,:),'LineWidth',2)
hold on
plot(omega,spectrum(2,:),'-.','LineWidth',2)
title('Spectral Density (in quarters per cycle)')
ylabel('Spectrum (s_x(\omega))'); 
xlabel('Quarters per cycle (\omega)');
legend('\rho_1=1.25; \rho_2= - 0.3','\rho_1=1.7; \rho_2= - 0.8');
print('Graph1', '-dpdf', '-r0');

%% Autocovarianzas y autocorrelaciones
jgrid = [-30:1:30];
for n=1:2
    for k = 1:length(jgrid)
        j = jgrid(k);
 %       tmp = zeros(omega,1);
        for m = 1:length(omega)
            w_m = omega(m);
            tmp2(n,m) = exp(i*w_m*j)*spectrum(n,m);
        end
        covars(n,k) = sum(tmp2(n,:));
    end
    for k = 1:size(jgrid,2)
        j = jgrid(k);
        corrs(n,k) = covars(n,k)/covars(jgrid==0);
    end
end


figure(2)
plot(jgrid,corrs(1,:),'LineWidth',2)
title('Auto-Correlation Function')
xlabel('Step-size (j)');
ylabel('\rho_j');
legend('\phi_1=1.25; \phi_2= - 0.3')    % eliminar


x0      = 0.15
y0      = 0.2
width   = 9
height  = 5.9
figure('Units','inches','Position',[x0 y0 width height],'PaperPositionMode','auto','PaperSize',[width (height-0.5)]);
plot(jgrid,corrs(1,:),'LineWidth',2)
%hold on
%plot(jgrid,corrs(2,:),'-.','LineWidth',2)
title('Auto-Correlation Function')
xlabel('Step-size (j)');
ylabel('\phi_j');
legend('\phi_1=1.25; \phi_2= - 0.3')
%legend('\rho_1=1.25; \rho_2= - 0.3','\rho_1=1.7; \rho_2= - 0.8')
print('Graph2', '-dpdf', '-r0');
















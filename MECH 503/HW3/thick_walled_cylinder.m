%% Compute the stress in a thick-walled cylinder under pressure

clear all, close all, clc

r2= 0.05; %m
r1= 0.025; %m
p1= 0;%5e6; % MPa
p2= -5e6; %MPa

r=r1:0.001:r2;

sigma_r = (r1^2*r2^2*(p2-p1)/(r2^2-r1^2) *(1./r.^2)) + ((r1^2*p1-r2^2*p2)/(r2^2-r1^2));

sigma_t = -r1^2*r2^2*(p2-p1)/(r2^2-r1^2) *(1./r.^2) + + (r1^2*p1-r2^2*p2)/(r2^2-r1^2);
figure(1)
plot((r-r1)/(r2-r1),sigma_r/p2,'-r','LineWidth',2.0)
hold on
plot((r-r1)/(r2-r1),sigma_t/p2,'-b','LineWidth',2.0)

grid on
grid minor
%title('Case 1: Internal Pressure')
title('Case 2: Remote Biaxial Pressure')
legend('\sigma_r','\sigma_\theta')
legend boxoff
xlabel('Dimensionless Radius')
ylabel('Dimensionless Stress ')
set(gca,'FontSize',18,'FontWeight', 'bold')
print('-f1','ThickWalledCylStresses','-dpng')


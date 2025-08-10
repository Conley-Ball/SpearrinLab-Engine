clc; close all; clear;

% Inputs
Pc = 400; % psi
Pa = 12.7;
OF = 0.75:0.05:4;

A_e = (pi*6.5^2/4)*0.0254^2; % m^2 (do not change)
T = 5000/3 * 4.448; % N

C_star_eff = 0.9;
C_F_eff = 0.95;

D_e = 4.1;
A_e = (pi*D_e.^2/4)*0.0254^2;

% Setup
for i=1:length(OF)
    data = CEA('problem','rocket','equilibrium','o/f',OF(i),'p(psi)',Pc+Pa,'reactants','fuel','RP-1','wt%',100,'t(k)',298.15,'oxid','O2(L)','wt%',100,'t(k)',90.0,'output','transport','mks','end');
    Tc_RP(i) = data.output.eql.temperature(2);
    data = CEA('problem','rocket','equilibrium','o/f',OF(i),'p(psi)',Pc+Pa,'reactants','fuel','C2H5OH(L)','wt%',75,'t(k)',298.15,'fuel','H2O(L)','wt%',25,'t(k)',298.15,'oxid','O2(L)','wt%',100,'t(k)',90.0,'output','transport','mks','end');
    Tc_eth(i) = data.output.eql.temperature(2);
    
    clc
    fprintf('%.0f%%\n',i/length(OF)*100)
end
figure(1)
hold on
plot(OF,Tc_RP,'LineWidth',2)
plot(OF,Tc_eth,'LineWidth',2)
% plot(OF,Tc_eth*0.9^2,'k','LineWidth',1,'LineStyle','--')
% plot(OF,Tc_eth*0.8^2,'k','LineWidth',1,'LineStyle','--')
% plot(OF,Tc_RP*0.9^2,'k','LineWidth',1,'LineStyle','--')
% scatter(0.85, 1100,75,'k','x','LineWidth',2)
% scatter(1.2, 1423,75,'k','x','LineWidth',2)
% % quiver(0.85,1100, 1.2,1423,'off')
% scatter(1, 2148,75,'k','x','LineWidth',2)
% scatter(1.5, 2481,75,'k','x','LineWidth',2)
% scatter(2.35, 2671,75,'k','LineWidth',2)
% scatter(1.5, 1790,75,'k','LineWidth',2)
hold off
legend('Kerosene','Ethanol')
xlabel('OF Ratio')
ylabel('Throat Temperature (K)')
grid on
clc; close all; clear;

% Inputs
Pc = 200:5:700; % psi
Pa = 12.7;
OF = 2.4;

% A_e = (pi*6.5^2/4)*0.0254^2; % m^2 (do not change)
T = 5000/3 * 4.448; % N

C_star_eff = 0.9;
C_F_eff = 0.95;

% Setup
for i=1:length(Pc)
    data = CEA('problem','rocket','equilibrium','o/f',OF,'p(psi)',Pc(i)+Pa,'pi/p',(Pc(i)+Pa)/Pa,'reactants','fuel','RP-1','wt%',100,'t(k)',298.15,'oxid','O2(L)','wt%',100,'t(k)',90.0,'output','transport','mks','end');
      isp(i) = data.output.eql.isp_vac(end);  
    % mdot_old = 0;
    % if i ==1
    % mdot(i) = 1;
    % else
    % mdot(i) = mdot(i-1);
    % end
    % A_t = 1e-3;
    % 
    % tol = 1e-5;
    % err = 1;
    % while abs(err)>tol
    % 
    %     data = CEA('problem','rocket','equilibrium','fac','ma,kg/s',mdot(i),'o/f',OF,'p(psi)',Pc(i)+Pa,'pi/p',(Pc(i)+Pa)/Pa,'reactants','fuel','RP-1','wt%',100,'t(k)',298.15,'oxid','O2(L)','wt%',100,'t(k)',90.0,'output','transport','mks','end');
    %     C_star = data.output.eql.cstar(1);
    %     C_F = data.output.eql.cf(end);
    %     isp(i) = data.output.eql.isp(end);
    %     mdot(i) = T/(C_star*C_star_eff*C_F*C_F_eff);
    %     A_t = C_star*C_star_eff*mdot(i)/((Pc(i)+Pa)*6894.7);
    %     err = mdot_old-mdot(i);
    %     mdot_old = mdot(i);
    % end
    % A_e(i) = data.output.eql.aeat(end)*A_t;
    % D_e(i) = sqrt(A_e(i)/pi)*2/0.0254
    clc
    fprintf('%.0f%%\n',i/length(Pc)*100)
end

figure(1)
plot(Pc,mdot)
xlabel('Pc (psi)')
ylabel('mdot (kg)')

figure(2)
plot(Pc,isp*0.9*0.95,'LineWidth',2)
xlabel('Pc (psi)')
ylabel('Specific Impulse (m/s)')
grid on

Pc = Pc';
mdot = mdot';
D_e = sqrt(A_e/pi)*2/0.0254
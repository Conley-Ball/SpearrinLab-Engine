clc; close all; clear;

Pa = 22.7;
Pc = 200:150:5000;

for j=1:length(Pc)
    step = 0.1;
    O_F_i = 2.5;
    for k=1:4
        O_F_i = O_F_i-5*step;
        for i=1:11
            O_F(i) = O_F_i + step*(i-1);
            RKT1=CEA('problem','rocket','equilibrium','o/f',O_F(i),'p(psi)',Pc(j),'pi/p',Pc(j)/Pa,'reactants','fuel','RP-1','wt%',100,'t(k)',298.15,'oxid','O2(L)','wt%',100,'t(k)',90.0,'output','transport','mks','end');
            C_star = RKT1.output.eql.cstar(1);
            % Isp(i) = RKT1.output.eql.isp(3);
            Isp(i) = RKT1.output.eql.cstar(1);
            Tc = RKT1.output.eql.temperature(1);
        end
        [~,index] = max(Isp);
        O_F_i = O_F(index);
        step = step/10;
    end
    peak_of(j) = O_F_i;
    clc
    fprintf('%.0f%%\n',j/length(Pc)*100)
end

figure(1)
plot(Pc,peak_of,'LineWidth',2)
xlabel('Chamber Pressure (psi)')
ylabel('Optimal O/F Ratio')
title('Optimal O/F vs. Chamber Pressure')
grid on

% for j=1:length(Pc)
%     for i=1:length(O_F)
%         RKT1(j)=CEA('problem','rocket','equilibrium','o/f',O_F(i),'p(psi)',Pc(j),'pi/p',Pc(j)/Pa,'reactants','fuel','RP-1','wt%',100,'t(k)',300,'oxid','O2(L)','wt%',100,'t(k)',90.0,'output','transport','mks','end');
%         C_star(i,j) = RKT1(j).output.eql.cstar(1);
%         Isp(i,j) = RKT1(j).output.eql.isp(3);
%         Tc(i,j) = RKT1(j).output.eql.temperature(1);
%         RKT2(j)=CEA('problem','rocket','equilibrium','o/f',O_F(i),'p(psi)',Pc(j),'pi/p',Pc(j)/Pa,'reactants','fuel','C2H5OH(L)','wt%',75,'t(k)',300,'fuel','H2O(L)','wt%',25,'t(k)',300,'oxid','O2(L)','wt%',100,'t(k)',90.0,'output','transport','mks','end');
%         C_star2(i,j) = RKT2(j).output.eql.cstar(1);
%         Isp2(i,j) = RKT2(j).output.eql.isp(3);
%         Tc2(i,j) = RKT2(j).output.eql.temperature(1);
%     end
%     j
% end
% 
% of_data = [1.5 1.48 0.85 1.2];
% cstar_data = [1456 1280 1050 1137];
% isp_data = [174 164 149 159];
% 
% figure(1)
% hold on
% colors = orderedcolors('gem');
% for i=1:length(Pc)
%     plot(O_F,C_star(:,i),'LineWidth',1.5,'Color',colors(i,:))
%     plot(O_F,C_star2(:,i),'LineWidth',1.5,'LineStyle','--','Color',colors(i,:))
% end
% plot(O_F,C_star2(:,1)*0.9,'LineWidth',1.0,'Color','k','LineStyle','--')
% plot(O_F,C_star2(:,1)*0.75,'LineWidth',1.0,'Color','k','LineStyle','--')
% plot(O_F,C_star(:,1)*0.9,'LineWidth',1.0,'Color','k','LineStyle','-')
% plot(O_F,C_star(:,1)*0.75,'LineWidth',1.0,'Color','k','LineStyle','-')
% scatter(of_data,cstar_data,50,'x','MarkerEdgeColor','k','LineWidth',1.75)
% scatter(2.0,1780*0.9,50,'MarkerEdgeColor','k','LineWidth',1.75)
% scatter(2.5,1777*0.9,50,'MarkerEdgeColor','k','LineWidth',1.75)
% hold off
% legend('Kerosene  400psi','Ethanol','Kerosene  500psi','Ethanol','Kerosene  5000psi','Ethanol','90% C*','80% C*')
% xlabel('O/F Ratio')
% ylabel('C* (m/s)')
% grid on
% 
% figure(2)
% hold on
% for i=1:length(Pc)
%     plot(O_F,Isp(:,i),'LineWidth',1.5,'Color',colors(i,:))
%     plot(O_F,Isp2(:,i),'LineWidth',1.5,'LineStyle','--','Color',colors(i,:))
% end
% plot(O_F,Isp2(:,1)*0.9,'LineWidth',1.0,'Color','k','LineStyle','--')
% plot(O_F,Isp2(:,1)*0.75,'LineWidth',1.0,'Color','k','LineStyle','--')
% plot(O_F,Isp(:,1)*0.9,'LineWidth',1.0,'Color','k','LineStyle','-')
% plot(O_F,Isp(:,1)*0.75,'LineWidth',1.0,'Color','k','LineStyle','-')
% scatter(of_data,isp_data,50,'x','MarkerEdgeColor','k','LineWidth',1.75)
% scatter(2.0,202.1,50,'MarkerEdgeColor','k','LineWidth',1.75)
% scatter(2.5,203.6,50,'MarkerEdgeColor','k','LineWidth',1.75)
% hold off
% legend('Kerosene  400psi','Ethanol','Kerosene  500psi','Ethanol','Kerosene  5000psi','Ethanol','90% C*','80% C*')
% xlabel('O/F Ratio')
% ylabel('Specific Impulse (s)')
% grid on
% 
% figure(3)
% hold on
% for i=1:length(Pc)
%     plot(O_F,Tc(:,i),'LineWidth',1.5,'Color',colors(i,:))
%     plot(O_F,Tc2(:,i),'LineWidth',1.5,'LineStyle','--','Color',colors(i,:))
% end
% hold off
% legend('Kerosene  400psi','Ethanol','Kerosene  500psi','Ethanol','Kerosene  5000psi','Ethanol')
% xlabel('O/F Ratio')
% ylabel('Chamber Temperature (K)')
% grid on

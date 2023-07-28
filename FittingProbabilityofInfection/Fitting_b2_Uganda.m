clear all
close all
clc



ViralData = [2 5.23  6.56  3.68  3.69  4.54 4.76]; %log Scale

Probabilty_Infection = [0 2/11  3/5   1/18  1/33  2/16 4/29];


lb = zeros(1,2);
ub = [10, 10];

k = [1e-5, 0.6];



[k,fval] =  fminsearchbnd(@err_in_data,k,lb,ub,optimset('Display','iter','TOLX', 1e-14, 'TOLFun', 1e-14));

a = k(1);
h = k(2);

V = 2:0.1:8;
Model_Prob_Infection = 1 - exp(-a*(V- 2).^h);

data_points = length(ViralData);
KK = 2;
AIC = data_points*log(fval/data_points) + 2*KK +(2*KK*(KK+1))/(data_points - KK -1)

figure
plot(ViralData, Probabilty_Infection,'r.','MarkerSize',25)
hold on
plot(V, Model_Prob_Infection, 'b','LineWidth',4)
ylim([0,1])
set(gca,'FontSize',15,'FontName','Arial','linewidth',3,'FontWeight','bold')
set(gca, 'YGrid', 'on', 'XGrid', 'off','LineWidth', 3,'fontsize',15)
title('Uganda Strain','FontSize',16,'FontName','Arial','FontWeight','bold')
xlabel('Viral Load (log scale)','FontSize',18,'FontName','Arial','FontWeight','bold')
ylabel('Probabilty of Infection','FontSize',18,'FontName','Arial','FontWeight','bold')


 function error_in_data = err_in_data(k) 

 ViralData = [2 5.23  6.56  3.68  3.69  4.54 4.76]; %log Scale

 Probabilty_Infection = [0 2/11  3/5   1/18  1/33  2/16 4/29];

 a = k(1);
 h = k(2);

 
 Model_Prbobability = 1 - exp(-a*(ViralData -2).^h);
 
 error_in_data = sum((Model_Prbobability - Probabilty_Infection).^2) ;           

 end

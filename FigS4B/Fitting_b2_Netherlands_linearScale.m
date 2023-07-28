clear all
close all
clc


ViralData =  [10^2  10^(2.7) 10^(3.04)  10^(3.49) 10^(3.66) 10^(5.99)  10^(7.18) ]; %linear Scale

Probabilty_Infection = [0  1/10  3/24  1/21  3/22  2/12  3/5];


lb = zeros(1,2);
ub = [10 10];

k = [1e-5, 0.6];



%[k,fval] =  fminsearchbnd(@err_in_data,k,lb,ub,optimset('Display','iter','TOLX', 1e-14, 'TOLFun', 1e-14));
[k,fval] = fmincon(@err_in_data, k,[], [], [], [], lb, ub,[],optimset('Display','iter','TOLX', 1e-14, 'TOLFun', 1e-14))
a = k(1);
h = k(2);


V = 100:1000:10^8;
Model_Prob_Infection = 1 - exp(-a*(V-100).^h);

data_points = length(ViralData);
KK = 2;
AIC = data_points*log(fval/data_points) + 2*KK +(2*KK*(KK+1))/(data_points - KK -1)

figure
semilogx(ViralData, Probabilty_Infection,'r.','MarkerSize',25)
hold on
semilogx(V, Model_Prob_Infection, 'b','LineWidth',4)
ylim([0,1])
set(gca,'FontSize',15,'FontName','Arial','linewidth',3,'FontWeight','bold')
set(gca, 'YGrid', 'on', 'XGrid', 'off','LineWidth', 3,'fontsize',15)
title('Netherlands Strain','FontSize',16,'FontName','Arial','FontWeight','bold')
xlabel('Viral Load (linear scale)','FontSize',18,'FontName','Arial','FontWeight','bold')
ylabel('Probabilty of Infection','FontSize',18,'FontName','Arial','FontWeight','bold')



 function error_in_data = err_in_data(k) 

ViralData =  [10^2  10^(2.7) 10^(3.04)  10^(3.49) 10^(3.66) 10^(5.99)  10^(7.18) ]; %linear Scale

Probabilty_Infection = [0  1/10  3/24  1/21  3/22  2/12  3/5];

a = k(1);
h = k(2);

 
 Model_Prbobability = 1 - exp(-a*(ViralData - 100).^h);
 
 error_in_data = sum((Model_Prbobability - Probabilty_Infection).^2) ;           

 end

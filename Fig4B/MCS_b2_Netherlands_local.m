clear all
close all
clc

   
numiter = 1000; 

     
true_params = [4.2927e-5,6.0582]; 

        
 X = zeros(length(true_params),numiter);




V_measure =  [2  2.7 3.04  3.49 3.66 5.99  7.18 ];
a = true_params(1);
h = true_params(2);


Model_Prob_Infection =  1 - exp(-a*(V_measure- 2).^h);
 

 noiselevel = 0.1;

 total_ARE =  zeros(length(noiselevel), length(true_params));
 actual_sigma = zeros(length(noiselevel),1);

 total_ARE_Table = {'a', 'h'};


for noisei = 1:1
    
rng default
noiselev = noiselevel(noisei)

    parfor i = 1:numiter
             i

  
Probabilty_Infection_Data = noiselev*randn(size(Model_Prob_Infection)) + Model_Prob_Infection;
sigma_d(i) = norm(Model_Prob_Infection - Probabilty_Infection_Data,2)/norm(Model_Prob_Infection,2)
 
            lb = [0; 0];
            ub = [10; 10];
             
             %k = fmincon(@(k)err_in_data(k,ViralData), true_params,[], [], [], [], lb, ub,[],optimset('MaxFunEvals', 1e+5,'MaxIter',1e+5,'TolX',1e-12, 'Display','iter'))

             k = fmincon(@(k)err_in_data(k,Probabilty_Infection_Data), true_params,[], [], [], [], lb, ub)
             % k =  fminsearchbnd(@(k)err_in_data(k,Probabilty_Infection_Data), true_params,lb,ub,optimset('TOLX', 1e-14, 'TOLFun', 1e-14));
             X(:,i) = k';
             
     end
        
        arescore = zeros(1,length(true_params));

    for i = 1:length(true_params)
        arescore(i) = 100*sum(abs(true_params(i) - X(i,:))/abs(true_params(i)))/numiter;
    end
    
    actual_sigma(noisei) = mean(sigma_d);
    total_ARE(noisei,:) = round(arescore,1);
    total_ARE_Table(noisei+1,:) = num2cell(total_ARE(noisei,:));

end

V = 2:0.1:8;
samples = zeros(length(V), 1000);
for i = 1:1000
    samples(:,i) = 1 - exp(-X(1,i)*(V- 2).^X(2,i));
end

plot(V, samples, 'Color', [0.9 1 0.8],'LineWidth',1) 
hold on
ViralData =  [2  2.7 3.04  3.49 3.66 5.99  7.18 ]; %log Scale
Probabilty_Infection = [0  1/10  3/24  1/21  3/22  2/12  3/5];
plot(ViralData, Probabilty_Infection,'r.','MarkerSize',25)
hold on
V = 2:0.1:8;
Model_Prob_Infection = 1 - exp(-4.2927e-5*(V - 2).^6.0582);
plot(V, Model_Prob_Infection, 'b','LineWidth',4)
ylim([0,1])
set(gca,'FontSize',12,'FontName','Arial','linewidth',3,'FontWeight','bold')
set(gca, 'YGrid', 'on', 'XGrid', 'off','LineWidth', 2,'fontsize',14)
title('Nertherlands Strain','FontSize',16,'FontName','Arial','FontWeight','bold')
xlabel('Viral Load (log scale)','FontSize',18,'FontName','Arial','FontWeight','bold')
ylabel('Probabilty of Infection','FontSize',18,'FontName','Arial','FontWeight','bold')
 
function error_in_data = err_in_data(k,Probabilty_Infection_Data) 


V_measure =  [2  2.7 3.04  3.49 3.66 5.99  7.18 ];
  
Model_Prbobability_fun = 1 - exp(-k(1)*(V_measure- 2).^k(2));

error_in_data = sum((Model_Prbobability_fun - Probabilty_Infection_Data).^2) ; 


end



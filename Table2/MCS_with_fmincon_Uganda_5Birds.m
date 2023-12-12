clear all
close all
clc

   
numiter = 1000; 

     
% true_params = [0.000717,8.13,41.22,43,26.94]; 

true_params = [1.4e-4,3.36,6.74,18.2,35.8]; 

        
 
X = zeros(length(true_params),numiter); 


dt = 0.1;

tforward = 0:dt:8;

t_v_measure =  [1, 2, 3, 4, 5, 6, 7]/dt + 1;
initial_cond = [4000000 0 0 10];


[~, y_trp] = ode23s(@(t,y)Model_USUTU_WithinHost(y,true_params),tforward,initial_cond);

 Model_Viral = log10(y_trp(t_v_measure(:),4));

 noiselevel = [0.01, 0.05, 0.1, 0.2];

 total_ARE =  zeros(length(noiselevel), length(true_params));


 total_ARE_Table = {'beta', 'd', 'delta',  'pi', 'c'};


for noisei = 1:4
    
rng default
noiselev = noiselevel(noisei)

    parfor i = 1:numiter
            i

  ViralData = (noiselev*randn(length(t_v_measure),1)) + Model_Viral;
  
 
            
            lb = [exp(-23) 0.5 0 0 0];  
            ub = [exp(2) 10 100 1e8 100];
             
             %k = fmincon(@(k)err_in_data(k,ViralData), true_params,[], [], [], [], lb, ub,[],optimset('MaxFunEvals', 1e+5,'MaxIter',1e+5,'TolX',1e-12, 'Display','iter'))

             k = fmincon(@(k)err_in_data(k,ViralData), true_params,[], [], [], [], lb, ub)
             X(:,i) = k';
             
     end
        
        arescore = zeros(1,length(true_params));

    for i = 1:length(true_params)
        arescore(i) = 100*sum(abs(true_params(i) - X(i,:))/abs(true_params(i)))/numiter;
    end
    
    total_ARE(noisei,:) = round(arescore,1);
    total_ARE_Table(noisei+1,:) = num2cell(total_ARE(noisei,:));

end

function error_in_data = err_in_data(k, ViralData) 
 
dt = 0.1;

tforward = 0:dt:8;

t_v_measure =  [1, 2, 3, 4, 5, 6, 7]/dt + 1;
initial_cond = [4000000 0 0 10];
 
 [~,y] = ode23s(@(t,y)Model_USUTU_WithinHost(y,k),tforward,initial_cond);
 

 
 Model_Viral_f = log10(y(t_v_measure(:),4));
 
 
 error_in_data = sum((Model_Viral_f - ViralData).^2);
                        

 end

function dy = Model_USUTU_WithinHost(y,k)

dy = zeros(4,1);

%params = [beta d delta pi c]
beta = k(1);
d = k(2);
delta = k(3);
pi = k(4);
c = k(5);



T = y(1);
E = y(2);
I = y(3);
V = y(4);


dy(1) = - beta* V.*T ;
dy(2) = beta* V.*T  - d*E;
dy(3) = d*E - delta*I;
dy(4) = pi*I - c*V;
 
end
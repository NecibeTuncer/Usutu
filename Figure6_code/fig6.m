%multiscale model output generated from fitted transmission model and
%with-inhost models using FDM 
%figure 6

clear all
clc

Final_Time = 100; %measured in days
Final_Infection_Age = 8; %measured in days
dt = 0.01; % dtau = dt;

initial_condition = [4e6 0 0 10]
infection_timeforward = 0:dt:Final_Infection_Age;

[t_wh, y_wh] = ode23s(@(t,y)Usutu_WithinHost(y),infection_timeforward,initial_condition);

%parameters;

V = y_wh(:,4);
V_log = log10(V);
indc = find(V_log <=2);%force all V below 2 to be at 2 (prob 0 in transmission model)
V_log(indc) = 2;

c=1;%contact rate
beta_b = 0.2; %bird infectivity
gb = 1/7;%bird recovery
mu_v = 1/60; %mosquito birth/death rate
Lambda_v = 1/60;
mu_b = 1/720;%bird birth/death rate
Lambda_b = 1/720;

%beta_v = c*(1 - exp(-0.1034*(V_log - 2 )));  %b1 transmission model Netherlands
beta_v = c*(1-exp(-0.0034*(V_log-2).^3.6697)); %b2 transmission model Usutu

tforward = 0:dt:Final_Time;

Sv = zeros(length(tforward),1);%initializing output storage
Iv = zeros(length(tforward),1);
Sb = zeros(length(tforward),1);
Ib = zeros(length(tforward),length(infection_timeforward));
Ibt = zeros(length(tforward),1);

Sv(1) = 0.95;%vector ICs
Iv(1) = 0.05;

Sb(1) = 0.99; %bird ICs
Ib(1,:) = 0.01/length(Ib(1,:));

Ibt(1) = dt*sum(Ib(1,:));%sum of infect birds over all age classes


for n=1:length(tforward)-1
    
    Int = dt*sum(beta_v'.*Ib(n,:));

    Sv(n+1) =  (Sv(n) + dt*Lambda_v)/(1 + dt*Int + dt*mu_v);%update step vectors
    Iv(n+1) = (dt*Sv(n+1)*Int + Iv(n))/(1 + dt*mu_v);

    Sb(n+1) = (Sb(n) +dt*Lambda_b)/(1 + dt*beta_b*Iv(n+1) + dt*mu_b);%update birds

    Ib(n+1, 2:end) = Ib(n,1:end-1)/( 1 + dt*gb+ dt*mu_b);%update Ib 
    Ib(n+1, 1) = beta_b*Sb(n+1)*Iv(n+1);%boundary cond

    Ibt(n+1) = dt*sum(Ib(n+1,:));%sum of infect birds over all age classes
end

figure(1)
subplot(2,2,2)
plot(tforward,Ibt,'r--','LineWidth',1)
hold on
ylabel('Infected Birds','FontSize',14)

subplot(2,2,4)
plot(tforward,Iv,'r--','LineWidth',1)
hold on
ylabel('Infected Mosquitos','FontSize',14)
xlabel('days','FontSize',14)


subplot(2,2,1)
plot(tforward,Sb,'b--','LineWidth',1)
hold on
ylabel('Susceptible Birds','FontSize',14)

subplot(2,2,3)
plot(tforward,Sv,'b--','LineWidth',1)
hold on
ylabel('Susceptible Mosquitos','FontSize',14)
xlabel('days','FontSize',14)

 function dy = Usutu_WithinHost(y) %inhost ode

dy = zeros(4,1);

 %true_params = [4.66e-5,7.07,6.95,7.49,48.8];  %Netherlands pop parms
 true_params = [1.4e-4,3.36,6.74,18.2,35.8]; %Usutu pop parms

beta = true_params(1);
d = true_params(2);
delta = true_params(3);
pi = true_params(4);
c = true_params(5);

T = y(1);
E = y(2);
I = y(3);
V = y(4);

dy(1) = - beta* V.*T ;
dy(2) = beta* V.*T  - d*E;
dy(3) = d*E - delta*I;
dy(4) = pi*I - c*V;

end
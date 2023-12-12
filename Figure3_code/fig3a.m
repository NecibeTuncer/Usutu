%monolix individual and population fits

load("netherlands_data.mat") %load netherlands birds

% load("uganda_data.mat") %load uganda birds
% all1=all(:,1:4);
% all1(:,5)=all(:,6);
% all=all1;

%ICs
T0=4e6; E0=0; I0=0; V0=10;
Y0=[T0 E0 I0 V0];

for i=1:10%solve ode with fit parms
    tdata=1:1:7;
    vdata=10.^all(:,i);
    figure(1)
    plot(tdata,(vdata),'ko')
    hold on
    plot(tdata,(vdata),'ko')
    hold on
end

parm_pop=[log(6.4e-5) 7.14 58.7 4.56 7.43];%median fit Netherlands
%parm_pop=[log(0.000292) 6.98 70.7 17.9 3.29];%median fit Uganda
sol=ode15s(@usuv_ode, [0 8], Y0, [], parm_pop);
% 
figure
plot(sol.x,(sol.y(4,:)),'Linewidth',1.5)
hold on






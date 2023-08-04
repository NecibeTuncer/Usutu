%monolix individual and population fits
%figures 2 and S1
%individual and pop parameters from monolix fitting

% load("netherlands_data.mat") %load netherlands birds
load("uganda_data.mat") %load uganda birds

%individual parameter estimates

% parm_table=[4.72E-05	7.05629	48.5969	6.51759	7.32115
% 1.48E-06	7.2646	20.0022	171.436	7.80424
% 1.88E-05	6.8829	68.855	17.232	6.95508
% 0.00191758	7.23835	188.436	0.40632	7.09338
% 0.000884623	7.13035	169.039	0.817999	7.62675
% 0.00124541	7.14016	178.665	0.560215	7.365
% 2.29E-06	7.38555	32.441	114.121	7.55569
% 0.00263763	7.06458	195.17	0.288624	7.10853
% 8.07E-05	7.02857	48.0758	2.6027	7.55352
% 1.92E-06	7.17623	30.6786	144.345	7.49968];%monolix parm fits
% % Netherlands

parm_table=[4.68E-04	7.72786	1274.09	11.6409	3.23048
3.93E-06	6.32984	0.941346	18.0595	3.28975
5.74E-04	7.40437	95.9084	17.2582	3.39439
4.81E-05	6.6322	15.3269	17.9154	3.88698
2.92E-04	6.98131	70.6545	18.0213	3.05369];%monolix parm fits
% Uganda
%remove bird 8
all1=all(:,1:4);all1(:,5)=all(:,6);all=all1;
%column 1=beta; column 2=delta; column 3=c; column 4=p; column 5=k


%ICs
T0=4e6; E0=0; I0=0; V0=10;
Y0=[T0 E0 I0 V0];
bird_count=length(parm_table(:,1));

for i=1:bird_count%solve ode with fit parms
    parms=[log(parm_table(i,1)) parm_table(i,2) parm_table(i,3) parm_table(i,4) parm_table(i,5)];
    sol_{i}=ode15s(@usuv_ode, [0 8], Y0, ['RelTol', 1e-9, 'AbsTol', 1e-12], parms);
    tdata=1:1:7;
    vdata=10.^all(:,i);
    figure(1) %individual fits
    subplot(4,3,i)
    semilogy(sol_{i}.x,sol_{i}.y(4,:),'Linewidth',1.5)
    xlim([0,7])
    ylim([1,10^8])
    hold on
    semilogy(tdata,vdata,'ko')
    hold on
    semilogy(tdata,vdata,'ko')
    hold on
    figure(2) %population plot
    semilogy(tdata,vdata,'ko')
    hold on
    semilogy(tdata,vdata,'ko')
    hold on
end

%Population fit
% parm_pop=[log(4.66E-05)	6.945069905	48.8153521	7.486610076	7.072362766];%Netherlands
parm_pop=[log(1.4E-04)	6.742574241	35.79268977	18.19087424	3.356925808];%Usutu
sol=ode15s(@usuv_ode, [0 8], Y0, ['RelTol', 1e-9, 'AbsTol', 1e-12], parm_pop);

figure(2)
semilogy(sol.x,sol.y(4,:),'Linewidth',1.5)
xlim([0,8])
ylim([1,10^8])


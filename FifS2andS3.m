%monolix scatter plot of parameter sampling
%figures S2 and S3

load("scatter_data_netherlands_monolix.mat") %netherlands
% load("scatter_data_uganda_5_monolix.mat") %uganda

figure(2)
%beta vs
subplot(5,5,2)
plot(scatter_data(:,3),scatter_data(:,1),'ko','MarkerFaceColor','k')
hold on
ylabel('\beta')
xlabel('c')
subplot(5,5,3)
plot(scatter_data(:,2),scatter_data(:,1),'ko','MarkerFaceColor','k')
hold on
ylabel('\beta')
xlabel('\delta')
subplot(5,5,4)
plot(scatter_data(:,5),scatter_data(:,1),'ko','MarkerFaceColor','k')
hold on
ylabel('\beta')
xlabel('k')
subplot(5,5,5)
plot(scatter_data(:,4),scatter_data(:,1),'ko','MarkerFaceColor','k')
hold on
ylabel('\beta')
xlabel('p')

%c vs
subplot(5,5,6)
plot(scatter_data(:,1),scatter_data(:,3),'ko','MarkerFaceColor','k')
hold on
ylabel('c')
xlabel('\beta')
subplot(5,5,8)
plot(scatter_data(:,2),scatter_data(:,3),'ko','MarkerFaceColor','k')
hold on
ylabel('c')
xlabel('\delta')
subplot(5,5,9)
plot(scatter_data(:,5),scatter_data(:,3),'ko','MarkerFaceColor','k')
hold on
ylabel('c')
xlabel('k')
subplot(5,5,10)
plot(scatter_data(:,4),scatter_data(:,3),'ko','MarkerFaceColor','k')
hold on
ylabel('c')
xlabel('p')

%delta vs
subplot(5,5,11)
plot(scatter_data(:,1),scatter_data(:,2),'ko','MarkerFaceColor','k')
hold on
ylabel('\delta')
xlabel('\beta')
subplot(5,5,12)
plot(scatter_data(:,3),scatter_data(:,2),'ko','MarkerFaceColor','k')
hold on
ylabel('\delta')
xlabel('c')
subplot(5,5,14)
plot(scatter_data(:,5),scatter_data(:,2),'ko','MarkerFaceColor','k')
hold on
ylabel('\delta')
xlabel('k')
subplot(5,5,15)
plot(scatter_data(:,4),scatter_data(:,2),'ko','MarkerFaceColor','k')
hold on
ylabel('\delta')
xlabel('p')

%k vs
subplot(5,5,16)
plot(scatter_data(:,1),scatter_data(:,5),'ko','MarkerFaceColor','k')
hold on
ylabel('k')
xlabel('\beta')
subplot(5,5,17)
plot(scatter_data(:,3),scatter_data(:,5),'ko','MarkerFaceColor','k')
hold on
ylabel('k')
xlabel('c')
subplot(5,5,18)
plot(scatter_data(:,2),scatter_data(:,5),'ko','MarkerFaceColor','k')
hold on
ylabel('k')
xlabel('\delta')
subplot(5,5,20)
plot(scatter_data(:,4),scatter_data(:,5),'ko','MarkerFaceColor','k')
hold on
ylabel('k')
xlabel('p')

%p vs
subplot(5,5,21)
plot(scatter_data(:,1),scatter_data(:,4),'ko','MarkerFaceColor','k')
hold on
ylabel('p')
xlabel('\beta')
subplot(5,5,22)
plot(scatter_data(:,3),scatter_data(:,4),'ko','MarkerFaceColor','k')
hold on
ylabel('p')
xlabel('c')
subplot(5,5,23)
plot(scatter_data(:,2),scatter_data(:,4),'ko','MarkerFaceColor','k')
hold on
ylabel('p')
xlabel('\delta')
subplot(5,5,24)
plot(scatter_data(:,5),scatter_data(:,4),'ko','MarkerFaceColor','k')
hold on
ylabel('p')
xlabel('k')


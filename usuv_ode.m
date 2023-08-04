%Usutu TIEV ODE
function dy=usuv_ode(t,y,parms)

beta=exp(parms(1));%parameters
delta=parms(2);
c=parms(3);
p=parms(4);
k=parms(5);

dy=zeros(4,1);

dy(1)=-beta*y(1)*y(4);
dy(2)=beta*y(1)*y(4)-k*y(2);
dy(3)=k*y(2)-delta*y(3);
dy(4)=p*y(3)-c*y(4);

if y(4)<=1e-4 %stop system once virus is gone
   dy(4)=0;
end

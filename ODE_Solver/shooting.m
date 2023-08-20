clc;
clear;
close all;
global m
%%
tic
h=35;
k=250;
Ro=0.02;
R1=0.1;
t=0.001;
To=100+273.15;
Tamb=25+273.15;
Tol=1e-5;
imax=1000;
m=(2*h/(k*t))^0.5;
%% B.C at r=Ro
Tetao=To-Tamb;
%% 2 initial guesses for dT at r=Ro
dTo(1)=-0.1;
dTo(2)=-1;
%%  solution
for i=1:250
    [r,y]=ode45(@odefun ,[Ro R1],[Tetao dTo(i)]);
    T=y(:,1)+Tamb;
    dT=y(:,2);

    %Dtermination of the Robins condition at r=R1
f(i)=dT(end)+(h/k)*(T(end)-Tamb);
if i>1
    if abs(f(i))>Tol
        dTo(i+1)=dTo(i)-f(i)*(dTo(i)-dTo(i-1))/(f(i)-f(i-1));
    else
        break
    end
end
end
%%  plotting results
plot(r*100,T-273.15,'r','LineWidth',2)
xlabel('r(cm)')
ylabel('T(\circC)')
toc








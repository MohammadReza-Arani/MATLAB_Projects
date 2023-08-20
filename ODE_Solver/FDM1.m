clc;
clear;
close all;
%%
tic
tinf=25+273.15;   %K
to=100+273.15;   %K
h=35;  %w/m^2.K
k=250;  %w/m.K
b0=0.001;  %m
Ri=0.02;  %m
Ro=0.1;  %m
m=sqrt(2*h/(b0*k));
n=1000;
h2=(Ro-Ri)/(n-1);
r=Ri:h2:Ro;
%%
for i=2:n
    %B.C
    A(1,1)=-2-(m^2)*h2^2;
    A(1,2)=(2*r(1)+h2)/(2*r(1));
    b(1)=((-2*r(1)+h2)/(2*r(1)))*(to-tinf);
    if i<n
        A(i,i)=-2-(m^2)*h2^2;
        A(i,i-1)=(2*r(i)-h2)/(2*r(i));
        A(i,i+1)=((2*r(i)+h2)/(2*r(i)));
        b(i)=0;
    elseif i==n
        A(i,i)=(1+(h*h2)/k);
        A(i,i-1)=-1;
        b(i)=0;
    end
end
b2=b';
x=A \ b2;
t=x+tinf;
rmm=r*100;
comet(rmm,t-273.15);
ho=plot(rmm,t-273.15,'b','LineWidth',2);
hold on;
xlabel('r(cm)')
ylabel('T(\circC)')
axis light
toc








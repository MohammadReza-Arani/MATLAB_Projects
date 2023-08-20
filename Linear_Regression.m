%% Linear Regression
clear; clc;
e0=8.85*10^(-12);
mu0=4*pi*10^(-7);
x=[-1 0 2];
y=[2 3 4 ];


[xx]=sum(x.^2,'all');
[x_s]=sum(x,'all');
[y_s]=sum(y,'all');
[xy]=sum(x.*y,'all');
N=length(x);
a=(N*xy-x_s*y_s)/(N*xx-x_s*x_s);
b=(y_s*xx-x_s*xy)/(N*xx-x_s*x_s);% if it passes through O(o,o) b is equal to zero;
a2=xy/xx;
dd=sum(abs(a2*x-y));
figure()
plot(x,y);
title("Linear Regression for Given Data Points ");
xlabel("dI/dt(A/s)");
ylabel("U(mv)(emf)");
hold on
xtemp=min(x):0.1:max(x);
plot(xtemp,a*xtemp+b,'black:');
hold on
plot(xtemp,a2*xtemp,'g');
hold on
for i=1:length(x)
    hold on
    plot(x(i),y(i),'r*');
end
legend("Original form","Linear Regression with b~=0","Linear Regression with b=0",'Given Points');
grid on

clear;
clc;

u = @( t ) double( t >= 0) ;
nx=0:14;  p = linspace(-5,30,length(nx));
x=(1/2).^(p-2).*u(p-2);

nh=0:24; p = linspace(-5,30,length(nh));
h=u(p);

ny=convindicesme(nx,nh);

y=conv(x,h);
title('plot of signalyn=xn*xn','interpreter','latex','fontsize',16);
xlabel('$n$','interpreter','latex');
ylabel('$y[n]$','interpreter','latex');
stem(ny,y,'linewidth',2);
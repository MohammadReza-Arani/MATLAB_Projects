close all
Fc=300e6;
C=3e8;
Landa=C/Fc;
K=2*pi/Landa;
dmax=10;
d=(0:Landa*0.5:dmax)';
M=length(d);

Theta=-90:0.1:90;
Theta_len=length(Theta);
MAP=exp(1j*K*d*sind(Theta));

zav=45;
a=exp(1j*K*d*sind(zav));
g=abs(a'*MAP);
g=g/max(g);
plot(Theta,(g),'b')
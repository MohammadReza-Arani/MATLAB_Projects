clear;
clc;
f0=10;
fs=1000;
fs1=25;
fs2=50;
t=linspace(-10/f0,10/f0,2*fs+1);
x=f0*(sinc(f0*t)).^2+15*sin(2*pi*f0*t);
% figure
% plot(t,x,'DisplayName','x(t)');
% xlabel('T');


x1=x(1:fs/fs1:end);
x2=x(1:fs/fs2:end);

% figure
% subplot(2,1,1);
% stem(x1);
% subplot(2,1,2);
% stem(x2);

f01=16000;
fs1=16000;
fs2=32000;
fs3=64000;

t1=-1000/f01:1/fs1:1000/f01;
t2=-1000/f01:1/fs2:1000/f01;
t3=-1000/f01:1/fs3:1000/f01;

x1=f01*(sinc(f01*t1)).^2+15*sin(2*pi*f01*t1);
x2=f01*(sinc(f01*t2)).^2+15*sin(2*pi*f01*t2);
x3=f01*(sinc(f01*t3)).^2+15*sin(2*pi*f01*t3);


%  figure
%  subplot(2,2,1);
% stem(t1,x1);
% xlim([-10/f01,10/f01]);
%  subplot(2,2,2);
% stem(t2,x2);
% xlim([-10/f01,10/f01]);
%  subplot(2,2,3);
% stem(t3,x3);
% xlim([-10/f01,10/f01]);

Fourierfunc(1,x1,fs1,'fs=16000');
Fourierfunc(1,x2,fs2,'fs=32000');
Fourierfunc(1,x3,fs3,'fs=64000');









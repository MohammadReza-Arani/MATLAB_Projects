clc;
clear;
ecg_data=load('ecg.dat');
fs=500;
t=0:1/fs:4;

ecg_short_data=ecg_data(1:4*fs+1);


% figure
% plot(t,ecg_short_data,'DisplayName','ecg_short_data');
% xlabel('T');
% ylabel('ecg short data');
% grid on
% legend('show');
% % %  
% % 
% 
%  Fourierfunc(2,ecg_short_data,fs,'ecg short data');
 
ecg_highpassed=filter([1 0 -1],[fs/2],ecg_short_data);
% figure
% subplot(2,1,1);
% plot(t,fs*ecg_highpassed);
% title('ecg highpassed');
% legend('show');
% subplot(2,1,2);
% plot(t,ecg_short_data);
% title('ecg short data');
% xlabel('T');
% ylabel('ecg');
% grid on
% legend('show')

% Fourierfunc(2,ecg_highpassed,fs,'ecg highpass');

fc=100;
wc=fc/(fs/2);
fc2=70;
wc2=fc2/(fs/2);
n1=2;
n2=4;
n3=8;
n4=8;

[b1 a1]= butter(n1,wc,'low');
[b2 a2]= butter(n2,wc,'low');
[b3 a3]= butter(n3,wc,'low');

[b4 a4]= butter(n4,wc2,'low');

ecg_lowpassed1=filter(b1,a1,ecg_highpassed);
ecg_lowpassed2=filter(b2,a2,ecg_highpassed);
ecg_lowpassed3=filter(b3,a3,ecg_highpassed);

ecg_lowpassed=filter(b4,a4,ecg_highpassed);

s=0:1/10:1;
H1=sqrt(1./(1+(s/(j*wc)).^(2*n1)));
H2=sqrt(1./(1+(s/(j*wc)).^(2*n2)));
H3=sqrt(1./(1+(s/(j*wc)).^(2*n3)));

% figure
% subplot(3,1,1)
% plot(s,H1,'DisplayName','N=2');
% subplot(3,1,2)
% plot(s,H2,'DisplayName','N=4');
% subplot(3,1,3)
% plot(s,H3,'DisplayName','N=8');

% Fourierfunc(2,ecg_lowpassed1,fs,'ecg_lowpassed1,N=2');
% 
% Fourierfunc(2,ecg_lowpassed2,fs,'ecg_lowpassed2,N=4');
% 
% Fourierfunc(2,ecg_lowpassed3,fs,'ecg_lowpassed3,N=8');


fpl = 50;
b = [1-2*cos(2*pi*fpl/fs) 1];
a = [1];
ecg_notched = filter(b,a,ecg_lowpassed);
figure;
plot(t,ecg_notched);
title('x(t) band-pass filter');
xlabel('T');
ylabel('amp');
grid on
% 
% [cor,lag] = xcorr(ecg_notched);
% figure;
% plot(lag,cor);
% title('Auto-Correlation Function of Filtered ECG');
% xlabel('Delay (ms)');






clear;
clc;
T=0.1;
Fs=1/T;
t= -2:1/Fs:2;
X1=heaviside(t)-heaviside(t-1);
Y1=fft(X1);
P1=abs(Y1);
L=length(X1);
f=(0:L-1)*(Fs/L); %frequency range
A1=unwrap(angle(Y1));
figure
subplot(2,1,1)
stem(f,P1,'LineStyle','-.', 'MarkerFaceColor','yellow','MarkerEdgeColor','blue');
xlabel('frequency');
ylabel('amp');
title('fft(X1(t))');
legend('fft(x1)');
grid on
subplot(2,1,2)
stem(f,A1,'LineStyle','-.', 'MarkerFaceColor','red','MarkerEdgeColor','green');
xlabel('frequency');
ylabel('phase');
title('fft(X1(t))');
legend('fft(x1)');
grid on

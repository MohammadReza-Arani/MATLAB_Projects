clear;
clc;
T=0.1;
Fs=1/T;
t= -2:1/Fs:2;
X1=heaviside(t)-heaviside(t-1);
Y1=fft(X1);
L=length(X1);
f=(0:L-1)*(Fs/L); %frequency range
P1=abs(Y1);
A1=unwrap(angle(Y1));
figure
subplot(2,1,1)
stem(f,P1);
xlabel('frequency');
ylabel('amp');
title('fft(X1(t))');
legend('fft(x1)');
grid on
subplot(2,1,2)
stem(f,A1);
xlabel('frequency');
ylabel('phase');
title('fft(X1(t))');
legend('fft(x1)');
grid on
Y2=fftshift(Y1);
A2=unwrap(angle(Y2));
P2=abs(Y2);
L2=linspace(-Fs/2,+Fs/2,size(Y2,2));
figure
subplot(2,1,1)
plot(L2,P2);
xlabel('frequency');
ylabel('amp');
title('fft(X1(t))');
legend('fft(x1)');
grid on
subplot(2,1,2)
plot(L2,A2);
xlabel('frequency');
ylabel('phase');
title('fft(X1(t))');
legend('fft(x1)');
grid on





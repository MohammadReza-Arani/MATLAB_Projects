clear;
clc;
T=0.1;
Fs=1/T;
t=-15:T:15;
F1=heaviside(t)-heaviside(t-4);
F2=heaviside(t-3)-heaviside(t-7);
F3=(cos(t*pi-pi/3)).^3;
[A1,f1]=my_fft(F1, Fs);
[f1]=fftshift(f1);

figure
 
subplot(2,2,1);
plot(f1,A1,'DisplayName','fft_F1');
xlabel('frequency');
ylabel('amplitude');
title('fft(X1(t))');
legend('fft(F1)');
grid on

[A2,f2]=my_fft(F2, Fs);
[f2]=fftshift(f2);

subplot(2,2,2);
plot(f2,A2,'DisplayName','fft_F2');
xlabel('frequency');
ylabel('amplitude');
title('fft(X2(t))');
legend('fft(F2)');
grid on

[A3,f3]=my_fft(F3, Fs);
[f3]=fftshift(f3);

subplot(2,2,3);
plot(f3,A3,'DisplayName','fft_F3');
xlabel('frequency');
ylabel('amplitude');
title('fft(X3(t))');
legend('fft(32)');
grid on



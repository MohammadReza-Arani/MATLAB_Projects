%% Hw-2 BSS Not in Live Script
clear; clc; close all;
%% Part-A:


fc = 150e+06; % 150 MHz

Target_1_theta = 10;
Target_2_theta = 20;

f1 = 20e03; % 20KHz
f2 = 10e03; % 10KHz



Trecord = 1e-03; % 1ms recording time
delta_f = 1/Trecord; % Freq Resolution

fs = 1e06; % 1MHz
ts = 1/fs;

t = 0 : ts : Trecord-ts;
freq = -fs/2 : delta_f : fs/2-delta_f;


c = 3e8;
Lambda = c/fc;
k = 2*pi*fc/c ; % 2pi/lambda = 2pif/c

M = 10;
dz = 1;
h = 1; % height of first element from Ground
D = h + dz*(0:M-1);

Tau  = 10e-6; % 10us
Tau_num = round(Tau/ts);

PRI  = 0.1e-3; % 0.1ms
PRF = 1/PRI;
PRI_num = round(PRI/ts);
pulse_num = round(Trecord/PRI); % Kolle Pulse ma chanta PRI hast!
sample_num = round(Trecord/ts); % Kolle Pulse ma chanta nemoone hast

s1 = exp(1j*2*pi*f1*t);
s2 = exp(1j*2*pi*f2*t);

SL_PRI=[ones(1,Tau_num) zeros(1,PRI_num-Tau_num)];
sl=repmat(SL_PRI,1,pulse_num); % Repeat The Baseband SIgnal for each PRI in total Pulse Length

% Target 1:
theta_1= 10*pi/180;
fd1 = 2e3; % 2KHz
R1 = 6e3; % 6Km


% Target 2:
theta_2  = 20*pi/180;
fd2 = 1e3; % 1KHz
R2 = 6e3; % 6Km



% Noise Generation:  Independent from Signals! -- Gaussian --- Zero Mean:
Noise = randn(M, length(t)) + 1j*randn(M, length(t));  % M*T Noise Matrix for M elements of Antenna and T samples
Noise = Noise/sqrt(2);


Y =  exp(1j*k*D'*sin(theta_1) ) .* circshift(sl , floor(2*R1/c) ) .* exp(1j*2*pi*fd1*t)  + ... % Target 1
     exp(1j*k*D'*sin(theta_2) ) .* circshift(sl , floor(2*R2/c) ) .* exp(1j*2*pi*fd2*t);


Y_noisy = Y + Noise;
%%
figure(1)
subplot(2,1,1)
stem(t,sl)
xlabel('Time(seconds)')
grid on
subplot(2,1,2)
slf=fftshift(fft(sl));
stem(freq,abs(slf));
xlabel('Frequency(Hz)');
grid on


figure(2)
y_in=Y(randi(10,1), :); % Choose a random Terminal to plot its received Signal
y_in_f=fftshift(fft(y_in));


subplot(2,1,1)
stem(t,abs(y_in))
xlabel('Time(seconds)')
ylabel("Signal Amplitude")
grid on

subplot(2,1,2)
stem(freq,abs(y_in_f))
xlabel('Frequency(Hz)')
ylabel("Signal Amplitude")
grid on

%%  BeamForming
[U, S, V] = svd(Y_noisy); % SVD Decomposition!
disp(S(:,1:20)) % M*T Matrix! --> Diagnal Elements have Values


%% 
K = 2;  % Number of Targets   ----  Known but chosen wisely from SVD result 

Usig = U(:, 1:K);  % M*K
Unull = U(:, K+1:end); % M* (M-K)

Vsig = V(:, 1:K); % T * K
Vnull = V(:, K+1:end); % T * (M-K)


%%  Beamforming
theta = 0 : 0.5e-1 : 90;
a = exp(1j*k*D'*sind(theta)); % Steering Vector
g = sum(abs(a'*Usig).^2, 2); % g(theta)

figure(4)
plot(theta, g)
xlabel('\theta(degree)')
ylabel('g(\theta)')
grid on
title("Beamforming Objective Function")
[Beamforming_peaks , Beamformingidx ] = findpeaks(g);
hold on
plot(theta(Beamformingidx),Beamforming_peaks,"r*")

%% MUISC:

theta = 0 : 0.5e-1 : 90;
a = exp(1j*k*D'*sind(theta)); % Steering Vector
f = 1./sum(  abs(a'*Unull).^2 ,2 ); % f(theta)

figure(5)
close all;
plot(theta, f)
xlabel('\theta(degree)')
ylabel('g(\theta)')
grid on
title("MUSIC Objective Function")
[Music_peaks , Music_idx ] = findpeaks(f);
hold on
plot(theta(Music_idx),Music_peaks,"r*")

%% Frequency Detection: Beamforming

fd = (-4e3 : 1 : 4e3)'; % The Value to be found is the doppler of the targets!
s = circshift(sl , floor(2*R1/c) ).* exp(1j*2*pi*fd*t); % Signal Space with given Range  
g = sum(  abs(s*Vsig).^2 ,2 ); % g(fd)

figure(6)
close all;
plot(fd , g)
xlabel('fd(Hz)')
ylabel('g(\theta)')
grid on
title("Beamforming Objective Function")
[Beamforming_peaks , Beamformingidx ] = findpeaks(g);
hold on
plot(fd(Beamformingidx),Beamforming_peaks,"r*")
%% Frequency Detection MUSIC
fd = (-4e3 : 1 : 4e3)'; % The Value to be found is the doppler of the targets!
s = circshift(sl , floor(2*R1/c) ).* exp(1j*2*pi*fd*t); % Signal Space with given Range  
f = 1./sum(  abs(s*Vnull).^2 ,2 ); % f(fd)

figure(6)
close all;
plot(fd , f)
xlabel('fd(Hz)')
ylabel('f(fd)')
grid on
title("MUSIC Objective Function")
[Music_peaks , Music_idx ] = findpeaks(f);
hold on
plot(fd(Music_idx),Music_peaks,"r*")

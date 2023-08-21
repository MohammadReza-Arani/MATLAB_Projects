clear; clc; close all;
% Generate the input signal in the Z-transform domain
    num = [1 2];
    den = [1 -1/2];
    Xz = filter(num, den, [1 zeros(1,99)]); % Z-transform of input signal
    Xn = ifft(Xz); % Time-domain signal

% Find the system response
% b = [0.45 0.4 -1];
% a = [1 0.4 -0.45];



    yic=[0,3];
    xic=[2,2];
    b =[1,-0.4,-0.45];
    a=[0.45,0.4,-1];
    
    zi = filtic(b, a, yic,xic);
    yzi =filter(b,a,Xn,zi);
%% Second Method

    [r, p, k] = residuez(b, a);
    hn = impz(b, a, 100);
    yn = conv(hn, Xn);
    yn = yn(1:100);

% Plot the input signal and system response
    n = 1:100;
    subplot(2,1,1);
    stem(n, Xn);
    title('Input Signal');
    xlabel('n');
    ylabel('x[n]');
    subplot(2,1,2);
    stem(n, abs(yn));
    title('System Response');
    xlabel('n');
    ylabel('y[n]');


%% 
    figure()
    stem(n, abs(yzi));
    title('System Response');
    xlabel('n');
    ylabel('y[n]');
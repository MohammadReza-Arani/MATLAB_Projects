% A CONVOLUTION COMPUTING CODE IN MATLAB WITHOUT USING MATLAB BUILTIN FUNCTION conv(x,h)
close all;
clear;clc;

x = [0 0 0 1 0 0 0 1 -1 2 12 8 6 5 3.5 2 1 0 0 ];
h = [1 2 4 8 16 32 64 128 256 512 1024];
Y=myCun(x,h);

figure()
subplot(2,2,1)
stem(Y);
ylabel('Y[n]');
xlabel('----->n');
title('Convolution of Two Signals without conv function');
grid on

subplot(2,2,2)
Y_matlab = conv(x,h);
stem(Y_matlab)
ylabel('Y[n]');
xlabel('----->n');
title('Convolution of Two Signals With MATLAB conv function');
grid on


subplot(2,2,3)
stem(x)
ylabel('x[n]');
xlabel('----->n');
title('Input Signal');
grid on
subplot(2,2,4)
stem(h)
ylabel('h[n]');
xlabel('----->n');
title('System Impulse Response');
grid on
%% Functions
function Y=myCun(x,h)
    m=length(x);
    n=length(h);
    X=[x,zeros(1,n)]; 
    H=[h,zeros(1,m)]; 
    for i=1:n+m-1
    Y(i)=0;
        for j=1:m
            if(i-j+1>0)
                Y(i)=Y(i)+X(j)*H(i-j+1);
            else
            end
        end
    end
end
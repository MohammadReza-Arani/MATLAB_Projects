%% noghte miani
clear; clc;
    a=0;
    b=2;
    h=2*sqrt(2)*10^(-3);
    % x0=a:h/5:b;
    L=ceil((b-a)/h);
    y=zeros(L,1);
    syms x
    f(x)=2*x/(1+3*x^2);
    for i=1:L
       y(i)=double(f(a+(i-1)*h+h/2));
    end
    Result=h*sum(y,'all');
    disp(Result);
    disp("MATLAB Symbolic SOlution:")
    disp(double(int(f(x),x,a,b)))
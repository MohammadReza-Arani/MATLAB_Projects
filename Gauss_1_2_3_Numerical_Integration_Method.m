%% Gauss
clear;clc;
    a=0;
    b=pi/2;
    
    syms x_sym
    f(x_sym)=sin(x_sym);
    g(x_sym)=(b-a)/2*f(1/2*((b-a)*x_sym+(b-a)));
    %1 Point
        x=0;
        y1=2*double(g(x));
        disp("1-Point Integration : "+y1)
    %2 Point
        x=sqrt(3)/3;
        y2=double(g(x))+double(g(-x));
        disp("2-Point Integration : "+y2)
    % 3Point
        x1=-sqrt(3/5);
        x2=0;
        x3=-x1;
        y3=1/9*(5*double(g(x1))+8*double(g(x2))+5*double(g(x3)));
        disp("3-Point Integration : "+y3);

        disp("Symbolic Integration : "+double(int(f(x_sym),x_sym,a,b)))
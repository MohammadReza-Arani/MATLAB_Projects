%% Zoozanagheyi
clear; clc;
    a=0;
    b=1;
    h=0.03;
    % x0=a:h/5:b;
    L=ceil((b-a)/h);
    y=zeros(L+1,1);
    syms x
    f(x)=exp(-2*x)/(2+sin(2*x));
    for i=1:L+1
       if(i==1) 
       y(i)=1/2*double(f(a+(i-1)*h));
       elseif(i==L+1)
           y(i)=1/2*double(f(a+(i-1)*h));
       else
           y(i)=double(f(a+(i-1)*h));
       end
    end
    Result=h*sum(y,'all')
    disp("Symbolic Integration Result: ")
    double(int(f(x),x,0,1))
clear;clc;
%Newton-Raphson Method: Find Roots of a function
    syms x
    f(x)=sin(x)-x/2;
    D(x)=diff(f(x));
    format long
    
    x0=1.75;
    i=0;
    flag=1;
    threshold=10^(-5);
    while(flag)
    
    X1=x0-double(f(x0)/D(x0));
    
    x0=X1;
    % va dobare hamin ravand
    if(abs(double(f(X1)))<threshold)
        break;
    elseif(i>10^7)
         break;
    end
    i=i+1;
    
    
    end
    
    a=x0;
    disp("x equals to :")
    disp(X1);
    disp("With number of repetitions : ")
    disp(i);
    F=matlabFunction(f);
    disp("Matlab solution")
    double(fzero(f,a))


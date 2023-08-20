clear; clc;
%Bisection Method for Root Finding:
    a=0.25;
    b=0.27;
    syms x
    f(x)=3*x-exp(-x);
    
    flag=1;
    threshold=10^(-5);
    while(flag)
    X1=(a+b)/2;
    if(f(a)*f(X1)<0)
        b=X1;
    else
        a=X1;
    end
    % va dobare hamin ravand
    if(abs(double(f(X1)))<threshold)
        break;
    elseif(i>10^7)
         break;
    end
    i=i+1;
    
    end
    
    
    disp("x equals to :")
    disp(X1);
    disp("With number of repetitions : ")
    disp(i);
    F=matlabFunction(f);
    disp("Matlab solution")
    fzero(f,a)

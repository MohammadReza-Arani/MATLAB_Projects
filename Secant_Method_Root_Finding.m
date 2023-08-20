%% Secant method
% Same as Newton-Raphson but we use numerical Differentiation:
clear; clc;
    x0=-2;
    x1=1;
    syms x
    f(x)=x^3-x+1;
    D(x)=diff(f(x),x);
    flag=1;
    threshold=10^(-3);
    i=0;
    while(flag)
        if(double(f(x1)-f(x0))~=0)
            x2=x1-double(f(x1))*(x1-x0)/double(f(x1)-f(x0));% if faced with 0/0 it should be derivated to remove ambiguity
        else
            x2=x1-double(f(x1))/double(D(x1));   
        end
        
        if(abs(double(f(x2)))<threshold)
            flag=0;
        elseif(i>10^7)
            flag=0;    
        end
        i=i+1;
        x0=x1;
        x1=x2;
    end

    disp(x2);
    % MATLAB roots solution:
    disp("MATLAB Solution: ")
    disp(roots([1 0 -1 1]))
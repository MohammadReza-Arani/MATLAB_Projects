%% Sympson
clear; clc;
    a=0;
    b=4;
    n=4;
    
    h=((b-a)/n);
    y=zeros(n,1);
    syms x
    f(x)=exp(x);
    for i=1:n+1
        if((i==1)||(i==n+1))
            y(i)= f(a+(i-1)*h);
        elseif(mod(i,2)==0)
            y(i)= 4*f(a+(i-1)*h);
        else
            y(i)= 2*f(a+(i-1)*h);
        end
    end
    Result=h/3*sum(y,'all');
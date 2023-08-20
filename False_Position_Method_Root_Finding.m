%% False Position Method for Root Finding:
%
% nabejayi shabih hamoon bisection faghat migim aval ye hads mizanim bad az
%ebteda o entehaye baze ye khat migzare un khat ehtemallan 0 o ghat mikone
    clear; clc;
    a=0;
    b=1;
    i=0;
    syms x
    f(x)=2*x^3-1+sin(2*x);
    flag=1;
    threshold=10^(-3);
    while(flag)
    X1=double((a*f(b)-b*f(a))/(f(b)-f(a)));
    if(double(f(a)*f(X1))<0)
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
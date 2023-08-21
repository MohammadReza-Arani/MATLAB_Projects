%% Matlab Code:
% Part 1:
% x = 0:0.1:10;
% f = x.^2+sqrt(x)+sin(x)+cos(x)+tan(x)+cot(x)+log(x)+exp(x);

% Using Secant Method: Getting 
% 1) Zarible Moteghayer
% 2) Hadse Avalie
% 3) Mizane Khataye Nesbi

% Gives:
% 1) The root
% 2) Num of Iterations

function [Root,Num_iter]=Secant_custom(x0,Purt,threshold)
% x0=-2;
% x1=1;


    syms x
    f(x) = x^2+sqrt(x)+sin(x)+cos(x)+tan(x)+cot(x)+log(x)+exp(x);
    D(x)=diff(f(x),x);
    flag=1;
    % threshold=10^(-3);
    i=0;
    x1= x0+Purt;
    while(flag)
        if(double(f(x1)-f(x0))~=0)
            x2=x1-(x1-x0)*double(f(x1))/(double(f(x1))- double(f(x0)) );% age 0/0 shod bayad rafe ebham she ke hamoon moshtaghe
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
    Root = x2; Num_iter =i;
end
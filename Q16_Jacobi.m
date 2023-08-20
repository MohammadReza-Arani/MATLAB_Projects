%% Solution of x in Ax=b using Jacobi Method
clear;
clc;

% A=[4 -1 1;4 -8 1;-2 1 5];
% B=[7;-21;15];
A=[44 -1 1 2 12; 4 68 1 3 -2;   -2 1 35 12 0;13 1 -1 33 5; 1 9 6 -2 53];
B=[7;-21;15;8;22];

x=[0 0 0 0 0]';
n=size(x,1);
normVal=Inf; 
tic
tol=1e-2; itr=0;

while normVal>tol
    xold=x;
    
    for i=1:n
        sigma=0;
        
        for j=1:n
            
            if j~=i
                sigma=sigma+A(i,j)*x(j);
            end
            
        end
        
        x(i)=(1/A(i,i))*(B(i)-sigma);
    end
    
    itr=itr+1;
    normVal=abs(xold-x);
end
toc
disp('Solution of the system is :');
disp(x)
disp('Converged in Number of iterations: ')
disp(itr)
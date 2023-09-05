clear; clc; close all;


syms x i j

f(x,i,j) = ( (i+1)*x^i - i*x^(i-1)  ) * ( (j+1)*x^j - j*x^(j-1)  ) ;

B(i,j) = int( f(x,i,j) , x , [0  1] )   ;


for q = 1:3
    
    B = zeros(q,q);
    F = zeros(1,q)';
    for i=1:q
        for j=1:q
        
            B(i,j) = ( (i+1)*(j+1)/(i+j+1) + i*j/(i+j-1) - ( j*(i+1)+(j+1)*i )/(i+j)  ) ;
        end
            F(i) = -1/( (i+2)*(i+3) ) ;
    end
    C = inv(B)*F ; 


end
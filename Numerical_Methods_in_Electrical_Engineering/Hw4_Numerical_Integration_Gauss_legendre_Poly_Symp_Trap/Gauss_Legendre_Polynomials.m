clear; clc;
a= 0 ; 
b=pi;
N = 10;

syms x
f(x)=((cos(pi/2*cos(x))).^2)./sin(x);

% g(x)=(b-a)/2*f(1/2*((b-a)*x+(b-a)));

% Legendre Polynomials  --> Usable in [-1,1]
% P_m = zeros(1,N);
P_m(1) = sym(1) ;
P_m(2) = x;
for m=2:N-1

P_m(m+1) = (2*(m-1)+1)/((m-1)+1)*x*P_m(m) - ((m-1))/((m-1)+1)*P_m(m-1)  ; 

end

X_roots = cell(1,N);
for m = 1:N
roots = vpasolve(legendreP(m,x) == 0);
X_roots{1,m} = double(roots);
end

%% Find H for m = N;
H = zeros(1,N);
p_roots  = X_roots{1,N};
Leg_polys(x) = P_m(1,N);
Sum = 0;

for m=1:N
    
    H(1,m) = 2*(1-p_roots(m)^2) / ( ( N )^2 *  ( double(Leg_polys(p_roots(m))) )^2  );

    if (   (b-a)/2*p_roots(m) + (b+a)/2 == 0 )
    Sum  = Sum +0; % Function Value at 0 is 0!    
    else
    Sum = Sum + H(1,m)*double(  f( (b-a)/2*p_roots(m) + (b+a)/2   )   ) ;
    end

end

Sum  = (b-a)/2 * Sum;

disp(Sum);


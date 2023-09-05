%% Part B - Q1:
clear; clc; close all;
%% Init:
c  = 3e+08;
fc = c ;
Lambda = fc/c;
K = 2*pi/Lambda;
T = 1;


a = Lambda/5;

z = (0.001:a/100:10*a/10)';
y = 5*a;
x = 5*a;

z_p = 8*a/10;
x_p = 0;
y_p = 0;




%% Modal Solution:

tic;
Max_iter = 3;
g_Modal = zeros(length(z),Max_iter);
X = zeros(length(z),1);
k_rho_n = zeros(1,Max_iter);
H_02 = k_rho_n;

for n=1:Max_iter
    k_rho_n(n) = sqrt( K^2-((n-1)*pi/a)^2 );
    Rho =  sqrt( (x-x_p)^2 + (y-y_p)^2 );
    H_02(n) = besselj(0,k_rho_n(n) * Rho) - 1j*bessely(0,k_rho_n(n) * Rho) ;  
    if(abs(H_02(n))>10^100)
       H_02(n) = 0;     
    end
    g_Modal(:,n) =  g_Modal(:,n) + H_02(n)*sin((n*pi/a)*z_p)*sin((n*pi/a)*z)/(1j*2*a);
    X = X  + H_02(n)*sin((n*pi/a)*z_p)*sin((n*pi/a)*z)/(1j*2*a);
end

Delta_Modal = toc;


%% Real Images Solution:

% Max_iter = 10;
tic;
g_real_images = zeros(length(z),Max_iter);

for n=1:Max_iter
    for l=0:3
        Z_L = [ abs(z-z_p) , z+z_p , 2*a-abs(z-z_p), 2*a-(z+z_p) ];
        R_nl =sqrt( (x-x_p)^2 + (y-y_p)^2 + (Z_L(:,l+1) + 2*(n-1)*a).^2    ) ;
        g_real_images(:,n) = g_real_images(:,n) + ((-1)^l )*   exp(-1j*K*R_nl)./(4*pi*R_nl) ;
    end
end
Delta_Real = toc;


%% Complex Images Solution:

g_Complex_images = zeros(length(z),Max_iter);

T = 10;
M = 3;
Beta = 0;

Solution_Complex = Solution_calc(M,T,a,K,Lambda,z,z_p,Beta);
a_m = Solution_Complex.A;
b_m = Solution_Complex.B_m;

tic;
for m=1:M
    for l=0:3
        Z_L_Complex = [ abs(z-z_p) , z+z_p , abs(-z+z_p) , -z-z_p  ] ;
        R_ml = sqrt( Rho^2 + (Z_L_Complex(l+1)+1j*b_m(m)).^2  );
        g_Complex_images(:,n) = g_Complex_images(:,n) + a_m(m)*((-1)^l)*exp(-1j*K*R_ml)./(4*pi*R_ml)  ;
    end


end
R_0 = sqrt( Rho^2 + (z-z_p).^2  ); 
R_1 = sqrt( Rho^2 + (z+z_p).^2  ); 

g_Complex_images(:,end) = g_Complex_images(:,end) + exp(-1j*K*R_0)./(4*pi*R_0) - exp(-1j*K*R_1)./(4*pi*R_1)   ;
Delta_Complex = toc;


%% Comparison:


figure()
subplot(3,1,1)
plot(z/a , abs(g_real_images(:,end))  )
title("Max_{iter} = "+ Max_iter+" and a = "+a/Lambda+" Lambda")
xlabel("z/a")
legend("Real Images Solution")
ylabel("")
grid on

subplot(3,1,2)
plot(z/a , abs(g_Modal(:,end)))
title("z_{source} = "+z_p/a+" a")
xlabel("z/a")
legend("Modal Solution")
ylabel("")
grid on


subplot(3,1,3)
plot(z/a , abs(g_Complex_images(:,end)))
title("M = "+M)
xlabel("z/a")
legend("Complex images Solution")
ylabel("")
grid on


%% Time Analysis:

figure()
stem(Delta_Complex);
hold on
stem(Delta_Real);
stem(Delta_Modal);
legend("Complex Images", "Real Images" , "Modal Solution")
grid on
title("Time Analysis of each Method")

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%


function F_approx = F_approx_calc(Beta,z_point,z_point_source,A,B_m)
syms Beta_z
F_tilde_approx(Beta_z) = reshape(A,1,[])*exp(B_m*Beta_z);

% a= 0.3429;
% syms Beta_z
% F(Beta_z)  =  exp(-1j*Beta_z*2*a) ./ ( 1- exp(-1j*Beta_z*2*a) )  ; 

F_approx = 1./(2j*Beta) .* ( exp(-1j*Beta*(z_point - z_point_source))-exp(-1j*Beta*(z_point + z_point_source))+double(F_tilde_approx(Beta)).*...
    ( exp(-1j*Beta*(z_point - z_point_source))-exp(-1j*Beta*(z_point + z_point_source))+exp(-1j*Beta*(-z_point + z_point_source))-...
    exp(-1j*Beta*(-z_point - z_point_source))  )     ) ;

end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function F_exact = F_exact_calc(Beta,z_point,z_point_source,a)

syms Beta_z
F(Beta_z)  =  exp(-1j*Beta_z*2*a) ./ ( 1- exp(-1j*Beta_z*2*a) )  ; 

F_exact = 1./(2j*Beta) .* ( exp(-1j*Beta*(z_point - z_point_source))-exp(-1j*Beta*(z_point + z_point_source))+double(F(Beta)).*...
    ( exp(-1j*Beta*(z_point - z_point_source))-exp(-1j*Beta*(z_point + z_point_source))+exp(-1j*Beta*(-z_point + z_point_source))-...
    exp(-1j*Beta*(-z_point - z_point_source))  )     ) ;

end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function Object_out  = Prony_calc(M,T,a,K)

syms Beta_z
F(Beta_z)  =  exp(-1j*Beta_z*2*a) ./ ( 1- exp(-1j*Beta_z*2*a) )  ; 

t  = linspace(0,T,2*M); 
L1 = length(t);

B_Z = K*( (1-t/T) - 1j*t  ) ;

for i=1:L1
    F_B_z  =  double( F(B_Z) ) ;
end

Out_object= Prony_Main(M,F_B_z,B_Z);

A = Out_object.A;
B_m = Out_object.B_m;
Mu = Out_object.Mu; 
Alpha = Out_object.Alpha ;

Object_out.K = K;
Object_out.Alpha = Alpha;
Object_out.Mu = Mu;
Object_out.M = M;
Object_out.T = T;
Object_out.A    =  A;
Object_out.B_m  = B_m;

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function Out_object= Prony_Main(M,F_B_z,X_vec)

Coeff_Alpha = zeros(M,M);

for p=1:M

    Coeff_Alpha(p,:) = flip(F_B_z(1+p-1:M+p-1)) ;

end
Sampled_Points = reshape(F_B_z(M+1:2*M),[],1) ;
Alpha          = -inv(Coeff_Alpha)*  Sampled_Points  ;


Mu = roots([1 ; Alpha]);
B_m  = log(Mu);

Coeff_A = zeros(M,M);

for p=1:M
    Coeff_A(p,:) = (Mu.^(p-1)).*ones(M,1) ;
end

Sampled_Points_2 = reshape(F_B_z(1:M),[],1); 
A = inv(Coeff_A) * Sampled_Points_2;

Out_object.A =A;
Out_object.B_m =B_m;
Out_object.Mu = Mu;
Out_object.Alpha = Alpha;


end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function Solution_1 = Solution_calc(M,T,a,K,Lambda,z_point,z_point_source,Beta)

Solution_1  = Prony_calc(M,T,a,K);
A = Solution_1.A ;
B_m = Solution_1.B_m;

% F_exact = F_exact_calc(Beta,z_point,z_point_source,a);
% F_approx = F_approx_calc(Beta,z_point,z_point_source,A,B_m);



Solution_1.z_point = z_point;
Solution_1.z_point_source = z_point_source;
Solution_1.Beta = Beta;
% Solution_1.F_exact = F_exact;
% Solution_1.F_approx = F_approx;
Solution_1.Lambda = Lambda;
Solution_1.a =a;


end








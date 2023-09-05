
clear; close all; clc;

% M = 3;
% T = 1;
% a  = 0.3;
% K = 10;
% 
% 
% 
% Object_out  = Prony_calc(M,T,a,K);

% M = 2;
% F_B_z  = [32,20,14,11];
% X_vec  = [0,1,2,3]; 
% 
% c = 3e+08;
% fc = c;
% Lambda = c/fc;
% M =5;
% T =10;
% a=Lambda/5;
% K = 2*pi/Lambda;
% 
% 
% syms Beta_z
% F(Beta_z)  =  exp(-1j*Beta_z*2*a) ./ ( 1- exp(-1j*Beta_z*2*a) )  ;
% 
% t  = linspace(0.0001,T,2*M); 
% L1 = length(t);
% 
% B_Z = K*( (1-t/T) - 1j*t  ) ;

% F_B_z  =  double( F(B_Z) ) ;
% X_vec = B_Z;

X_vec = [0.5, 2 , 2.5, 2 , 10 , 15];
F_B_z = exp(-X_vec)./(1-exp(-X_vec));
M = 3;



Out_object= Prony_Main(M,F_B_z,X_vec);
Estimated  = conj(Out_object.A)'*exp(Out_object.B_m.*X_vec);
figure()
plot(abs(Estimated))
hold on
plot(abs(F_B_z))
legend("Estimated","Original")
ylabel("Amplitude (Abs)")
xlabel("Index")
grid on


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
clear; clc; close all;




%% Init:
c = 3e+08;   % Electromagnetics wave's  speed in  Air
f = 10e+06; % 10MHz 



Lambda = c/f; % Wavelength


a = 1e-3*Lambda; 
r = a;    % Radius of Wires
D = 2*a;  % Diameter of each Wire

%% Preprocessing:
N = 8;
L = Lambda/2 ; % Antenna Length (sum of Dipole's length) 

delta_l = L/(N+1);
delta_l_index = 100; % Each Triangle Delta_l is equal to 100 indexes in l_vec



l_vec = linspace(0,L,(N+1)*delta_l_index);

[W,W2] = W_calc(N,delta_l_index,l_vec);

%% Processing:
L1 = Lambda;
N1= 80;
Antenna_L1_N1 = Total_Worker(N1,L1,a,f,c,0);



%% Input Impedance Value:

disp(Antenna_L1_N1.Z_in1)

disp(Antenna_L1_N1.Z_in2)

%% Pattern Draw:


Antenna_L1_N1  = Pattern_draw(Antenna_L1_N1);

%% Plotting:

% V:
% for incident Wave we have:
figure()
stem(abs(Antenna_L1_N1.V2))
grid on
title("Gap Voltage Fed")
xlabel("Nodes Index")
ylabel("Amplitude")
legend("Voltage")


figure()
stem(abs(Antenna_L1_N1.V1))
grid on
title("Gap Voltage Fed E_i based")
xlabel("Nodes Index")
ylabel("Amplitude")
legend("Voltage")

figure()
stem(abs(Antenna_L1_N1.I1))
grid on
title("I obtained from E_i")
xlabel("Nodes Index")
ylabel("Amplitude")
legend("I_1")



figure()
stem(abs(Antenna_L1_N1.I2))
grid on
title("I obtained from V Gap")
xlabel("Nodes Index")
ylabel("Amplitude")
legend("I_2")


figure()
imagesc(abs(Antenna_L1_N1.Z1))
grid on
title("Z_1 obtained from E_i")
xlabel("Nodes Index")
ylabel("Amplitude")



figure()
imagesc(abs(Antenna_L1_N1.Z2))
grid on
title("Z_2 obtained from V Gap")
xlabel("Nodes Index")
ylabel("Amplitude")
colorbar

%% Convergence Test:


L = Lambda/2;
N_vec = [20,50,60,70,80,100,110,130,150];
N = N_vec(1);
A1_N1 = Total_Worker(N,L,a,f,c,0);

N = N_vec(2);
A1_N2 = Total_Worker(N,L,a,f,c,0);

N = N_vec(3);
A1_N3 = Total_Worker(N,L,a,f,c,0);

N = N_vec(4);
A1_N4 = Total_Worker(N,L,a,f,c,0);

N = N_vec(5);
A1_N5 = Total_Worker(N,L,a,f,c,0);

N = N_vec(6);
A1_N6 = Total_Worker(N,L,a,f,c,0);

N = N_vec(7);
A1_N7 = Total_Worker(N,L,a,f,c,0);

N = N_vec(8);
A1_N8 = Total_Worker(N,L,a,f,c,0);

N = N_vec(9);
A1_N9 = Total_Worker(N,L,a,f,c,0);



A1_vec = {  A1_N1 , A1_N2 , A1_N3 , A1_N4 ,A1_N5 ,A1_N6 , A1_N7 ,A1_N8 ,A1_N9  };
Z_A1 = zeros(length(N_vec),1);

for i=1:length(N_vec)
    Temp = A1_vec{i};
    Z_A1(i) =  Temp.Z_in2;

end

figure()
stem( N_vec, abs(Z_A1))
title("Convergence Analysis for Input Impedance")
xlabel("Number of Nodes")
grid on





%% MATLAB Calculation using Antenna Tool BOX:
Antenna_Length = [Lambda/4 , Lambda/2 , 3*Lambda/4 , Lambda ];

Z_in = zeros(length(Antenna_Length) , 1);

for i=1:length(Antenna_Length)
    Z_in(i) = Dipole_Antenna_exact_Z(Antenna_Length(i),a,Lambda);
    disp("For Dipole with length equal to: "+Antenna_Length(i) )
    disp(Z_in(i));

end

%% Functions:

function V_m  = G_m_calc(W , E_i , m ,l_vec , delta_l_index)

    Axis   = zeros(1,length(l_vec));
    Axis(1:2*delta_l_index) = W(1:2*delta_l_index);
   
    T_m    = circshift( Axis , (m-1)*delta_l_index  ) ;

    V_m    = sum( T_m.*E_i , 'all') ; 

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function  M = M_calc( m , n , delta_l , d )

delta = delta_calc(m,n,delta_l);    

alpha = 2*delta_l + delta;
Beta  =  delta_l  + delta;
Gamma =  delta_l  + delta;

if(m==n)% self Term:
         M = 2*delta_l*asinh(delta_l/d) - 2*sqrt(delta_l^2+ d^2) + 2*d;
else
         M = alpha*asinh(alpha/d) - Beta*asinh(Beta/d)  - Gamma*asinh(Gamma/d) + delta*asinh(delta/d) ...
             - sqrt(alpha^2 +d^2) + sqrt(Beta^2 +d^2) + sqrt(Gamma^2 +d^2) - sqrt(delta^2 +d^2) ;
end

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function delta = delta_calc(m,n,delta_l)
    delta = (abs(m-n)-1)*delta_l;
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function  Z_in = Dipole_Antenna_exact_Z(l,a,Lambda)
    eta = 120*pi;
    k = 2*pi/Lambda;
    C = 0.5772;

    X  = eta/(4*pi) * ( 2*fresnels(k*l) + cos(k*l)*(2*fresnels(k*l) - fresnels(2*k*l) ) ...
                  -sin(k*l)*( 2*fresnelc(k*l)-fresnelc(2*k*l)-fresnelc(2*k*a^2/l)  ) ) ;
    
    Rr  = eta/(2*pi) * (  ...
         C + log(k*l) - fresnelc(k*l) ...
         + 1/2*sin(k*l)* (fresnels(2*k*l)-2*fresnels(k*l)) ...
         + 1/2*cos(k*l)* (C + log(k*l/2)+ fresnelc(2*k*l) - 2*fresnelc(k*l) ) ...
                      )  ; %  C = 0.5772 (Euler s constant)
    R_in = Rr/(1e-3+sin(k*l/2))^2;
    X_in = X/(sin(k*l/2))^2;

    Z_in = R_in + 1j*X_in;



end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [W,W2] = W_calc(N,delta_l_index,l_vec)

W = zeros(1,length(l_vec));

for i=1 : N+1
    if(mod(i,2)==1)
        W(1 + (i-1)*delta_l_index:i*delta_l_index) = 0.5*linspace(0,1,delta_l_index);
    else
        W(1 + (i-1)*delta_l_index:(i)*delta_l_index) =  -0.5*( linspace(1,2,delta_l_index) )+1;
    end
end
W2 =  circshift([W(1:end-delta_l_index),zeros(1,delta_l_index) ],delta_l_index);
if(mod(N,2)==0)
 W(end-delta_l_index:end) = 0;   
else
 W2(end-delta_l_index:end) = 0;
end

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function Object_Antenna= Pattern_draw(Object_Antenna)

L =Object_Antenna.L ;
I = Object_Antenna.I2;
Lambda = Object_Antenna.Lambda;
delta_l =Object_Antenna.delta_l;
k =Object_Antenna.k;

theta = -180: 0.1 :180 ;
zn = linspace(-L/2,L/2,length(I))';
Pattern = sind(theta).*sum( delta_l*I.*exp(1j*k*zn*cosd(theta)) ) ;


Object_Antenna.theta = theta;
Object_Antenna.Pattern_theta = Pattern;

figure()
polarplot(pi*theta/180, abs(Pattern))
% plot( abs(Pattern) , theta )
% hold on
% plot( -abs(Pattern) , theta )
title("Pattern of Antenna for \lambda  = "+Lambda+" and L = "+L/Lambda+"*Lambda")
grid on

end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function Z2 = Z_calc(N,w,mu0,delta_l,eps0,PSAI )


Z2 = zeros(N,N);
Z2(1,1) = 1j*w*mu0*(delta_l^2) *PSAI(1,1) + 1/(1j*w*eps0) * ( 2*PSAI(1,1)-2*PSAI(1,2) ) ; % Self Term

m=1;
for n=2:N
        Z2(m,n) = 1j*w*mu0*(delta_l^2)*PSAI(1,abs(m-n)+1) + 1/(1j*w*eps0) * ( 2*PSAI(1,abs(m-n)+1) - PSAI(1,abs(m-n)+2) - PSAI(1,abs(m-n)) )  ;
end
   

for m=2:N
    for n=1:N
        % Self Terms:
        if(m==n)
            Z2(m,n) = Z2(1,1);
        else
            Z2(m,n) = Z2(1,abs(m-n)+1) ;
        end


    end
end


end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function  PSAI = PSAI_calc(N,delta_l , d,k)
    PSAI =zeros(N+1,N+1);
    M = PSAI;

    for m=1:N+1
        for n=1:N+1
            M(m,n)    = M_calc( m , n , delta_l , d );
            Reff      =  (delta_l*delta_l)/M(m,n);
            PSAI(m,n) = exp(-1j*k*Reff)/(4*pi*Reff);
        end
    end

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function Total_Object = Total_Worker(N,L,a,f,c,draw)



Lambda = c/f ;
delta_l = L/(N+1);
delta_l_index = 100; % Each Triangle Delta_l is equal to 100 indexes in l_vec
l_vec = linspace(0,L,(N+1)*delta_l_index);

Total_Object = struct();
Total_Object.Lambda        = Lambda;
Total_Object.delta_l       = delta_l;
Total_Object.delta_l_index = delta_l_index ; 
Total_Object.l_vec = l_vec;
Total_Object.draw = draw;
Total_Object.N = N;
Total_Object.L = L;
Total_Object.f = f;


[W,W2] = W_calc(N,delta_l_index,l_vec);

if(draw==1)
    figure()
    plot(l_vec,W);
    
    grid on
    hold on
    plot(l_vec , W2)
    for i=1:N+1
        plot( i*delta_l*ones(1,10) , linspace(0,max(W),10),'r--');
    end
    legend("W","W2")
end


Total_Object.W = W;
Total_Object.W2 = W2;



V1 = zeros(N+1,1);
E_i = zeros(1,(N+1)*delta_l_index)  ;

mid_point = floor(length(E_i)/2) ; % Tahrik az vasat



E_i(mid_point) = 1  ;
for m = 1:N+1 
    V1(m)  = G_m_calc(W , E_i , m ,l_vec , delta_l_index);
end

V2 = zeros(N+1,1);
mid_point2 = floor(length(V2)/2);
V2(mid_point2) = 1;

Total_Object.mid_point = mid_point;
Total_Object.E_i = E_i;
Total_Object.V1 = V1;
Total_Object.V2 = V2;


M = zeros(N+1,N+1);
Z1 = M;
PSAI1 = M;
% PSAI_f = PSAI;

d = a;
k = 2*pi/Lambda;  % wave number

w = 2*pi*f; % Rad/m

mu0  = 4*pi*1e-07; % H/m
eps0 = 8.85*1e-12; % F/m

Total_Object.eps0 = eps0;
Total_Object.mu0 = mu0;
Total_Object.w = w;
Total_Object.d = a;
Total_Object.k = k;



for m = 1:N+1
    for n=1:N+1
        M(m,n)    = M_calc( m , n , delta_l , d );
        Reff      =  (delta_l*delta_l)/M(m,n);
        PSAI1(m,n) = exp(-1j*k*Reff)/(4*pi*Reff);
        % PSAI_f(m,n) =   ;
        
        if( (m==N+1) || (n==N+1)  )
            Z1(m,n) = 1j*w*mu0*delta_l*delta_l*PSAI1(m,n) + ...
                (1/(1j*w*eps0)*( 0+PSAI1(m,n)- 0 - 0 ) ) ;
        else
            Z1(m,n)    = 1j*w*mu0*delta_l*delta_l*PSAI1(m,n) +...
                (1/(1j*w*eps0)*( PSAI1(m+1,n+1)+PSAI1(m,n)- PSAI1(m+1,n) - PSAI1(m,n+1) ) ) ;
        end

    end
end

PSAI2 = PSAI_calc(N,delta_l , d,k);
Z2 = Z_calc(N,w,mu0,delta_l,eps0,PSAI2 );

% First Check whether the Z2 is illconditioned:
dZ2 = decomposition(Z2);
is_ILL_Cond = isIllConditioned(dZ2);

if(is_ILL_Cond)
    disp("Z2 is ill Conditioned!!!")
end

I2 = inv(Z2)*V2(1:end-1);


Total_Object.PSAI2 = PSAI2;
Total_Object.Z2 = Z2;



Total_Object.M =M;
Total_Object.PSAI1 = PSAI1;
Total_Object.Z1 = Z1;


I1 = inv(Z1)*V1;

Z_in1 = V1(floor(mid_point/delta_l_index))/I1(floor(mid_point/delta_l_index));

Z_in2 = V2((mid_point2))/I2((mid_point2));

% disp("Impedance for Antenna with L = "+L/Lambda+"*Lambda: ")
% disp((Z_in))

Total_Object.Z_in1 = Z_in1;
Total_Object.Z_in2 = Z_in2;
Total_Object.I1 = I1;
Total_Object.I2 = I2;

end


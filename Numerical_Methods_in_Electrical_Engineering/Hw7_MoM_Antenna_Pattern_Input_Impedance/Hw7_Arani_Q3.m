%% HW7_ARani_Q3:

clear all; close all; clc;

%% Init:

c = 3e+08;   % Electromagnetics wave's  speed in  Air
f = 10e+06; % 10MHz 



Lambda = c/f; % Wavelength


a = 1e-3*Lambda; 
r = a;    % Radius of Wires
D = 2*a;  % Diameter of each Wire


%% Preprocessing:
N = 40;
L = Lambda/2 ; % Antenna Length (sum of Dipole's length) 

delta_l = L/(N+1);
delta_l_index = 100; % Each Triangle Delta_l is equal to 100 indexes in l_vec



l_vec = linspace(0,L,(N+1)*delta_l_index);

% [W,W2] = W_calc(N,delta_l_index,l_vec);


%% Processing:
L1 = Lambda;
N1= 80;
Feed_portion = 2;

L_tot            = [0.457 ; 0.46 ; 0.44 ]*Lambda;
Feed_portion_tot = [2 ; 2 ; 2 ];
% Antenna_O1_N1 = Total_Worker_Multi_object(N1,L1,a,f,c,0,Feed_portion);
Total_Object = Total_Worker_Multi_object(N1,L_tot ,a ,f,c,Feed_portion_tot );

%%

I_TOT = Total_Object.I_TOT;

I1 = I_TOT(1:N1);
I2 = I_TOT(N1+1:2*N1);
I3 = I_TOT(2*N1+1:3*N1);


Total_Object.Super_I = [I1 ; I2 ; I3] ; 

%% Pattern Draw:

Total_Object = Pattern_draw_Total(Total_Object);




  %% Functions:


function Total_Object =Pattern_draw_Total(Total_Object)


L1 =Total_Object.L1 ;
L2 =Total_Object.L2 ;
L3 =Total_Object.L3 ;

Super_I = Total_Object.Super_I;
N1 = Total_Object.N;

I1 = Super_I(1:N1);
I2 = Super_I(N1+1:2*N1);
I3 = Super_I(2*N1+1:3*N1);

D = Total_Object.D;

Lambda = Total_Object.Lambda;
delta_l =Total_Object.delta_l;
k =Total_Object.k;

theta = -180: 0.1 :180 ;



zn1 = linspace(-L1/2,L1/2,length(I1))';
Pattern1 = sind(theta).*sum( delta_l*I1.*exp(1j*k*zn1*cosd(theta)) ) ;

zn2 = linspace(-L2/2,L2/2,length(I2))';
Pattern2 = sind(theta).*sum( delta_l*I2.*exp(1j*k*zn2*cosd(theta)) ) ;

zn3 = linspace(-L3/2,L3/2,length(I3))';
Pattern3 = sind(theta).*sum( delta_l*I3.*exp(1j*k*zn3*cosd(theta)) ) ;


Pattern_TOT = exp(1j*k*D(1,1)*sind(theta)).*Pattern1 + ...
              exp(1j*k*D(2,1)*sind(theta)).*Pattern2 + ...
              exp(1j*k*D(3,1)*sind(theta)).*Pattern3 ;


Total_Object.theta = theta;
Total_Object.Pattern1_theta = Pattern1;
Total_Object.Pattern2_theta = Pattern2;
Total_Object.Pattern3_theta = Pattern3;
Total_Object.Pattern_TOT = Pattern_TOT  ;

figure()
polarplot(pi*theta/180, abs(Pattern_TOT))
title("Pattern of 3 Antennas")
grid on


  end



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
function Total_Object = Total_Worker(N,L,a,f,c,draw,Feed_portion)



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

mid_point = floor(length(E_i)/Feed_portion) ; % Tahrik az vasat



E_i(mid_point) = 1  ;
for m = 1:N+1 
    V1(m)  = G_m_calc(W , E_i , m ,l_vec , delta_l_index);
end

V2 = zeros(N+1,1);
mid_point2 = floor(length(V2)/Feed_portion);
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
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function Total_Objects = Total_Worker_Multi_object(N,L_tot ,a ,f,c,Feed_portion_tot )

L1 = L_tot(1,1);
L2 = L_tot(2,1);
L3 = L_tot(3,1);

Total_Objects.L1 = L1;
Total_Objects.L2 = L2;
Total_Objects.L3 = L3;

Feed_portion1 = Feed_portion_tot(1,1);
Feed_portion2 = Feed_portion_tot(2,1);
Feed_portion3 = Feed_portion_tot(3,1);

Total_Objects.Feed_portion1 = Feed_portion1;
Total_Objects.Feed_portion2 = Feed_portion2;
Total_Objects.Feed_portion3 = Feed_portion3;



Antenna_element_1 = Total_Worker(N,L1,a,f,c,0,Feed_portion1);
Antenna_element_2 = Total_Worker(N,L2,a,f,c,0,Feed_portion2);
Antenna_element_3 = Total_Worker(N,L3,a,f,c,0,Feed_portion3);


Total_Objects.Antenna_element_1 = Antenna_element_1;
Total_Objects.Antenna_element_2 = Antenna_element_2;
Total_Objects.Antenna_element_3 = Antenna_element_3;


Total_Objects.delta_l  = Antenna_element_3.delta_l;
Total_Objects.f = f;
Total_Objects.k = Antenna_element_3.k;




Z_11  = Antenna_element_1.Z2;
Z_22  = Antenna_element_2.Z2;
Z_33  = Antenna_element_3.Z2;

Total_Objects.Z_11 = Z_11;
Total_Objects.Z_22 = Z_22;
Total_Objects.Z_33 = Z_33;


Lambda = Antenna_element_1.Lambda; 

Total_Objects.Lambda = Lambda;
% Mutual Term between Antenna 1 and 2:
d1_2 = 0.25*Lambda;
Antenna_element_1_2 = Total_Worker(N,L1,a+d1_2,f,c,0,Feed_portion1);
Z_12 = Antenna_element_1_2.Z2;

Total_Objects.d1_2 = d1_2;
Total_Objects.Antenna_element_1_2 = Antenna_element_1_2;
Total_Objects.Z_12 = Z_12;

% Mutual Term between Antenna 1 and 3:
d1_3 = 0.25*Lambda + 0.31*Lambda ;
Antenna_element_1_3 = Total_Worker(N,L1,a+d1_3,f,c,0,Feed_portion1);
Z_13 = Antenna_element_1_3.Z2;


Total_Objects.d1_3 = d1_3;
Total_Objects.Antenna_element_1_3 = Antenna_element_1_3;
Total_Objects.Z_13 = Z_13;

% Mutual Term between Antenna 2 and 3:
d2_3 =  0.31*Lambda;
Antenna_element_2_3 = Total_Worker(N,L2,a+d2_3,f,c,0,Feed_portion2);
Z_23 = Antenna_element_2_3.Z2;


Total_Objects.d2_3 = d2_3;
Total_Objects.Antenna_element_2_3 = Antenna_element_2_3;
Total_Objects.Z_23 = Z_23;


Total_Objects.D = [0 ; d1_2 ; d1_3];

Z_TOT = [ Z_11 , Z_12 , Z_13 ;...
          Z_12 , Z_22 , Z_23 ;...
          Z_13 , Z_23 , Z_33 ] ;

V1 = Antenna_element_1.V2;
V2 = Antenna_element_2.V2;
V3 = Antenna_element_3.V2;

Total_Objects.Z_TOT = Z_TOT;
Total_Objects.V1 = V1;
Total_Objects.V2 = V2;
Total_Objects.V3 = V3;



V_TOT = [ V1(1:end-1) ; V2(1:end-1) ; V3(1:end-1)  ];

I_TOT  =  inv(Z_TOT) * V_TOT     ;


Total_Objects.V_TOT = V_TOT;
Total_Objects.I_TOT = I_TOT;
Total_Objects.N = N;

end






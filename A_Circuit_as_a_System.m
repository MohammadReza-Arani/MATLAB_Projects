clear;
clc;

% Step response of a system:
syms x

P(x) = x^3  +  9/4*x^2   +  7/4*x  +9/8 ;

disp("eigen value and vectors of P :")
[U1,V1]=eig(P)


P2 =[x+3/4 0 0.5; 0 x -3; -0.5 0.5 x];
disp("eigen value and vectors of P2 :")
[U2,V2]=eig(P2)

disp("Natural frequencies of the circuit: ")
roots([1 9/4 7/4 9/8])
%% 
s= tf('s');

sys = s^3 + 9/4*s^2 + 7/4*s + 9/8 ;
try
    step(sys)
catch
    disp("Not a Causal System ==> No STep response in Time")
end

%% 
A = [-3/4 0 -1/2; 0 0 3; 1/2 -1/2 0] ;
B= [0.5 0 0]' ;

C= [0.5 0 0 ; 0 1 0] ;

D= [0 ; 0] ;

%% solve the problem:
M=A-s;

X= M\B;

%% step response based on s(variable) in laplace transform:

step(X(3))
title('Step Response of ix')
legend('Response')
grid on


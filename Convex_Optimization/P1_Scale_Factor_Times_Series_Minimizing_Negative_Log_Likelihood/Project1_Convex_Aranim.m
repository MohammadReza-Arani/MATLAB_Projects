%% Convex Optimization Project 1
% Mohammad Reza Arani
% 810100511
% 13/May/2022
%% Covar_series_data.m
clear; clc; 
rng('default');
rng(2);

n = 10;
F = rand(10, 10);
Sigma = F'*F;
T = 10;
a = [0.2, 0.1, 0.2, 0.4, 0.8, 1.0, 1.0, 0.8, 0.7, 0.8];
x = zeros(n, T);
y = zeros(n, T);
for t = 1:T
    x(:, t) = mvnrnd(zeros(1,n), a(t)*Sigma);
    y(:, t) = mvnrnd(zeros(1,n), a(t)*Sigma);
end
%%

X=zeros(1,T);
X(1)=2; X(2)=-2;
P = toeplitz(X);

P(1,1) =1; P(end,end) =1;

for t =1:1:T-1
    P(t+1,t)=0;
end

% x and Sigma is defined in initialization:
c = zeros(1,T);
for t = 1 : T
    c(t) = 1/2*x(:,t)'*pinv(Sigma)*x(:,t);
end

%% CVX


Lambda = [0.01 , 1 , 100] ;
c = c';
figure()
S = zeros(length(Lambda),T) ;

for i=1:length(Lambda)
cvx_begin
        variables s(T);
        minimize( ones(1,T)*(T/2*s+c.*exp(-s)) + Lambda(i) * s'*P*s );
cvx_end
disp("Solution for Lambda = "+Lambda(i)+ ": ")
disp(s)
plot(1:length(s),exp(s));
hold on
S(i,:) = s'; 
end
grid on
title("a_t VS T samples");
xlabel("Samples");
ylabel("Amplitude");
legend("a_\lambda_1","a_\lambda_2","a_\lambda_3")


%% Validation with Given y
J1 = 0; J2 = J1; J3 = J2; 
for t = 1 : T
    J1 = J1 + T/2*log(exp(S(1,t))) + 0.5*y(:, t)'*pinv(Sigma)*y(:, t)./exp(S(1,t));
    J2 = J2 + T/2*log(exp(S(2,t))) + 0.5*y(:, t)'*pinv(Sigma)*y(:, t)./exp(S(2,t));
    J3 = J3 + T/2*log(exp(S(3,t))) + 0.5*y(:, t)'*pinv(Sigma)*y(:, t)./exp(S(3,t));
end
%%
clc;
formatSpec = 'J1 is %4.2f meters and J2 is %4.2f and J3 is %4.2f \n';
fprintf(formatSpec,J1,J2,J3)

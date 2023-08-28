%% Hw-2 BSS Not in Live Script
clear; clc; close all;
%% Part-1

% Uniform Disttribution:
T = 1000;

a1 =-3*ones(1,T) ; b1 = 3*ones(1,T) ;
s1 = unifrnd(a1,b1);

a2 =-2*ones(1,T) ; b2 = 2*ones(1,T) ;
s2 = unifrnd(a2,b2);

% removing the bias from the distribution:
s1 = s1 - mean(s1);
s2 = s2 - mean(s2);

S = [s1 ; s2]; % 2*T

A = [1 -2; 2 -1; 3 -2]; % 3*2

Y = A*S ; % 3*T

figure()
scatter(S(1,:),S(2,:))
grid on
title("Scatter Plot of Data Distribution")
legend("Generated Data")

figure(1)
scatter3(Y(1,:) , Y(2,: ) , Y(3,:))% scatter3(X,Y,Z)
grid on
title("Scatter Plot of Data Distribution")
legend("Received Data")

%% EIG decomposition:
Rx = Y*Y' ; % YY^T
[U,Lambda] =  eig(Rx) ;

% [V,D] = eig(A) returns diagonal matrix D of eigenvalues and matrix V whose columns are the corresponding right eigenvectors, so that A*V = V*D.
disp("Eigen-Vectors")
disp(U) % Eigen-Vectors
disp("Eigen-Values")
disp(Lambda) % Eigen-Values

disp("Trace of D: "+ trace(Lambda)+" and Trace of Rx: "+ trace(Rx)) % Trace of eigne-Values = Trace of Rx

%%

u3 = U(:,1);  Lambda_3 = Lambda(1,1);
u2 = U(:,2);  Lambda_2 = Lambda(2,2);
u1 = U(:,3);  Lambda_1 = Lambda(3,3);

disp("u1:")
disp(u1); % The best Vector to map Data to --> Based on Lambda values
disp("u2:")
disp(u2); % Second Best Vector
disp("u3:")
disp(u3); % The worst!

disp("Lambda_3")
disp(Lambda_3) % we can see that the value of this parameter is 0 --> No Variation in this direction (u3)

% u3TX=0:
disp(sum(abs(u3'*Y)))
% u3TA=0:
disp((abs(u3'*A)))


C = pinv([u1 , u2])*A; % Minimized Data Distance from X using 2 describing Vectors u1,u2
disp("C:")
disp(C)

%% Part-3:
% Before we had:   Y = A*S % 3*T
% Now we will map these data to a new space -->
Z = [u1 , u2]'*Y ; % New Space --> Image SPace

figure(2)
scatter(Z(1,:) , Z(2,:));
grid on
xlabel("u1")
ylabel("u2")
title("Image Space of Data(S)")
%% Part-4:
% Data Properties in New Space:
Rz = Z*Z';
[U_z , Lambda_z]  = eig(Rz); 
% We can see the eigen vectors of the New Space:
disp("New eigen Vectors:")
disp(U_z)

% It is clear that we have 2 dimension and the U shall be rank-2 based on
% data dimension and 0 nullity --> It is satisfied! --> meaning we can span
% the data space using these 2 orthogonal vectors! Also means that the
% range of U matrix covers all the dimensions(2 dimensions)!
disp("New eigen Values:")
disp(Lambda_z)

% All eigen values are non zero! --> As expected

% whiten the Z data:
W = chol(inv(Rz));

Z_w = W*Z ;
Rz_w = Z_w*Z_w';

disp("Rz_whitenned")
disp(Rz_w) % It is obvious that we are dealing with I matrix! 


figure(3)
scatter(Z_w(1,:),Z_w(2,:))
grid on
title("Whitenned Data Scatter Plot")
xlabel("u1");
ylabel("u2")

%% [Q,G,V] = svd(Y)
clc;
[Q,G,V] = svd(Y);  % Y = Q*G*V'

% Previously we had U from eig
disp("U was:")
disp(U)

disp("SQRT of Lambda: ") % The Same as S
disp(abs(sqrt(Lambda)))

disp(G(:,1:3))
% We can see G(:,1:3) = abs(sqrt(Lambda)) -- >> eigen values of Ry
% Other columns of G are all 0!

disp(G(find(G>1e-10))) % Non-zero elements of diagnal Matrix G:


% Number of nonzero columns of G is 2 as can be observed!
rank_Y = rank(Y); % Non-zero Columns of diagnal Matrix G_d --> actually here, G_d is equal to G
% G is a diagnal matrix of singular values in descending form!
disp("Rank of Y is: "+rank_Y)

% 2 First Rows of V' show Z
Check = Z*V;
disp(sum(sum(abs(Check)>1e-10)));
%% Part5:
% With respect to above relationship between S,F and Z:
F_z = S*pinv(Z);
F_z_w =  S*pinv(Z_w);
Temp =(S*V)';
disp(Temp(1:10,:)) % S is orthogonal to all rows of V but the first two of them!

%% Part-6:

% Total Energy of Received Signal:
Etot = trace(Lambda);

Proportional_E = [Lambda_1 , Lambda_2, Lambda_3 ] / Etot;
disp(Proportional_E); % with respect to Energy Values -- >> to contain at least 90% of energies 

% in received Signal
% It is wise to Map Data over u1 solely! == >>> This preserves 96% of the
% Energy , More than desired 90%
Y_PCA = U(:,1)'*Y; % Mapping Over u1
figure(5)
scatter(Y_PCA(1,:), zeros(length(Y_PCA)))
grid on
title("Data over u1")
xlabel("Time Samples")
ylabel("Data")


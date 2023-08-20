%% Gauss-Seidel:  Solve Ax=B
clear;
clc;

A=[44 -1 1 2 12; 4 68 1 3 -2;   -2 1 35 12 0;13 1 -1 33 5; 1 9 6 -2 53];
B=[7;-21;15;8;22];

N=5;
n = length(B);
X = zeros(n,1);
e = ones(n,1);
%%Check if the matrix A is diagonally dominant
tic
for i = 1:n
    j = 1:n;
    j(i) = [];
    C = abs(A(i,j));
    Check(i,1) = abs(A(i,i)) - sum(B); % Is the diagonal value greater than the remaining row values combined?
    if Check(i) < 0
        fprintf('The matrix is not strictly diagonally dominant at row %2i\n\n',i)
    end
end

C_n=10^-2;
iteration = 0;
while max(e) > C_n%check error
    iteration = iteration + 1;
    Z = X;  % save current values to calculate error later
    for i = 1:N
        j = 1:N; % define an array of the coefficients' elements
        j(i) = [];  % eliminate the unknow's coefficient from the remaining coefficients
        Xtemp = X;  % copy the unknows to a new variable
        Xtemp(i) = [];  % eliminate the unknown under question from the set of values
        X(i,1) = (B(i,1) - sum(A(i,j) * Xtemp)) / A(i,i);
    end
    Xs = X;
    e = sqrt((X - Z).^2);
end
%%Display Results
Xs
iteration
toc

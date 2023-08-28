%% Hw-1 BSS Not in Live Script
clear; clc; close all;
%% Part-2:

[X,Y] = meshgrid(-10:.1:10);
Z = X.^2 + Y.^2 -4*X - 6*Y+13 + X.*Y;
figure()
mesh(X,Y,10*log10(Z))
grid on
title("db of f(x)")
colorbar


%% Show Contours


figure()
contour(X,Y,Z,'ShowText','on')
grid on
title("Contours of db(f(x))")
colorbar


%% Hessian

H = [2,1;1,2];
[V,D] = eig(H);
disp(D);

%% Solve g=0

A  = [2,1;1,2];
B = [4;6];

Glob_Point = A\B;

%% Show the Point

figure()
contour(X,Y,Z,'ShowText','on')
grid on
title("Contours of db(f(x))")
colorbar
hold on
plot(Glob_Point(1),Glob_Point(2),'r*')


figure()
mesh(X,Y,10*log10(Z))
grid on
title("db of f(x)")
colorbar
hold on
X1 = Glob_Point(1);
Y1 = Glob_Point(2);
Z1 =  X1.^2 + Y1.^2 -4*X1 - 6*Y1+13 + X1.*Y1;
plot3(X1,Y1,10*log10(Z1),'r*')


%% Steepest Descent Method
Points_old = [6,6];
Points_new = Points_old;
syms a b
g(a,b) = [ 2*a-4+b , 2*b-6+a]  ; 
fun(a,b) = a.^2+b^2-4*a-6*b+13+a*b;

figure()
contour(X,Y,Z,'ShowText','on')
grid on
title("Contours of db(f(x))")
colorbar
hold on

cost = zeros(1,100);
miu = 1e-2;
for iter=1:100
    x = Points_new(1) ;
    y = Points_new(2) ;

    cost(iter) = double(fun(x,y)) - Z1;

    Points_new =  Points_new  - miu*(g(x,y)) ;% *double(fun(x,y))
    Points_new = double(Points_new);
    z = double(fun(x,y));
    plot3(x,y,10*log10(z),'r*')


end

hold off


figure()
plot(1:iter,cost)
title("Cost function VS Iters")
grid on
xlabel("Iter")
%% Newton Method:



Points_old = [6;6];
Points_new = Points_old;
syms a b
g(a,b) = [ 2*a-4+b , 2*b-6+a]  ; 
fun(a,b) = a.^2+b^2-4*a-6*b+13+a*b;

figure()
contour(X,Y,Z,'ShowText','on')
grid on
title("Contours of db(f(x))")
colorbar
hold on

cost_Newton = zeros(1,100);
H = [2,1;1,2];
miu = 1e0;
for iter=1:100
    x = Points_new(1) ;
    y = Points_new(2) ;

    cost_Newton(iter) = double(fun(x,y)) - Z1;

    Points_new =  Points_new  - miu*inv(H)*(g(x,y)') ;% *double(fun(x,y))
    Points_new = double(Points_new);

    z = double(fun(x,y));
    plot3(x,y,10*log10(z),'r*')


end

hold off


figure()
plot(1:iter,cost_Newton)
title("Cost_{Newton} function VS Iters")
grid on
xlabel("Iter")
%% Alternation Minimization:

Points_old = [6;6];
Points_new = Points_old;

figure()
contour(X,Y,Z,'ShowText','on')
grid on
title("Contours of db(f(x))")
colorbar
hold on

cost_Alter = zeros(1,100);
for i=1:100


    x = Points_new(1);
    y = Points_new(2);

    z = double(fun(x,y));
    plot3(x,y,10*log10(z),'r*')
    cost_Alter(i) = double(fun(x,y)) - Z1;

    Points_new(1) = (4-Points_new(2))/2;
    Points_new(2) = (6-Points_new(1))/2;

end


figure()
plot(1:iter,cost_Alter)
title("Cost_{Alter} function VS Iters")
grid on
xlabel("Iter")

%% Steepest Descent Method
Points_old = [6,6];
Points_new = Points_old;
syms a b
g(a,b) = [ 2*a-4+b , 2*b-6+a]  ; 
fun(a,b) = a.^2+b^2-4*a-6*b+13+a*b;

figure()
contour(X,Y,Z,'ShowText','on')
grid on
title("Contours of db(f(x))")
colorbar
hold on

cost_constrained = zeros(1,100);
miu = 1e-1;

for iter=1:100
    x = Points_new(1) ;
    y = Points_new(2) ;

    cost_constrained(iter) = double(fun(x,y)) ;

    Points_new =  Points_new  - miu*(g(x,y)) ;% *double(fun(x,y))
    Points_new = double(Points_new)/norm(double(Points_new)); % implementation of Constraint

    z = double(fun(x,y));
    plot3(x,y,10*log10(z),'r*')


end

hold off


figure()
plot(1:iter,cost_constrained)
title("Cost_{constrained} function VS Iters")
grid on
xlabel("Iter") 


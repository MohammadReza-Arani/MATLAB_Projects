%% Hw-1 BSS Not in Live Script
clear; clc; close all;

%% Part-1:
% X = rand(n) returns an n-by-n matrix of uniformly distributed random numbers.
T = 1000;
S1 = 6*rand(1,T)-3;
S1 = S1';
S2 = 4*rand(1,T)-2;
S2  = S2';


S = [S1,S2];
A = [0.6,0.8;0.8,-0.6];

X = A*S';
x1 = X(1,:);
x2 = X(2,:);

%% Scatter Plot:
figure()
scatter(x1,x2);
grid on
title(" x2 vs x1 Scatter plot ")
xlabel("  x1 ");
ylabel("  x2 ");

figure()
hist(x1)
title("Histogram of x1")
grid on


figure()
hist(x2)
title("Histogram of x2")
grid on
%% 

Ts=1;
fs=1/Ts;
fsample=100;

Tsample=1/fsample;
t=-5: Tsample :5; % time vector

TAU1 = 6 ;
TAU2 = 4 ;

h1=heaviside(t+TAU1/2)-heaviside(t-TAU1/2);
h2=heaviside(t+TAU2/2)-heaviside(t-TAU2/2);

New_Dist=conv(0.8*h1,-0.6*h2);
figure()
plot(linspace(-5,5,length(New_Dist)),abs(New_Dist));
grid on
title("x2 expecting Dist.")

%%
% P = zeros(1,1);
% 
% for index_x1=1:length(x1)
%     targets_x1 = find( x1 == x1(index_x1) );
%     [~,targets_x2 ]= min( x2(targets_x1) );
% 
%     P(end+1,1) = targets_x2 ;
%     
% end
% 
% figure()
% scatter(x1(P),x2(P))
% grid on


Point_1 = [max(x1),x2(find(x1==max(x1)))];

Point_2 = [(x1( find(x2==max(x2)) )),max(x2)] ;

Point_3 = [min(x1),x2(find(x1==min(x1)))];

Point_4 = [(x1( find(x2==min(x2)) )),min(x2)] ;

P_points =  [Point_1;Point_2;Point_3;Point_4 ] ;
X_points =   P_points(:,1);
Y_points =   P_points(:,2);

figure()
scatter( X_points,Y_points ,"filled");
grid on
hold on
scatter( x1,x2 );
legend("Certain Points","Data Dist")

Diff_a1_2 = Point_1-Point_3;
Sum_a1_2 = Point_2-Point_4;

a_1 = 0.5*( Sum_a1_2 +  Diff_a1_2);
a_2 = 0.5*( Sum_a1_2 -  Diff_a1_2);

a_1_orginal =  Point_1 -  Point_4;
a_2_orginal =  Point_3 -  Point_4;

A_hat = [a_1' , a_2'];
Compare = A_hat./A;



%% Part-3:
%  Noisy Data
Sigma = 0.5;
Noise = Sigma.^2*randn(T,2);

X_noisy = A*S'+Noise';

x1_noisy = X_noisy(1,:);
x2_noisy = X_noisy(2,:);

%% Scatter Plot:
figure()
scatter(x1_noisy,x2_noisy);
grid on
title(" x2_{Noisy} vs x1_{Noisy} Scatter plot ")
xlabel("  x1_{Noisy} ");
ylabel("  x2_{Noisy} ");


Point_1_Noisy = [max(x1_noisy),x2_noisy(find(x1_noisy==max(x1_noisy)))];

Point_2_Noisy = [(x1_noisy( find(x2_noisy==max(x2_noisy)) )),max(x2_noisy)] ;

Point_3_Noisy = [min(x1_noisy),x2_noisy(find(x1_noisy==min(x1_noisy)))];

Point_4_Noisy = [(x1_noisy( find(x2_noisy==min(x2_noisy)) )),min(x2_noisy)] ;

P_points_Noisy =  [Point_1_Noisy;Point_2_Noisy;Point_3_Noisy;Point_4_Noisy ] ;
X_points_Noisy =   P_points_Noisy(:,1);
Y_points_Noisy =   P_points_Noisy(:,2);

figure()
scatter( X_points_Noisy,Y_points_Noisy ,"filled");
grid on
hold on
scatter(x1_noisy,x2_noisy);
legend("Certain Points","Data Dist")


a_1_orginal_Noisy =  Point_1 -  Point_4;
a_2_orginal_Noisy =  Point_3 -  Point_4;
%% Fit the Curve:

x = linspace(-5,5,T);

m_1 =  a_1_orginal_Noisy(2)/a_1_orginal_Noisy(1);
Curve_1 = m_1*(x-Point_4(1)) + Point_4(2);

m_2 =  a_2_orginal_Noisy(2)/a_2_orginal_Noisy(1);
Curve_2 = m_2*(x-Point_4(1)) + Point_4(2);

figure()
scatter(Curve_2,x,'filled');
scatter(Curve_1,x,'filled');
scatter(x1_noisy,x2_noisy);
hold on
grid on
title("Noisy Data & estimated vectors")


%% Part 4 and 5

figure()
hist(x1,20)
title("Histogram of x1")
grid on


figure()
hist(x2,20)
title("Histogram of x2")
grid on

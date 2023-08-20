clear; clc; close all;

%% Init:
N=100;
Max_iter=10;
Phi = 0.5;
h = 1e-2;


A = zeros(N,N) ;
Y = zeros(Max_iter,N);
%% Boundary Condition
B = zeros(N,1);
B(1,1) = -1;
cntr =0;


for i=1:N
        for j=1:N

            if i==N
                A(i,:) = [zeros(1,N-2) , 2 , -2-(h*Phi)^2 ] ;
            elseif i==1 
                A(i,:) = [ -2-(h*Phi)^2,1,zeros(1,N-2)  ]  ;
            else 
                A(i,:) = [ zeros(1,i-2) , 1, -2-(h*Phi)^2  ,1    , zeros(1,N-(i-2)-3)    ] ;
            end
        
        end
end

Y(:,:,1) = 1;
Y(end,:,:)  = inv( A )*B;



%% Plotting the result:

X = linspace(0,1,N);

figure()

plot(X,Y(end,:))
title("Implicit FDM for N="+N)
grid on
xlabel("")
ylabel("Amplitude")
legend("Solution with Implicit FDM")


%% Second Method

Y_dir = zeros(Max_iter,N);

while(1)

    cntr = cntr+1;
    for i=1:N
        
            if i==N
                Y_dir(N) = 2*Y_dir(N-1)/(2+(h*Phi)^2 ) ;
            elseif i==1
                Y_dir(1) = (1+Y_dir(2))/(2+(h*Phi)^2 ) ;
            else
                Y_dir(i) = (Y_dir(i-1)+Y_dir(i+1))/(2+(h*Phi)^2 ) ;
            end
        
    end
    
%     Y(cntr,:) = inv( A )*B;

    if(cntr==Max_iter)
        break;
    end

end

%% Plotting the result: Y_dir

X = linspace(0,1,N);

figure()

plot(X,Y(end,:))
title("Direct FDM Method for N="+N)
grid on
xlabel("")
ylabel("Amplitude")
legend("Solution with Direct FDM")



%% Analytical Solution:
c1= 0.2689 ; c2 = 0.7311 ;
x_anal = linspace(0,1,N) ;
Y_anal = c1*exp(x_anal/2)+c2*exp(-x_anal/2) ;

figure()
plot(x_anal,Y_anal)
grid on
legend("Analytical Solution")
ylabel("Amplitude")
xlabel("X")
title("Analytical Solution of given Problem")

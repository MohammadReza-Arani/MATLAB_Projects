%% HW-4 BSS-Not In LiveScript: Arani-810100511
clear; clc; close all;
%% Generate S <Source Signal>:
fs = 20; %Hz
ts = 1/fs;

K = 5;

T_rec = K-ts;
t= 0:ts:T_rec;

c = [0.2,0.4, 0.6,-0.1,-0.3];
d = [0.1 0.3 -0.2  0.5 -0.3];

C = repmat(c,fs,1);
D = repmat(d,fs,1);

s1 = (C(:))'.*sin(2*pi*t);
s2 = (D(:))'.*sin(4*pi*t);

S = [s1;s2];

%% Generate X:
A = [0.8,-0.6; 0.6 , 0.8];
X = A*S;

%% Plot the Source (S) & Observation Signal(X):

figure()
subplot(2,1,1)
plot(t,S)
grid on
legend("Source Signal-1","Source Signal-2")
title("Source Signal")
ylabel("Amplitude")

subplot(2,1,2)
plot(t,X)
grid on
legend("Observ Signal-1","Observ Signal-2")
title("Observ Signal")
xlabel("Time [s]")



%% Data Whitening & BSS SOlve -Part B


RX = X*X';
[U_X , Lambda_X] = eig(RX);
% Method 1:
% W_whit  = chol(inv(RX));
% Z = W_whit'*X;
% Method 2:
Whitener =(Lambda_X^(-1/2))*U_X.';
% Z = Whitener*X;
Z=X;
RZ = Z*Z';

[U_Z , Lambda_Z] = eig(RZ);

k=1;
RZ_k1 = Z(:, (k-1)*fs+1:(k)*fs  )*Z(:,(k-1)*fs+1:(k)*fs  )';
k=2;
RZ_k2 = Z( :,(k-1)*fs+1:(k)*fs  )*Z(:,(k-1)*fs+1:(k)*fs  )';

[Q,Lambda_Z_K] = eig(inv(RZ_k2)*RZ_k1);
B = Q';

S_Amp = S/norm(S,"fro");
[~ , B_Hat_Chosen ] = Perm_AMP_Disamb(B,Z,S_Amp);
S_hat = B_Hat_Chosen*Z;

S_hat_perm = S_hat;
% permutation disambiguation: 
% S_hat_perm(1,:) = S_hat(2,:);
% S_hat_perm(2,:) = S_hat(1,:);

% Amplitude disambiguation: 
% S_hat_perm = [1/max(S_hat_perm(1,:)) , 0; 0 ,  1/max(S_hat_perm(2,:))]*S_hat_perm;
% S_hat_perm_AMp = [max(S(1,:)) , 0; 0 ,  max(S(2,:))]*S_hat_perm ;
S_hat_perm_AMp = S_hat_perm/norm(S_hat_perm,'fro'); 
Error = (norm(S_hat_perm_AMp-S_Amp,'fro'))^2/(norm(S_Amp,'fro'))^2 ; 

disp("Error = "+Error);

figure()
subplot(2,1,1)
plot(t,S_hat_perm_AMp(1,:),'r-*')
hold on
plot(t,S_Amp(1,:))
title("Comparison of S^{hat} & S")
grid on
ylabel("Amplitude")
legend("S^{hat} Source 1","S_{True} Source 1" )

subplot(2,1,2)
plot(t,S_hat_perm_AMp(2,:),'r-*')
hold on
plot(t,S_Amp(2,:))
title("Comparison of S^{hat} & S")
grid on
legend("S^{hat} Source 2","S_{True} Source 2" )
xlabel("Time [s]")


%%  JOINT Diagonalization:


B_Hat = orth(rand(size(B)));
% R_1 = zeros(size(RZ));
% R_2 = R_1; 
% Whitener:
% Whitener = 1;
Whitener =(Lambda_X^(-1/2))*U_X.';
Z = Whitener*X;

for iter = 1:10
     R_2 = 0;
     R_1 =0;
    for k=1:K
        RZ_k  = Z(:, (k-1)*fs+1:(k)*fs  )*Z(:,(k-1)*fs+1:(k)*fs  )';

        for j=1:length(B_Hat)
            if (j==1)
                b_j = B_Hat(j,:)' ;
                R_2 = R_2 + (RZ_k*b_j)*(RZ_k*b_j)';
    
            else
                b_j = B_Hat(j,:)' ; % Doubt Here!!!! Choose Columns or Rows of B???
                R_1 = R_1 + (RZ_k*b_j)*(RZ_k*b_j)';
    
            end
        end
    
    end

[V_R1,Lambda_R1] = eig(R_1);
[~, ind] = sort(diag(Lambda_R1),'descend');
V_R1 = V_R1(:, ind);
b_1 =  V_R1(:,end);


% For b_2:
[V_R2,Lambda_R2] = eig(R_2);
[~, ind] = sort(diag(Lambda_R2),'descend');
V_R2 = V_R2(:, ind);
b_2 =  V_R2(:,end); % According to lowest eigen Value

% b_2 = ( eye(size(b_2))  - b_1*b_1' )*b_2 ;
b_2 = b_2 - (b_1'*b_2)*b_2 ;
b_2 = b_2/norm(b_2,2);


B_Hat  = [b_1' ; b_2'];
% B_Hat  = B_Hat; 
[~,B_Hat]  = Perm_AMP_Disamb(B_Hat,Whitener*Z,S_Amp);

end


[Erroes_B ,  B_Hat_CHosen]  = Perm_AMP_Disamb(B_Hat,Whitener*Z,S_Amp);
S_hat_2 = B_Hat_CHosen*(Whitener*Z);
S_hat_perm_AMp_2 = S_hat_2/norm(S_hat_2,'fro'); 
Error = (norm(S_hat_perm_AMp_2-S_Amp,'fro'))^2/(norm(S_Amp,'fro'))^2 ; 
disp("Error = "+Error);

%%
figure()
subplot(2,1,1)
plot(t,S_hat_perm_AMp_2(1,:),'r-*')
hold on
plot(t,S_Amp(1,:))
title("Comparison of S^{hat} & S For k = "+K)
grid on
ylabel("Amplitude")
legend("S^{hat} Source 1","S_{True} Source 1" )

subplot(2,1,2)
plot(t,S_hat_perm_AMp_2(2,:),'r-*')
hold on
plot(t,S_Amp(2,:))
title("Comparison of S^{hat} & S")
grid on
legend("S^{hat} Source 2","S_{True} Source 2" )
xlabel("Time [s]")


%%   PART-D:

W = randn(size(X));
W_normalised = W/norm(W,"fro") ; 

disp("Normlaised W  with Freb. Norm = "+norm(W_normalised,"fro"));

SNR = 20;
SNR_lin = 10^(SNR/10);
sigma = sqrt( norm(X,"fro")^2/SNR_lin );

Y  = X + sigma*W_normalised ;

Window_Num = K;
Num_iter = 1000;
% B_Hat_K_Y = Diag_Joint_Func(Y,fs,Window_Num,Num_iter);
B_Hat_K_Y = Diag_Joint_Func(Y,fs,Window_Num,Num_iter,S_Amp);
[~,B_Hat_Chosen]  = Perm_AMP_Disamb(B_Hat_K_Y,Z,S_Amp);

% S_hat_Y = B_Hat_K_Y*Z;
% S_hat_perm_Y = S_hat_Y;
% 
% S_hat_perm_AMp_Y = -S_hat_perm_Y/norm(S_hat_perm_Y,'fro'); 
S_hat_Y = B_Hat_Chosen*Z;
S_hat_perm_AMp_Y = S_hat_Y/norm(S_hat_Y,'fro'); 
S_Amp = S/norm(S,"fro");
Error_Y = (norm(abs(S_hat_perm_AMp_Y)-abs(S_Amp),'fro'))^2/(norm(S_Amp,'fro'))^2 ; 

disp("Error = "+Error_Y);

%%

figure()
subplot(2,1,1)
plot(t,S_hat_perm_AMp_Y(1,:),'b-*')
hold on
plot(t,S_Amp(1,:))
title("Comparison of S^{hat} & S For k = "+Window_Num)
grid on
ylabel("Amplitude")
legend("S^{hat} Source 1","S_{True} Source 1" )

subplot(2,1,2)
plot(t,S_hat_perm_AMp_Y(2,:),'b-*')
hold on
plot(t,S_Amp(2,:))
title("Comparison of S^{hat} & S")
grid on
legend("S^{hat} Source 2","S_{True} Source 2" )
xlabel("Time [s]")


%% Generate Noise for 100 times then use the Expectation:
Z =X;
Num_of_Trails = 100;
Error_Y = zeros(length(2:K),Num_of_Trails);

S_Amp = S/norm(S,"fro");
SNR = 20;
SNR_lin = 10^(SNR/10);
sigma = sqrt( norm(X,"fro")^2/SNR_lin );
Num_iter = 100;

for k=2:K
    for i=1:Num_of_Trails
    
        W = randn(size(X));
        W_normalised = W/norm(W,"fro") ;  
        %disp("Normlaised W  with Freb. Norm = "+norm(W_normalised,"fro"));
        Y  = X + sigma*W_normalised ;
        Window_Num = k; 
        B_Hat_K_Y = Diag_Joint_Func(Y,fs,Window_Num,Num_iter,S_Amp);
        [~,B_Hat_Chosen]  = Perm_AMP_Disamb(B_Hat_K_Y,Z,S_Amp);

        S_hat_Y = B_Hat_Chosen*Z;
%         S_hat_perm_Y = S_hat_Y;
        S_hat_perm_AMp_Y = S_hat_Y/norm(S_hat_Y,'fro');      
        Error_Y(k-1,i) = (norm((S_hat_perm_AMp_Y)-(S_Amp),'fro'))^2 ; 
        %disp(Error_Y(k-1,i))
        if(isnan(Error_Y(k-1,i)))
            Error_Y(k-1,i) = Error_Y(k-1,i-1);
        end
    
    end
end

Error_Y_Mean = mean(Error_Y,2);
disp(Error_Y_Mean)
%% Plot The E versus k:
figure()
plot(2:K,Error_Y_Mean,'r-*')
grid on
xlabel("k")
ylabel("Mean of Error")
title("Error Value Versus Number of Windows")



%%




















%% Perm_AMP_Disamb
function [Erroes_B,B_Hat_Chosen]  = Perm_AMP_Disamb(B_Hat,Z,S_Amp)
Perm_1 = B_Hat;
Perm_2 = -B_Hat;
Perm_3 = B_Hat([end, 1:end-1],:);
Perm_4 = -B_Hat([end, 1:end-1],:);
Perms = { Perm_1,Perm_2,Perm_3,Perm_4  };
% Estimate S:
S_hat_1 = Perm_1*Z;
S_hat_1_Normalised = S_hat_1/norm(S_hat_1,'fro');
S_hat_2 = Perm_2*Z;
S_hat_2_Normalised = S_hat_2/norm(S_hat_2,'fro');
S_hat_3 = Perm_3*Z;
S_hat_3_Normalised = S_hat_3/norm(S_hat_3,'fro');
S_hat_4 = Perm_4*Z;
S_hat_4_Normalised = S_hat_4/norm(S_hat_4,'fro');
% Calc Error:
Error_1 = (norm((S_hat_1_Normalised)-(S_Amp),'fro'))^2 ; 
Error_2 = (norm((S_hat_2_Normalised)-(S_Amp),'fro'))^2 ; 
Error_3 = (norm((S_hat_3_Normalised)-(S_Amp),'fro'))^2 ; 
Error_4 = (norm((S_hat_4_Normalised)-(S_Amp),'fro'))^2 ; 

Erroes_B = [Error_1,Error_2,Error_3,Error_4];
[ ~, index] = min(Erroes_B);
B_Hat_Chosen = Perms{index};
end



%% Joint Diag Function:



function B_Hat = Diag_Joint_Func(Z,fs,Window_Num,Num_iter,S_Amp)

RZ = Z*Z';
B_Hat = orth(rand(size(RZ)));



    for iter = 1:Num_iter
        R_1 = zeros(size(RZ));
        R_2 = R_1; 
        for k=1:Window_Num
            RZ_k  = Z(:, (k-1)*fs+1:(k)*fs  )*Z(:,(k-1)*fs+1:(k)*fs  )';
    
            for j=1:length(B_Hat)
                if (j==1)
                    b_j = B_Hat(j,:)' ;
                    R_2 = R_2 + (RZ_k*b_j)*(RZ_k*b_j)';
        
                else
                    b_j = B_Hat(j,:)' ; % Doubt Here!!!! Choose Columns or Rows of B???
                    R_1 = R_1 + (RZ_k*b_j)*(RZ_k*b_j)';
        
                end
            end
        
        end
    
    [V_R1,Lambda_R1] = eig(R_1);
    [~, ind] = sort(diag(Lambda_R1),'descend');
    V_R1 = V_R1(:, ind);
    b_1 =  V_R1(:,end);
    
    % For b_2:
    [V_R2,Lambda_R2] = eig(R_2);
    [~, ind] = sort(diag(Lambda_R2),'descend');
    V_R2 = V_R2(:, ind);
    b_2 =  V_R2(:,end); % According to lowest eigen Value
    % b_2 = ( eye(size(b_2))  - b_1*b_1' )*b_2 ;
    b_2 = b_2 - (b_1'*b_2)*b_2 ;
    b_2 = b_2/norm(b_2,2);
    
    B_Hat = [b_1' ; b_2'];
    [~,B_Hat]  = Perm_AMP_Disamb(B_Hat,Z,S_Amp);
    end



end





















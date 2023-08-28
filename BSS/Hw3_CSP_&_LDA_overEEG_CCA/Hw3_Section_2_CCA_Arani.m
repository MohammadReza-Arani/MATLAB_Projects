%% Hw-3-BSS:
%% Section-2: 
clear; clc; close all;
load("hw3-2.mat");
% And Also fs is given:
fs = 250; %Hz

Num_Rows = length(data(:,1,1));
Num_Cols = length(data(1,:,1));
Num_Pages = length(data(1,1,:));

Trec = Num_Cols/fs;
ts = 1/fs;
t = 0:ts:Trec-ts;


%% Generate X for given frequencies:

X_F = cell(length(freq),1);

for i=1:length(X_F)
    % The frequency itself:
    Temp = [sin(2*pi*freq(i)*t); cos(2*pi*freq(i)*t)];
    % Consider Harmonics till 40Hz:
    for j=2:10
        if(j*freq(i)<40)
            Temp = [Temp ; sin(2*pi*j*freq(i)*t); cos(2*pi*j*freq(i)*t)];
        end
    end
    X_F{i,1} = Temp;
end
%% GEVD for each freq:

Ray_Ratio_for_Each_Freq =zeros(length(X_F) , Num_Pages);
for i=1:length(X_F)
    XF = X_F{i,1};
    RX = XF*XF';

    for j=1:Num_Pages
    Y  = data(:,:,j) ;    
    RY = Y*Y';
    RXY = XF*Y';
    RYX = Y*XF';
    COV_Matrix = (RX^(-1/2))*RXY*(RY^(-1))*RYX*(RX^(-1/2));
    I = eye(size(COV_Matrix));
    [U,Lambda] = eig(COV_Matrix,I);
    % Now Like previous Part let's Sort this descendingly:
    [Lambda_2, ind] = sort(diag(Lambda),'descend');
    U = U(:, ind);
    % Choosing the maximum eigne Value and its corresponding eigen vector:
    U_Chosen = U(:,1);
    Ray_Ratio_for_Each_Freq(i,j) = Lambda_2(1);


    end

end

 [Val , pos]= max(Ray_Ratio_for_Each_Freq,[],1);
 Freq_Labels = freq(pos);

 %% COmparison:

 figure();
 plot(Freq_Labels,"r*","LineWidth",2)
 hold on
 plot(label,"b--*")
 grid on
 title("Comparison of True and predicted frequencies");
 legend("Predicted","True",'Location','north');
 xlabel("Frequency index")
 ylabel("Frequency")




%% HW-3----------Section-1:
clear; clc; close all;

load("hw3-1.mat");
% Removing Data Bias:
Num_of_Rows  = length(TrainData_class1(:,1,1));
Num_of_Cols  = length(TrainData_class1(1,:,1));
Num_of_Pages = length(TrainData_class1(1,1,:));

TrainData_class1 = TrainData_class1 - repmat( mean(TrainData_class1, 2) ,1,Num_of_Cols) ;
TrainData_class2 = TrainData_class2 - repmat( mean(TrainData_class2, 2) ,1,Num_of_Cols) ;
TestData         =  TestData        - repmat( mean(TestData, 2) ,1,Num_of_Cols) ;

%% Covariance Calc:
R_class_1 = cell(Num_of_Pages,1);
R_class_2 = R_class_1;
R_class_1_Sum = zeros(Num_of_Rows);
R_class_2_Sum = R_class_1_Sum;


for i=1:Num_of_Pages
R_class_1{i,1} = TrainData_class1(:,:,i)*TrainData_class1(:,:,i)';
R_class_2{i,1} = TrainData_class2(:,:,i)*TrainData_class2(:,:,i)';

R_class_1_Sum = R_class_1{i,1} + R_class_1_Sum ;
R_class_2_Sum = R_class_2{i,1} + R_class_2_Sum ; 

end

% Calculate the Mean of every Page:
R_class_1_Mean  = R_class_1_Sum/Num_of_Pages;
R_class_2_Mean  = R_class_2_Sum/Num_of_Pages;

% Use GEVD (Generalized Eigen-Value Decomposition) over R_class_1_Mean on
% top and R_class_2_Mean in denominator:

[W,Lambda] =  eig(R_class_1_Mean,R_class_2_Mean) ;
[Lambda_2, ind] = sort(diag(Lambda),'descend');
W = W(:, ind);

% Maximizing the Rayleigh Ratio requires the biggest Lambda and its
% corresponding Eigen Vector to be chosen!

% Normalizing W Columns:
Norm_Matrix_of_Cols_W = W'*W;
for j=1:Num_of_Rows
    W(:,j) = W(:,j)/(Norm_Matrix_of_Cols_W(j,j));
end




%% Part A):

 

Z_Class1 = tensorprod(W',TrainData_class1,2,1);
Z_Class2 = tensorprod(W',TrainData_class2,2,1);


%% Plot Filtered Data:
Page_Num = 49;
Filter_Num = [1,30];


Sample_Data_1_Filtered_Class1 = Z_Class1(Filter_Num(1),:,Page_Num);
Sample_Data_2_Filtered_Class1 = Z_Class1(Filter_Num(2),:,Page_Num);

Sample_Data_1_Filtered_Class2 = Z_Class2(Filter_Num(1),:,Page_Num);
Sample_Data_2_Filtered_Class2 = Z_Class2(Filter_Num(2),:,Page_Num);


figure(1)
subplot(4,1,1)
stem(Sample_Data_1_Filtered_Class1)
title("Class 1: Trial Number:"+Page_Num+" Using Filters: "+Filter_Num(1)+" and "+Filter_Num(2))
legend("Filter Number:"+Filter_Num(1))
grid on
subplot(4,1,2)
stem(Sample_Data_2_Filtered_Class1)
xlabel("Data sample")
ylabel("Amplitude")
legend("Filter Number:"+Filter_Num(2))
grid on



subplot(4,1,3)
stem(Sample_Data_1_Filtered_Class2)
title("Class 2: Trial Number:"+Page_Num+" Using Filters: "+Filter_Num(1)+" and "+Filter_Num(2))
legend("Filter Number:"+Filter_Num(1))
grid on
subplot(4,1,4)
stem(Sample_Data_2_Filtered_Class2)
xlabel("Data sample")
ylabel("Amplitude")
legend("Filter Number:"+Filter_Num(2))
grid on
%% Calculate Variance for Given Data:

Var_Sample_Data_1_Filtered_Class1 = var(Sample_Data_1_Filtered_Class1);
Var_Sample_Data_2_Filtered_Class1 = var(Sample_Data_2_Filtered_Class1);

Var_Sample_Data_1_Filtered_Class2 = var(Sample_Data_1_Filtered_Class2);
Var_Sample_Data_2_Filtered_Class2 = var(Sample_Data_2_Filtered_Class2);

figure()
plot(Var_Sample_Data_1_Filtered_Class1,Var_Sample_Data_2_Filtered_Class1,'r*')
hold on
plot(Var_Sample_Data_1_Filtered_Class2,Var_Sample_Data_2_Filtered_Class2,'g*')
legend("Class 1 Data","Class 2 Data",'Location','north')
grid on
title("Scatterplot of Data Points for each class")
xlabel("Feature 1")
ylabel("Feature 2")


%% Part B:


figure()
subplot(2,1,1)
plot(W(:,1))
title("CSP W filter First and Last Comparison")
grid on
subplot(2,1,2)
plot(W(:,Num_of_Rows))
grid on
xlabel("Sensor Index")
ylabel("Importance")


%% Part C:

% Calculate LDA (Linear Dsicriminator Analysis) Classifier Weights:

% Our Data After Choosing 14 Filters from 30:
Chosen_Num_Rows = 7;
Z_Class1_Chosen = zeros(Chosen_Num_Rows*2,Num_of_Cols,Num_of_Pages);
Z_Class2_Chosen = Z_Class1_Chosen;

Z_Class1_Chosen(1:7,:,:) = Z_Class1(1:Chosen_Num_Rows,:,:) ;
Z_Class1_Chosen(8:14,:,:) = Z_Class1(Num_of_Rows-Chosen_Num_Rows+1:Num_of_Rows,:,:);

Z_Class2_Chosen(1:7,:,:) = Z_Class2(1:Chosen_Num_Rows,:,:) ;
Z_Class2_Chosen(8:14,:,:) = Z_Class2(Num_of_Rows-Chosen_Num_Rows+1:Num_of_Rows,:,:);

% Calculate the Mean:

MIU_Class_1 = mean(Z_Class1_Chosen, 2);
MIU_Class_2 = mean(Z_Class2_Chosen, 2);

% Calculate the Co-Variance:

Z_Class1_Chosen_unbiased = Z_Class1_Chosen   - repmat( MIU_Class_1 ,1,Num_of_Cols);
Z_Class2_Chosen_unbiased = Z_Class2_Chosen   - repmat( MIU_Class_2 ,1,Num_of_Cols);

% COV_Matrix_Class_1 =  tensorprod( Z_Class1_Chosen_unbiased , Z_Class1_Chosen_unbiased, 2,2 ) ;

COV_Matrix_Class_1=cell(Num_of_Pages,1);
COV_Matrix_Class_2 =COV_Matrix_Class_1;

COV_Matrix_Sum_Class_1 = 0;
COV_Matrix_Sum_Class_2 = COV_Matrix_Sum_Class_1;
for i=1:Num_of_Pages
     COV_Matrix_Class_1{i,1} = (Z_Class1_Chosen_unbiased(:,:,i)* Z_Class1_Chosen_unbiased(:,:,i)')/Num_of_Cols;
     COV_Matrix_Class_2{i,1} = (Z_Class2_Chosen_unbiased(:,:,i)* Z_Class2_Chosen_unbiased(:,:,i)')/Num_of_Cols;
    
     COV_Matrix_Sum_Class_1 = COV_Matrix_Sum_Class_1 +  COV_Matrix_Class_1{i,1};
     COV_Matrix_Sum_Class_2 = COV_Matrix_Sum_Class_2 +  COV_Matrix_Class_2{i,1};
end

COV_Matrix_Mean_Class_1 = COV_Matrix_Sum_Class_1/Num_of_Pages;
COV_Matrix_Mean_Class_2 = COV_Matrix_Sum_Class_2/Num_of_Pages;

MIU_Class_1_Mean_Trial = mean(MIU_Class_1,3);
MIU_Class_2_Mean_Trial = mean(MIU_Class_2,3);

[W_LDA , LAMBDA_LDA] = eig( (MIU_Class_1_Mean_Trial-MIU_Class_2_Mean_Trial)*(MIU_Class_1_Mean_Trial-MIU_Class_2_Mean_Trial)' ...
                          , ( COV_Matrix_Mean_Class_1+COV_Matrix_Mean_Class_2 )  );

[LAMBDA_LDA_2, ind] = sort(diag(LAMBDA_LDA),'descend');
W_LDA = W_LDA(:, ind);

%% Normalization of W_LDA:

% Normalizing W Columns:
Norm_Matrix_of_Cols_W_LDA = W_LDA'*W_LDA;
for j=1:Chosen_Num_Rows*2
    W_LDA(:,j) = W_LDA(:,j)/(Norm_Matrix_of_Cols_W(j,j));
end

% Choosing the Maximum Discriminator!
W_LDA_Chosen = W_LDA(:,1);

New_Mean_Class_1 = W_LDA_Chosen'*MIU_Class_1_Mean_Trial;
New_Mean_Class_2 = W_LDA_Chosen'*MIU_Class_2_Mean_Trial;

C = mean([New_Mean_Class_1 , New_Mean_Class_2]) ; % The Threshold



%% Part D:


Cntr1 = 0;
Cntr2 = Cntr1;

New_Data_Projected_LDA_Class_1 = tensorprod( W_LDA_Chosen, Z_Class1_Chosen , 1,1 );
New_Data_Projected_LDA_Class_2 = tensorprod( W_LDA_Chosen, Z_Class2_Chosen , 1,1 );

for i=1:Num_of_Pages
    Data_Mean_Check_Class_1 = mean(New_Data_Projected_LDA_Class_1(1,:,i)) ;
    Data_Mean_Check_Class_2 = mean(New_Data_Projected_LDA_Class_2(1,:,i)) ;
    if( 100*C < Data_Mean_Check_Class_1   ) % Belongs to Class 1
        Cntr1 = Cntr1+1;
    end
    if( 0.0001*C > Data_Mean_Check_Class_2   ) % Belongs to Class 1
        Cntr2 = Cntr2+1;
    end
end

Accuracy_over_Train_Data_Class_1 = Cntr1/Num_of_Pages;
Accuracy_over_Train_Data_Class_2 = Cntr2/Num_of_Pages;

%% Now Test Data:

Z_Test = tensorprod(W',TestData,2,1);
Num_of_Pages_Test = length( TestData(1,1,:));

Z_Test_Chosen = zeros(Chosen_Num_Rows*2,Num_of_Cols,Num_of_Pages_Test);
Z_Test_Chosen(1:7,:,:) = Z_Test(1:Chosen_Num_Rows,:,:) ;
Z_Test_Chosen(8:14,:,:) = Z_Test(Num_of_Rows-Chosen_Num_Rows+1:Num_of_Rows,:,:);


New_Data_Projected_LDA_Test = tensorprod( W_LDA_Chosen, Z_Test_Chosen , 1,1 );

My_Labels = zeros(1,Num_of_Pages_Test);
for i=1:Num_of_Pages_Test
    Data_Mean_Check_Test = mean(New_Data_Projected_LDA_Test(1,:,i)) ;
    
    if( C < Data_Mean_Check_Test   ) % Belongs to Class 1
       My_Labels(i) = 1;
    else
       My_Labels(i) = 2;
    end

end
Error = (TestLabel - My_Labels);
Accuracy_over_Test_Data = sum(Error==0)/length(Error);

%%

figure();
plot(My_Labels,'*');
hold on;
plot(TestLabel,'+');
legend('predicted','True');
xlabel('index');
ylabel('label');
grid on


clear; clc;

load("hw6-part3.mat");

S_hat = OMP_Mine_CAlc( X(:,1)  , D , N0);
[s_OMP] = OMP(X(:,1) ,D , N0);

[S_hat , s_OMP]


function [s_OMP] = OMP(x ,D , N0)
    xr = x;
    Dr = D;
    idx_OMP = zeros(1, N0);
    for i = 1 : N0
        [val, idx_OMP(i)] = max(abs(xr'*Dr),[],'omitnan');
        xr = xr - x'*Dr(:, idx_OMP(i))*Dr(:,idx_OMP(i));
        if i > 1
            xr = xr - D(:,idx_OMP(1:i))*(pinv(D(:,idx_OMP(1:i)))*xr);
        end
        Dr(:,idx_OMP(i)) = nan;
    end

    s_OMP = zeros(size(D, 2), 1);
    s_OMP(idx_OMP) = pinv(D(:, idx_OMP))*x;
end




function  S_hat = OMP_Mine_CAlc( X  , D , N0)


xr = X;
[~ , Col_D] = size(D);
Chosen_IDX = inf+zeros(1,Col_D);
S_hat = zeros(Col_D,1); 
for  i=1:N0
    
    B = D;

    if(i>1)
        B(:,Chosen_IDX(1:i-1)) = []; % Removing previously chosen Indices!
    end

    Corr_matrix = abs(repmat(xr,1,Col_D-i+1).*B);
    Corr_sum = sum(Corr_matrix,1); % Summation for each Column -->> a Row Vector
    [Value,Idx] = max(Corr_sum); % Choosing Best Fitted

    if(sum(Idx>Chosen_IDX)>0)
        Idx  = Idx + sum(Idx>Chosen_IDX) ; % Update the Index with respect to the Original D 
    end

    Chosen_IDX(i) =  Idx;
    xr = xr - X'*D(:,Idx)*D(:,Idx); % xr = x  - <x,di>di

    % Update Coefficients:
    if(i>1)
        D_sub_omp = D(:,Chosen_IDX(1:i)) ;
        xr = xr - D(:,Chosen_IDX(1:i))*(pinv(D_sub_omp)*xr);
    end
    % Stop Criteria:
%     if(norm(x - D(:,Chosen_IDX_3(1:i)) ,2)<thresh )
%         break;
%     end
    S_hat(Chosen_IDX(i)) = pinv(D(:,i))*X;
end

% delta_t_OMP = toc;
% disp("Elapsed Time (OMP): "+ delta_t_OMP+"(s)");
end

%% LU Decomposition:
clear; clc;
    n = input("n = : ");
    
    A = input("A = : ") ;
    
    if(size(A)~=[n,n])
        disp("A must be n*n")
        return
    else
    disp("Matlab LU Decompose Output: ")    
    [L1,U1] = lu(A); % Matlab Solution
    disp([L1,U1]);
    disp("Our LU Decompose Output: ")    
    [L2,U2,P2] = LU_Decomp(A); % Our Function
    disp([L2,U2]);
    
    
    
    
    end

%%
function [L,U,P] = LU_Decomp(A)
NA = size(A,1);
AP = [A eye(NA)]; 
for k = 1:NA - 1
    % Partial Pivoting at AP(k,k) -> Mehvar Giri
    [akx, kx] = max(abs(AP(k:NA,k)));
    if akx < eps
        error('Singular matrix and No LU decomposition')
    end
    mx = k+kx-1;
    if kx > 1 % Row change if necessary
        tmp_row = AP(k,:);
        AP(k,:) = AP(mx,:);
        AP(mx,:) = tmp_row;
    end
    % LU decomposition
    for m = k + 1: NA
        AP(m,k) = AP(m,k)/AP(k,k); % Eq.(2.4.8.2)
        AP(m,k+1:NA) = AP(m,k + 1:NA)-AP(m,k)*AP(k,k + 1:NA); % Eq.(2.4.9)
    end
end
P = AP(1:NA, NA + 1:NA + NA); % Permutation matrix
for m = 1:NA
    for n = 1:NA
        if m == n, L(m,m) = 1.; U(m,m) = AP(m,m);
        elseif m > n, L(m,n) = AP(m,n); U(m,n) = 0.;
        else L(m,n) = 0.; U(m,n) = AP(m,n);
        end
    end
end

end
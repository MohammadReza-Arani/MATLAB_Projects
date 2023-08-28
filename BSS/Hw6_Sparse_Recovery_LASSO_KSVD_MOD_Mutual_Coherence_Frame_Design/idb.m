function [F, cohv] = idb(m, n, desired_coh, lambda, K, gamma0, rho, search_it, tol_stop)

% Minimize frame coherence using a distance barrier function


if nargin < 9
  tol_stop = 0.0001;
end

%F = gen_rand_untf(m,n);
F = normc(randn(m,n));

cohv = zeros(1,K);
coh_min = 1;

for k = 1 : K
  perm = randperm(n-1)+1;
  for j = perm  % a round of the atoms in random order
    d = F(:,j);
    v = F'*d;
    v(j) = 0;
    % what atoms are too close?
    i_minus = find(v>desired_coh);
    i_plus = find(v<-desired_coh);
    % weighted least squares
    w = max(abs(v)/desired_coh, 1);
    % gradient
    %g = F*v;
    g = F*(w.*v);
    %g = zeros(m,1);
    if ~isempty(i_minus)
      g = g + lambda*sum(F(:,i_minus) - repmat(d,1,length(i_minus)), 2);
    end
    if ~isempty(i_plus)
      g = g - lambda*sum(F(:,i_plus) + repmat(d,1,length(i_plus)), 2);
    end
    
    % loop trying to find better atom
    g_now = gamma0;
    c_max = max(abs(v));
    for i = 1 : search_it   % a number of iterations to avoid getting stuck
      dn = d - g_now*g;
      dn = dn / norm(dn);
      vn = F'*dn;
      vn(j)=0;
      if max(abs(vn)) < c_max
        i = i-1;
        break;
      end
      g_now = g_now/2;
    end
    F(:,j) = dn;
  end
  gamma0 = gamma0*rho;
  
  cohv(k) = max(max(abs(F'*F - eye(n))));
  if cohv(k) < coh_min
    coh_min = cohv(k);
    F_best = F;
  end
  if coh_min - desired_coh < tol_stop  % if coherence is near enough the target, stop
    cohv = cohv(1:k);
    break
  end
  %coh_min
end

F = F_best;

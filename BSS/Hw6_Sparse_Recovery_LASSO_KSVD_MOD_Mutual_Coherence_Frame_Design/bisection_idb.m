function [Fbest, coh_best, iter_total_count] = bisection_idb(mu_min, mu_max, nr_it_bis, m, n, ...
                                             K, lambda, gamma0, rho, search_it, tol_stop)

% Compute low coherence frames by bisection and the IDB algorithm


iter_total_count = 0;
coh_best = 1;
for i_bis = 1:nr_it_bis
  %i_bis;
  mu = (mu_min + mu_max)/2;
  [F, coh] = idb(m, n, mu, lambda, K, gamma0, rho, search_it, tol_stop);
  iter_total_count = iter_total_count + length(coh);
  coh_crt = min(coh);
  if coh_crt < coh_best
    coh_best = coh_crt;
    Fbest = F;
  end
  if coh_crt < mu + 10*tol_stop   % successful design
    if coh_crt > mu_max  % if current coherence is larger than current best, stop
      break              % further improvement is impossible, since the interval stays the same
    end
    mu_max = coh_best;
  else
    mu_min = mu;
    mu_max = min(coh_crt, mu_max); % maybe the upper bound can be also improved
  end
  
  %coh_best
  %fprintf("Interval [%f,%f]\n", mu_min, mu_max);

  if mu_max - mu_min < tol_stop  % stop if search interval is very small
    break
  end
end

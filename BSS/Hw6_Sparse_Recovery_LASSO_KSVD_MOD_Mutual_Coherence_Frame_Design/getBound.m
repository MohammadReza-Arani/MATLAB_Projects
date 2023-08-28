function bound = getBound(N,M)
    if M <= N^2
            bound = sqrt((M-N)/(N*(M-1)));
    elseif N^2 < M && M <= 2*(N^2-1)
            bound = max(max(sqrt(1/N), sqrt((2*M-N^2-N)/((N+1)*(M-N)))), 1 - 2*M^(-1/(N-1)));
    elseif M > 2*(N^2-1)
            bound = max(sqrt((2*M-N^2-N)/((N+1)*(M-N))), 1 - 2*M^(-1/(N-1)));
    end
end


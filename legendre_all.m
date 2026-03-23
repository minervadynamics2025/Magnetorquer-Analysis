function P = legendre_all(N, x)
%% Helper: compute Legendre polynomials P_0..P_N at x
    % Returns vector P(1..N+1) where P(k+1) = P_k(x)
    P = zeros(N+1,1);
    P(1) = 1;            % P0
    if N == 0
        return;
    end
    P(2) = x;            % P1
    for n = 2:N
        % recurrence: (n+1) P_{n+1} = (2n+1) x P_n - n P_{n-1}
        P(n+1) = ((2*n-1) * x * P(n) - (n-1) * P(n-1)) / n;
    end
end


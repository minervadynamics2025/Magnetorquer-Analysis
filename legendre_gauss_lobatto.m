function [x, w] = legendre_gauss_lobatto(N)
% LEGENDRE_GAUSS_LOBATTO  Compute LGL nodes and weights
%   [x,w] = legendre_gauss_lobatto(N) returns N+1 LGL nodes x in [-1,1]
%   and corresponding quadrature weights w for polynomial order N (N>=1).
%
%   The interior nodes are the zeros of P_N'(x) (derivative of Legendre P_N),
%   found by Newton-Raphson. We use stable recurrences to evaluate Legendre
%   polynomials P_k(x) for k=0..N and compute derivatives analytically.
%
% Example:
%   [x,w] = legendre_gauss_lobatto(6);
%   disp(x'); disp(w');

    if N < 1
        error('N must be >= 1');
    end

    % Allocate storage
    x = zeros(N+1,1);
    w = zeros(N+1,1);

    % Endpoints
    x(1)   = -1;
    x(end) =  1;

    % If N == 1, trivial 2-point LGL
    if N == 1
        w(1) = 1;
        w(2) = 1;
        return;
    end

    % Initial guesses for interior nodes: Chebyshev-Gauss-Lobatto points
    k = (1:N-1)';
    x0 = -cos(pi * k / N);  % interior initial guesses

    % Newton-Raphson to find zeros of P_N'(x)
    tol = 1e-14;
    maxiter = 100;
    interior = zeros(N-1,1);

    for idx = 1:length(x0)
        xi = x0(idx);
        for iter = 1:maxiter
            % Evaluate Legendre polynomials P_0..P_N at xi
            P = legendre_all(N, xi);   % P(1)=P0, P(N+1)=PN in MATLAB indexing
            PN   = P(N+1);
            PNm1 = P(N);    % P_{N-1}

            % Compute P_N'(xi) using standard identity
            denom = 1 - xi^2;
            if abs(denom) < 1e-16
                break;  % avoid division by zero; xi near endpoints
            end
            PNp = N/denom * (PNm1 - xi * PN);  % P_N'(xi)

            % Compute derivative of P_N' (i.e., P_N'')
            % Need P_{N-2} and derivative P_{N-1}
            if N >= 2
                PNm2 = P(N-1); % P_{N-2}
                % P_{N-1}' using identity
                PN1p = (N-1)/denom * (PNm2 - xi * PNm1);
            else
                PN1p = 0;
            end

            % Analytical derivative of PNp = N*(PNm1 - xi*PN)/denom
            % d/dx PNp = N * ( (PN1p - PN - xi*PNp)*denom - (PNm1 - xi*PN)*(-2*xi) ) / denom^2
            numerator = (PN1p - PN - xi * PNp) * denom + 2*xi * (PNm1 - xi * PN);
            PNpp = N * numerator / (denom^2);

            % Newton update for root of PNp = 0
            dx = - PNp / PNpp;
            xi = xi + dx;

            if abs(dx) < tol
                break;
            end
        end
        interior(idx) = xi;
    end

    % Place interior nodes into x (sorted)
    x(2:N) = interior;

    % Ensure monotonic order (small numerical unsorting can occur)
    x = sort(x);

    % Compute PN(x_i) for weights and then weights
    for i = 1:N+1
        xi = x(i);
        P = legendre_all(N, xi);
        PN = P(N+1);
        w(i) = 2 / ( N * (N+1) * (PN^2) );
    end

    % Ensure symmetry and numerical cleanup
    % (enforce x(1)=-1, x(end)=1 exactly, and symmetric interior)
    x(1) = -1; x(end) = 1;
    % enforce symmetry of nodes and weights to reduce numerical noise
    for i = 1:floor((N+1)/2)
        j = N+2-i;
        xm = 0.5*(x(i) - x(j));
        if abs(xm) > 1e-12
            % nothing
        end
        % enforce symmetry explicitly
        x(i) = -abs(x(i));
        x(j) =  abs(x(i));
        w(j) = w(i);  % weights symmetric
    end
end







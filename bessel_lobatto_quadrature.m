function [r, w] = bessel_lobatto_quadrature(N, nu, R)
% BESSEL_LOBATTO_QUADRATURE  Bessel-Lobatto nodes and quadrature weights
%
%   [r,w] = bessel_lobatto_quadrature(N, nu, R)
%
% Inputs:
%   N  - total number of radial points (including endpoints)  -> size(r) = N
%        (so there are N-1 interior points)
%   nu - order of Bessel function J_nu (typical: nu = 0)
%   R  - radius (domain is [0,R]). Default R = 1.
%
% Outputs:
%   r  - column vector of radial nodes (r(1)=0, r(end)=R)
%   w  - column vector of quadrature weights such that
%        integral_0^R f(r) r dr  ?  sum_k w(k) * f(r(k))
%
% Notes:
%   - Interior nodes are scaled zeros of J_nu'(x).
%   - We compute weights by enforcing discrete orthogonality for the
%     first (N) Bessel modes (using the first N zeros of J_nu).
%
% Example:
%   [r,w] = bessel_lobatto_quadrature(10,0,1);
%   sum(w)   % should be close to R^2/2 * (some consistency), not necessarily 2

    if nargin < 3 || isempty(R), R = 1; end
    if N < 2
        error('N must be >= 2 (two endpoints at least).');
    end

    % number of interior nodes
    Ni = N - 2;

    % --------------------------
    % 1) compute zeros of J_nu'(x)  (need Ni+1 zeros for scaling)
    % --------------------------
    % We'll compute zeros of derivative using fzero with asymptotic guesses.
    Mder = Ni + 1; % we want Ni interior zeros plus one largest for scaling
    beta = zeros(Mder,1);
    % initial asymptotic guesses for derivative zeros
    for k = 1:Mder
        beta_guess = (k + nu/2 - 3/4) * pi;  % good asymptotic initial guess
        % root solve for J_nu'(x) = 0. Use formula J_nu'(x) = (J_{nu-1}-J_{nu+1})/2
        f = @(x) (besselj(nu-1,x) - besselj(nu+1,x))/2;
        % choose bracket around guess to improve fzero
        a = max(1e-8, beta_guess - pi/2);
        b = beta_guess + pi/2;
        try
            beta(k) = fzero(f, [a b]);
        catch
            % fallback single initial guess
            beta(k) = fzero(f, beta_guess);
        end
    end

    % --------------------------
    % 2) form radial nodes r including endpoints 0 and R
    % --------------------------
    r = zeros(N,1);
    r(1) = 0;
    r(end) = R;
    if Ni > 0
        % interior scaled by beta(end)
        for k = 1:Ni
            r(1+k) = R * beta(k) / beta(end);
        end
    end

    % --------------------------
    % 3) compute first N zeros of J_nu (alpha_m) -- need N zeros
    % --------------------------
    Mz = N;  % first N zeros of J_nu
    alpha = zeros(Mz,1);
    for k = 1:Mz
        guess = (k + nu/2 - 1/4) * pi;
        try
            alpha(k) = fzero(@(x) besselj(nu,x), [max(1e-8,guess-pi/2) guess+pi/2]);
        catch
            alpha(k) = fzero(@(x) besselj(nu,x), guess);
        end
    end

    % --------------------------
    % 4) build matrix A and RHS b for linear system A*w = b
    %    A(m,k) = J_nu(alpha_m * r_k / R)^2
    %    b(m)   = (R^2/2) * J_{nu+1}(alpha_m)^2
    % --------------------------
    B = zeros(N, Mz);   % B(k,m) = J_nu(alpha_m * r_k / R)
    for k = 1:N
        arg = alpha * (r(k)/R);   % vector of arguments for m=1..Mz
        B(k, :) = besselj(nu, arg).';
    end

    A = zeros(Mz, N);
    for m = 1:Mz
        for k = 1:N
            A(m,k) = B(k,m)^2;
        end
    end

    bvec = (R^2/2) * (besselj(nu+1, alpha).^2);

    % Solve linear system A * w = bvec  for weights w (size N)
    % A is Mz x N and here Mz==N, so square system
    % Use regularization if matrix is ill-conditioned.
    if cond(A) > 1e12
        warning('A is ill-conditioned (cond=%g). Using regularized solve.', cond(A));
        % Tikhonov regularization small parameter
        lambda = 1e-12;
        w = (A'*A + lambda*eye(N)) \ (A'*bvec);
    else
        w = A \ bvec;
    end

    % convert to column vector just in case
    w = w(:);

    % small numerical cleanup: ensure weights non-negative small tolerance
    w(abs(w) < 1e-16) = 0;

end





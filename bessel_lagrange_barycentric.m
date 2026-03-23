%% bessel_lagrange_barycentric.m
% Represent Bessel function J_nu(r) by Lagrange interpolant using barycentric formula
% Usage: [r_eval, f_interp] = bessel_lagrange_barycentric(nu, r_nodes, r_eval)
function f_interp = bessel_lagrange_barycentric(nu, r_nodes, r_eval)
% Inputs:
%   nu      - order of Bessel function J_nu
%   r_nodes - column vector of interpolation nodes (Nx x 1)
%   r_eval  - row or column vector of evaluation points (1 x M or M x 1)
% Output:
%   f_interp - interpolated values at r_eval (same shape as r_eval)

    % ensure column vector for nodes and column for eval
    r_nodes = r_nodes(:);
    r_eval  = r_eval(:);     % we'll reshape back at end
    N = numel(r_nodes);

    % function values at nodes
    f_nodes = besselj(nu, r_nodes);

    % --- compute barycentric weights w_i = 1 / prod_{j != i} (r_i - r_j) ---
    % Build difference matrix (r_i - r_j)
    Rmat = repmat(r_nodes, 1, N);
    D = Rmat - Rmat.';      % D(i,j) = r_i - r_j
    D(1:(N+1):end) = 1;     % set diagonal to 1 so prod ignores it

    denom = prod(D, 2);     % product over columns for each row (i)
    w = 1 ./ denom;         % barycentric weights (Nx x 1)

    % --- interpolate at r_eval using barycentric formula ---
    M = numel(r_eval);
    f_interp = zeros(M,1);

    % vectorized evaluation: compute matrix (r_eval_k - r_nodes_i)
    % A(k,i) = r_eval(k) - r_nodes(i)
    X = repmat(r_eval, 1, N);
    RN = repmat(r_nodes.', M, 1);
    A = X - RN;   % M x N

    % Detect eval points that coincide with nodes (avoid division by zero)
    tol = 1e-14;
    for k = 1:M
        row = A(k,:);
        idx = find(abs(row) < tol, 1);
        if ~isempty(idx)
            % exact node match: value is f_nodes(idx)
            f_interp(k) = f_nodes(idx);
        else
            % barycentric formula: f(x) = sum_i (w_i/(x - r_i) * f_i) / sum_i (w_i/(x - r_i))
            invdiff = w.' ./ row;        % 1 x N
            numerator = invdiff * f_nodes; % scalar
            denominator = sum(invdiff);
            f_interp(k) = numerator / denominator;
        end
    end

    % return same orientation as input r_eval
    f_interp = reshape(f_interp, size(r_eval));
end



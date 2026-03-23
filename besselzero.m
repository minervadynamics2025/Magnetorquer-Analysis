function z = besselzero(m, N)
% Computes first N zeros of BesselJ_m using fzero

    z = zeros(1, N);
    guess = (1:N) * pi + m*pi/2 - pi/4;  % asymptotic initial guess

    for k = 1:N
        z(k) = fzero(@(x) besselj(m, x), guess(k));
    end
end



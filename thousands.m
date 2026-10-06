function s = thousands(x)
%THOUSANDS Integer with a LaTeX thousands separator, e.g. 2166 -> '2{,}166'.
n = round(x);
if n >= 1000
    s = sprintf('%d{,}%03d', floor(n/1000), mod(n,1000));
else
    s = sprintf('%d', n);
end
end

function s = infOr( x )
%INFOR Integer string for the LaTeX macros, or \infty for Inf.
if isinf(x), s = '\infty'; else, s = sprintf('%.0f', x); end
end

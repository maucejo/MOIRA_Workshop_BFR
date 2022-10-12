function [tau_n, tau_f] = estimate_precision(A, b, W, lambda)
    
N = length(b);
denom = b'*((A*(W\A') + lambda*eye(N))\b);
tau_f = real(N/denom);
tau_n = tau_f/lambda;
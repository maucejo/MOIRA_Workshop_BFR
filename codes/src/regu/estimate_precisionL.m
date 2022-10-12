function [tau_n, tau_f] = estimate_precisionL(A, b, L, lambda)

Ac = A/L;

N = length(b);
denom = b'*((Ac*Ac' + lambda*eye(N))\b);
tau_f = real(N/denom);
tau_n = tau_f/lambda;
function [X, tau_n, tau_f, q] = OBR(A, b, q0, varargin)

% Initialization
an = 1;
bn = 1e-18;
af = 1;
bf = 1e-18;
tol = 1e-8;
[N, M] = size(A);

[X, lambda, W] = reg_add(A, b, q0);
[~, tau_f] = estimate_precision(A, b, W, lambda);
X0 = X;

crit = 1;
while crit > tol
    q = optim_q(X, tau_f);
    tau_f = (M + q*(af - 1))/(q*(bf + norm(X, q)^q));
    tau_n = (N + an - 1)/(bn + norm(b - A*X)^2);
    X = reg_(A, b, q, tau_f/tau_n, X0);
    
    crit = norm(X - X0, 1)/norm(X0, 1);
    X0 = X;
end

end

%%
function qopt = optim_q(X, tau_f)
aq = 1;
bq = 1e-18;
M = length(X);

q = linspace(0.05, 2, 500);
fq = M*log(gamma(1./q)) - M*log(tau_f)./q - (M + aq - 1).*log(q) + bq*q + tau_f*sum(abs(X).^q);

qopt = q(fq == min(fq)); 

end


function [X, tau_n, tau_f, q] = OBR_L(A, b, q0, L)

% Initialization
an = 1;
bn = 1e-18;
af = 1;
bf = 1e-18;
tol = 1e-8;
[M, N] = size(L);

[X, lambda, W] = reg_addL(A, b, q0, L);
[~, tau_f] = estimate_precisionL(A, b, sqrt(W)*L, lambda);
X0 = X;

crit = 1;
while crit > tol
    LX = L*X;
    q = optim_q(LX, tau_f);
    tau_f = (M + q*(af - 1))/(q*(bf + norm(LX, q)^q));
    tau_n = (N + an - 1)/(bn + norm(b - A*X)^2);
    X = regL_(A, b, q, L, tau_f/tau_n, X0);
    
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


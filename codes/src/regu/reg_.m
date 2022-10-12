function [X, WR] = reg_(A, b, q, lambda, X0)
%% Initialization
[~, N] = size(A);
tol = 1e-8;
maxit = 200;
X = X0;

%% Iteration
epsr = 1e-3; 
crit = 1;
ee = 2;
AtA = A'*A;
Atb = A'*b;
while (crit > tol) && (ee <= maxit)
    wr = (q/2)*max(abs(X), epsr).^(q-2);

    WR = sparse(1:N, 1:N, wr(:), N, N);

    X = (AtA + lambda*WR)\Atb;

    crit = norm(X - X0, 1)/norm(X0, 1);
    X0 = X;

    ee = ee + 1;
end
function [X, WR] = regL_(A, b, q, L, lambda, X0)
%% Initialization
N = size(L, 1);
tol = 1e-8;
maxit = 200;
X = X0;

%% Iteration
epsr = 1e-3;
crit = 1;
ee = 2;
while (crit > tol) && (ee <= maxit)
    wr = (q/2)*max(abs(L*X), epsr).^(q-2);
    WR = sparse(1:N, 1:N, wr(:), N, N);
    LambWR = sparse(1:N, 1:N, lambda.*wr, N, N);

    X = (A'*A + (L'*LambWR*L))\(A'*b);
    
    crit = norm(X - X0, 1)/norm(X0, 1);
    X0 = X;
    
    ee = ee + 1;
end
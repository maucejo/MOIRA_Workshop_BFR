function [X, lambda, WR] = reg_addL(A, b, q, L)
%% Initialization
[M, N] = size(L);
tol = 1e-8;
maxit = 200;

%% Initial solution
[U, sm, V] = cgsvd(A, L);

lambda = reg_param(U, sm, b, 1, 0);
X = tikhonov(U, sm, V, b, lambda);

if all(q == 2) && all(all(L'*L == speye(N)))
    WR = eye(N);
    return;
end

X0 = X;

%% Iteration
epsr = 1e-3;
crit = 1;
ee = 2;
while (crit > tol) && (ee <= maxit)
    wr = (q/2)*max(abs(L*X), epsr).^(q-2);
    
    WR = sparse(1:M, 1:M, wr(:), M, M);
    
    [U, sm, V] = cgsvd(A, sqrt(WR)*L);
    lambda = reg_param(U, sm, b, 1, 0);
    X = tikhonov(U, sm, V, b, lambda);
    
    crit = norm(X - X0, 1)/norm(X0, 1);
    X0 = X;
    
    ee = ee + 1;
end
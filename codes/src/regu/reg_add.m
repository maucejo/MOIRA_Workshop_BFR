function [X, lambda, WR] = reg_add(A, b, q)
%% Initialization
[~, N] = size(A);
tol = 1e-8;
maxit = 200;

%% Initial solution
[U, sm, V] = csvd(A);

lambda = reg_param(U, sm, b, 1, 0);
X = tikhonov(U, sm, V, b, lambda);

if q == 2
    WR = eye(N);
    return;
end

X0 = X;

%% Iteration
epsr = 1e-3;
crit = 1;
ee = 2;
while (crit > tol) && (ee <= maxit)
    wr = (q/2)*max(abs(X), epsr).^(q-2);
    
    WR = sparse(1:N, 1:N, wr(:), N, N);
    
    [U, sm, V] = cgsvd(A, sqrt(WR));
    lambda = reg_param(U, sm, b, 1, 0);
    X = tikhonov(U, sm, V, b, lambda);
    
    crit = norm(X - X0, 1)/norm(X0, 1);
    X0 = X;
    
    ee = ee + 1;
end
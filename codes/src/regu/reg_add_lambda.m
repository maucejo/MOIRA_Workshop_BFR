function [X, rho, eta, WR] = reg_add_lambda(A, b, q, lambda)

%% Initialization
M = size(A, 2);
tol = 1e-8;
maxit = 200;

%% Initial solution 
X = (A'*A + lambda*eye(M))\(A'*b);

if q == 2
    rho = norm(b - A*X);
    eta = norm(X);
    WR = eye(M);
    return;
end

X0 = X;

%% Iteration
epsr = 1e-3;
crit = 1;
ee = 2;
while (crit > tol) && (ee <= maxit)
    wr = (q/2)*max(abs(X), epsr).^(q-2);
    
    WR = sparse(1:M, 1:M, wr(:), M, M);
    
    X = (A'*A + lambda*WR)\(A'*b);
    
    crit = norm(X - X0, 1)/norm(X0, 1);
    X0 = X;
    
    ee = ee + 1;
end

rho = norm(b - A*X);
eta = norm(X, q);

end
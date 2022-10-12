function [F, M, res] = lcurve_tsvd(A, b)

[U, S, V] = svd(A);
Sd = diag(S);

Ns = length(Sd);

rho = zeros(Ns, 1);
eta = zeros(Ns, 1);

for ee = 1:Ns
    Ured = U(:, 1:ee);
    Sred = diag(Sd(1:ee));
    Vred = V(:, 1:ee);
    invSred = diag(1./Sd(1:ee));
    
    Ared = Ured*Sred*Vred';
    F = (Vred*invSred*Ured')*b;
    
    rho(ee) = norm(b - Ared*F);
    eta(ee) = norm(F);
end

% Curvature computation
k = kappa([rho, eta], 0);

% Number of singular values retained
M = find(k == max(k));

% Computation of the solution
Ured = U(:, 1:M);
Vred = V(:, 1:M);
invSred = diag(1./Sd(1:M));

F = (Vred*invSred*Ured')*b;

res.eta = eta;
res.rho = rho;
res.k = k;


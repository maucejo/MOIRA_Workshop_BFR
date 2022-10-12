function [X, taun, tauf] = RVR_L(A, b, L)
warning off
%Initialization
an = 1;
bn = 1e-18;
af = 1;
bf = 1e-18;
tol = 1e-2;
N = size(A, 1);

X = reg_addL(A, b, 2, L);
X0 = X;

crit = 1;
while crit > tol
    tauf = (2*af - 1)./(2*bf + abs(L*X).^2);
    tauf(tauf > 1e10) = 1e10;
    taun = (N + an - 1)/(bn + norm(b - A*X)^2);
    
    Lambda = diag(tauf/taun);
    X = (A'*A + L'*Lambda*L)\(A'*b);
    
    crit = norm(X - X0, 1)/norm(X0, 1);
    X0 = X;
end
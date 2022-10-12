function [X, tau_n, tau_f] = RVRL_sampling(A, b, L, Ns)
warning off
%Initialization
an = 1;
bn = 1e-18;
af = 1;
bf = 1e-18;
[M, N] = size(L);

tau_n = zeros(1, Ns);
tau_f = zeros(M, Ns);
X = zeros(N, Ns);

[X(:, 1), tau_n(1), tau_f(:, 1)] = RVR_L(A, b, L);

for ee = 2:Ns
    if mod(ee, 500) == 0
        disp([num2str(ee) '/' num2str(Ns)])
    end
    
    LX = L*X(:, ee - 1);
    for gg = 1:M
        tau_f(gg, ee) = gamrand(af + 0.5, bf + abs(LX(gg))^2);
    end
    
    tau_f(tau_f(:, ee) > 1e10, ee) = 1e10;
    tau_n(ee) = gamrand(an + N, bn + norm(b - A*X(:, ee - 1))^2);
    X(:, ee) = draw_F(A, b, L, tau_n(ee), tau_f(:, ee));   
end
end

%%
function X = draw_F(A, b, L, taun, tauf)
M = size(A, 2);

Tf = diag(tauf);

muX = (A'*A + L'*(Tf/taun)*L)\(A'*b);     
Sigma_f = (taun*(A'*A) + L'*Tf*L)\eye(M);

didx = (0:M-1)*M + (1:M); % Index of diagonal elements
Sigma_f(didx) = real(diag(Sigma_f));

% Compute Cholesky decomposition
[U, S] = svd(Sigma_f);
S = U*sqrt(S);

% Draw samples from standard normal
z = (randn(M, 1) + 1i*randn(M, 1))/sqrt(2);

X = muX + S*z;
end
function [X, tau_n, tau_f] = RVR_sampling(A, b, Ns)
warning off
%Initialization
an = 1;
bn = 1e-18;
af = 1;
bf = 1e-18;
[N, M] = size(A);

tau_n = zeros(1, Ns);
tau_f = zeros(M, Ns);
X = zeros(M, Ns);

[X(:, 1), tau_n(1), tau_f(:, 1)] = RVR(A, b);

for ee = 2:Ns
    if mod(ee, 500) == 0
        disp([num2str(ee) '/' num2str(Ns)])
    end
    for gg = 1:M
        tau_f(gg, ee) = gamrand(af + 0.5, bf + abs(X(gg, ee - 1))^2);
    end
    tau_n(ee) = gamrand(an + N, bn + norm(b - A*X(:, ee - 1))^2);
    X(:, ee) = draw_F(A, b, tau_n(ee), tau_f(:, ee));   
end
end

%%
function X = draw_F(A, b, taun, tauf)
M = size(A, 2);

Tf = diag(tauf);

muX = (A'*A + Tf/taun)\(A'*b);     
Sigma_f = (taun*(A'*A) + Tf)\eye(M);

didx = (0:M-1)*M + (1:M); % Index of diagonal elements
Sigma_f(didx) = real(diag(Sigma_f));

% Compute Cholesky decomposition
S = chol(Sigma_f, 'lower');

% Draw samples from standard normal
z = (randn(M, 1) + 1i*randn(M, 1))/sqrt(2);

X = muX + S*z;
end
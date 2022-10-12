function [Y, X] = MBFL_sampling(A, b, L, q, Nsamples)
% Minimal Bayesian Formulation
M = size(A, 2);

% Compute the mean
[X, lambda, W] = reg_addL(A, b, q, L);

% Estimate tau_n and tau_f
[tau_n, tau_f] = estimate_precisionL(A, b, sqrt(W)*L, lambda);

% Compute coviariance matrix
Sigma_f = (tau_n*(A'*A) + tau_f*L'*W*L)\eye(M);

didx = (0:M-1)*M + (1:M); % Index of diagonal elements
Sigma_f(didx) = real(diag(Sigma_f));

% Compute Cholesky decomposition
S = chol(Sigma_f, 'lower');

% Draw samples from standard normal
z = (randn(M, Nsamples) + 1i*randn(M, Nsamples))/sqrt(2);

X_samples = X + S*z;

Y = prctile(real(X_samples), [2.5, 50, 97.5], 2);
end
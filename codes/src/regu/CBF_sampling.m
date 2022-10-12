function [X, tau_n, tau_f, q] = CBF_sampling(A, b, q0, Ns)

%Initialization
an = 1;
bn = 1e-18;
af = 1;
bf = 1e-18;
[N, M] = size(A);

tau_n = zeros(1, Ns);
tau_f = zeros(1, Ns);
q = zeros(1, Ns);
X = zeros(M, Ns);

[X(:, 1), tau_n(1), tau_f(1), q(1)] = OBR(A, b, q0);

for ee = 2:Ns
    if mod(ee, 500) == 0
        disp([num2str(ee) '/' num2str(Ns)])
    end
    q(ee) = draw_q(q(ee - 1), X(:, ee - 1), tau_f(ee - 1));
    tau_f(ee) = gamrand(af + M/q(ee), bf + norm(X(:, ee - 1), q(ee))^q(ee));
    tau_n(ee) = gamrand(an + N, bn + norm(b - A*X(:, ee - 1))^2);
    X(:, ee) = draw_F(A, b, tau_n(ee), tau_f(ee), q(ee), X(:, ee - 1));   
end
end

%%
function X = draw_F(A, b, taun, tauf, q, X0)
M = size(A, 2);

[muX, W] = reg_(A, b, q, tauf/taun, X0);
     
Sigma_f = (taun*(A'*A) + tauf*W)\eye(M);

didx = (0:M-1)*M + (1:M); % Index of diagonal elements
Sigma_f(didx) = real(diag(Sigma_f));

% Compute Cholesky decomposition
S = chol(Sigma_f, 'lower');

% Draw samples from standard normal
z = (randn(M, 1) + 1i*randn(M, 1))/sqrt(2);

X = muX + S*z;
end

%%
function x = draw_q(x0, Y, tau)
% Initialisation
alpha_p = 1;
beta_p = 1e-18;
lb = 0.05;
ub = 2.05;
delta = 5e-3;
L = 20;
M = length(Y);

U = @(x) M*log(gamma(1/x)) - M*log(tau)/x - (M + alpha_p - 1)*log(x) + beta_p*x + tau*sum(abs(Y).^x);
dU = @(x) beta_p + tau*sum(log(abs(Y)).*abs(Y).^x) - (alpha_p + M - 1)/x + (M*log(tau) - M*psi(1/x))/x^2;

K = @(q) q.'*q/2;

%% Sampling - HMC

% SAMPLE RANDOM MOMENTUM
p0 = randn;

xStar = x0;
pStar = p0;

% SIMULATE HAMILTONIAN DYNAMICS
% FIRST 1/2 STEP OF MOMENTUM
pStar = pStar - delta/2*dU(xStar);

% FULL STEPS
for jL = 1:L - 1
    % POSITION/SAMPLE
    xStar = xStar + delta*pStar;
    
    if (xStar > ub) || (xStar < lb)
        while (xStar > ub) || (xStar < lb)
            if xStar > ub
                xStar = ub - (xStar - ub);
                pStar = -pStar;
            elseif xStar < lb
                xStar = lb + (lb - xStar);
                pStar = -pStar;
            end
        end
        
    else
        %  MOMENTUM
        pStar = pStar - delta*dU(xStar);
    end
end

% LAST HALP STEP
pStar = pStar - delta/2*dU(xStar);

% Negate momentum at end of trajectory to make the proposal symmetric
pStar = -pStar;

% EVALUATE ENERGIES AT
% START AND END OF TRAJECTORY
U0 = U(x0);
UStar = U(xStar);

K0 = K(p0);
KStar = K(pStar);

% ACCEPTANCE/REJECTION CRITERION
% alpha = min(1, exp((U0 + K0) - (UStar + KStar)));
log_alpha = min(0, (U0 + K0) - (UStar + KStar));

% u = rand;
log_u = log(rand);
if log_u < log_alpha
    x = xStar;
else
    x = x0;
end
end
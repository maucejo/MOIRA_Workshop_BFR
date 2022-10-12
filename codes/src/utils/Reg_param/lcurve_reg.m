function [X, M, res] = lcurve_reg(A, b, q, npoints)

smin_ratio = 16*eps;  % Smallest regularization parameter
s = svd(A);

reg_param = zeros(npoints, 1);

reg_param(npoints) = max([s(end), s(1)*smin_ratio]);
ratio = (s(1)/reg_param(npoints))^(1/(npoints-1));

for ee = npoints-1:-1:1
    reg_param(ee) = ratio*reg_param(ee+1);
end

lambda = sort(reg_param.^2);

rho = zeros(npoints, 1);
eta = zeros(npoints, 1);
Xtemp = zeros(size(A, 2), npoints);

for ee = 1:npoints
   if mod(ee, 50) == 0
       disp([num2str(ee) '/' num2str(npoints)])
   end
   [Xtemp(:, ee), rho(ee), eta(ee)] = reg_add_lambda(A, b, q, lambda(ee)); 
end

% Computation of the curvature
k = kappa([rho, eta], 0);

% Index of the maximum curvature of the L-curve
M = find(k == max(k));

X = Xtemp(:, M);

res.lambda = lambda;
res.rho = rho;
res.eta = eta;
res.k = k;



    
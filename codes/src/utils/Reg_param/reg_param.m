function lambda = reg_param(U, sm, b, method, disp)

% lambda = reg_param(U, sm, b, method, disp) calcule le parametre de
% regularisation optimal
%
% Entrees:
%   - U: Vecteur singulier gauche
%   - sm: Vecteur des valeurs singulieres
%   - b: Vecteur des donnees mesurees ou simulees
%   - method: Methode de calcul du parametre de regularisation
%             * 1: L-curve
%             * 2: Estimateur bayesien
%             * 3: GCV
%
% Sortie:
%   - lambda: parametre de regularisation optimal

npoints = 200;  % Number of points on the L-curve
smin_ratio = 16*eps;  % Smallest regularization parameter

[m, n] = size(U);
[p, ps] = size(sm);

y = U'*b;
y2 = norm(b)^2 - norm(y)^2;

if (ps==1)
    s = sm;
    y = y(1:p);
else
    s = sm(p:-1:1,1)./sm(p:-1:1,2);
    y = y(p:-1:1);
end

xi = y(1:p)./s;

reg_param = zeros(npoints, 1);
s2 = s.^2;

reg_param(npoints) = max([s(p), s(1)*smin_ratio]);
ratio = (s(1)/reg_param(npoints))^(1/(npoints-1));

for ee = npoints-1:-1:1
    reg_param(ee) = ratio*reg_param(ee+1);
end

if method == 1 % Lcurve
    eta = zeros(npoints, 1);
    rho = zeros(npoints, 1);

    for ee = 1:npoints
        f = s2./(s2 + reg_param(ee)^2);
        eta(ee) = norm(f.*xi);
        rho(ee) = norm((1-f).*y(1:p));
    end

    if (m > n && y2 > 0)
        rho = sqrt(rho.^2 + y2);
    end

    [reg_corner, rho_c, eta_c] = l_corner(rho, eta, reg_param, U, sm, b, 'Tikh');

    lambda = reg_corner^2;

    if disp == 1
        figure
        loglog(rho, eta)
        hold on
        loglog(rho_c, eta_c, 'rx')
        xlabel('Residual norm ||b - Ax ||_2')
        ylabel('Solution semi-norm || Lx ||_2')
    end

elseif method == 2 % Bayes

    Jmap = zeros(npoints, 1);

    M = length(y);

    for ee = 1:npoints
        alpha2 = sum(abs(y).^2./(s2 + reg_param(ee)^2))/M;
        a = sum(log(s2 + reg_param(ee)^2));
        Jmap(ee) = a + (M - 2)*log(alpha2);
    end

    f = find(Jmap == min(Jmap));

    lambda = reg_param(f(end))^2;

    if disp == 1
       figure
       plot(reg_param.^2, Jmap)
       hold on
       plot(lambda, Jmap(f(end)), 'rx')
       xlabel('\lambda')
       ylabel('J_m_a_p')
    end

else % GCV

    Jgcv = zeros(npoints, 1);

    delta0 = 0;
    if (m > n && y2 > 0)
        delta0 = y2;
    end

    for ee = 1:npoints
        al = reg_param(ee)^2./(s2 + reg_param(ee)^2);
        Jgcv(ee) = (norm(al.*y)^2 + delta0)/(m-n + sum(al)^2);
    end

    f = find(Jgcv == min(Jgcv));

    lambda = reg_param(f(end))^2;

    if disp == 1
       figure
       plot(reg_param.^2, Jgcv)
       hold on
       plot(lambda, Jgcv(f(end)), 'rx')
       xlabel('\lambda')
       ylabel('J_g_c_v')
    end
end

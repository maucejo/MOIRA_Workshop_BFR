%% Beam parameters
Lp = 1;         % Length of the beam [m]
b = 0.03;       % Width of the beam [m]
h = 0.01;       % Thickness of the beam [m]
S = b*h;        % Area of the beam section (IPN) [m^2]
I = b*h^3/12;   % Second moment of area [m^4]

E = 2.1e11;     % Young modulus (acier) [Pa]
rho = 7850;     % Density [kg/m^3]
eta = 1e-2;     % Damping factor [dimensionless]

freq = 350;      % Frequency of interest [Hz]

% Mesh definition
Nelem = 20;
[K, M, bmesh] = beam_matrix(0, Lp, S, I, E, rho, Nelem, 2);

Nexc = 13;                     % ID of the excitation node
Xexc = bmesh.Nodes(Nexc, 2);   % Coordinate of the excitation point
Xi = bmesh.Nodes(:,2);         % Reconstruction mesh
Ni = length(bmesh.Nodes(:,2)); % Number of reconstruction points

% Signal-to-Noise ratio of the Gaussian white noise 
SNR = 20; 

%% Generation of measured noisy data
load('data_point_force.mat');

rng('default');
X = agwn(Xref.', SNR).';

%% Visualization
F = zeros(Ni, 1);
F(Nexc) = 1;

figure(1)
plot(bmesh.Nodes(:, 2), F, 'b', 'linewidth', 2)
xlabel('x (m)')
ylabel('Force (N)')

figure(2)
subplot(2, 1, 1)
plot(bmesh.Nodes(:, 2), real(Xref), 'b', 'linewidth', 2, 'displayname', 'Reference')
hold on
plot(bmesh.Nodes(:, 2), real(X), '--r', 'linewidth', 2, 'displayname', 'Noisy')
hold off
ylabel('Real part - Acc. (m/s^2)')
legend

subplot(2, 1, 2)
plot(bmesh.Nodes(:, 2), imag(Xref), 'b', 'linewidth', 2, 'displayname', 'Reference')
hold on
plot(bmesh.Nodes(:, 2), imag(X), '--r', 'linewidth', 2, 'displayname', 'Noisy')
hold off
xlabel('x (m)')
ylabel('Imag. part - Acc. (m/s^2)')

%% Computation of the transfer function matrix
Ndof = size(K, 1);
pos_dof_mes = 1:2:Ndof;    % Measured dofs (transverse displacement motion)
pos_dof_nmes = 2:2:Ndof;   % Rotation of the section 

Dc = dynamic_condensation(K*(1+1i*eta), M, freq, pos_dof_mes, pos_dof_nmes);
H = -(2*pi*freq)^2*(Dc\eye(Ni));

%% Naive reconstruction
Fnaive = H\X;

figure(3)
plot(bmesh.Nodes(:, 2), F, 'b', 'linewidth', 2, 'displayname', 'Reference')
hold on
plot(bmesh.Nodes(:, 2), real(Fnaive) ,'--r', 'linewidth', 2, 'displayname', 'Naive')
hold off
xlabel('x (m)')
ylabel('Real part - Force (N)')
axis tight
legend

%% OLS
Fols = (H'*H)\(H'*X);

figure(4)
plot(bmesh.Nodes(:, 2), F, 'b', 'linewidth', 2, 'displayname', 'Reference')
hold on
plot(bmesh.Nodes(:, 2), real(Fols) ,'--r', 'linewidth', 2, 'displayname', 'OLS')
hold off
xlabel('x (m)')
ylabel('Real part - Force (N)')
axis tight
legend

%% TSVD
[Ftsvd, Ns, res] = lcurve_tsvd(H, X);

figure(5)
plot(bmesh.Nodes(:, 2), F, 'b', 'linewidth', 2, 'displayname', 'Reference')
hold on
plot(bmesh.Nodes(:, 2), real(Ftsvd) ,'--r', 'linewidth', 2, 'displayname', 'TSVD')
hold off
xlabel('x (m)')
ylabel('Real part - Force (N)')
axis tight
legend

figure(6)
plot(res.rho, res.eta, 'b', 'linewidth', 2)
hold on
plot(res.rho(Ns), res.eta(Ns), 'rd', 'linewidth', 3, 'markersize', 6)
xlabel('Residual norm')
ylabel('Solution norm')
hold off
axis tight

figure(7)
plot(1:Ni, res.k, 'b', 'linewidth', 2)
hold on
plot(Ns, res.k(Ns), 'rd', 'linewidth', 3, 'markersize', 6)
xlabel('Singular value index')
ylabel('Curvature')
xlim([1 Ni])
hold off

ftsvd = filter_factor_tsvd(H, Ns);

%% Tikhonov regularization
npoints = 200;
[Ftik, pos_tik, res_tik] = lcurve_reg(H, X, 2, npoints);

figure(8)
plot(bmesh.Nodes(:, 2), F, 'b', 'linewidth', 2, 'displayname', 'Reference')
hold on
plot(bmesh.Nodes(:, 2), real(Ftik) ,'--r', 'linewidth', 2, 'displayname', 'l2-reg.')
hold off
xlabel('x (m)')
ylabel('Real part - Force (N)')
axis tight
legend

figure(9)
plot(res_tik.rho, res_tik.eta, 'b', 'linewidth', 2)
hold on
plot(res_tik.rho(pos_tik), res_tik.eta(pos_tik), 'rd', 'linewidth', 3, 'markersize', 6)
xlabel('Residual norm')
ylabel('Solution norm')
hold off
axis tight

figure(10)
semilogx(res_tik.lambda, res_tik.k, 'b', 'linewidth', 2)
hold on
semilogx(res_tik.lambda(pos_tik), res_tik.k(pos_tik), 'rd', 'linewidth', 3, 'markersize', 6)
xlabel('Regularization parameter value')
ylabel('Curvature')
ylim([min(res_tik.k), max(res_tik.k) + 0.5])
xlim([min(res_tik.lambda), max(res_tik.lambda)])
hold off

ftik = filter_factor_Tik(H, res_tik.lambda(pos_tik));

%% LASSO regularization
npoints = 200;
[Fl1, pos_l1, res_l1] = lcurve_reg(H, X, 1, npoints);

figure(11)
plot(bmesh.Nodes(:, 2), F, 'b', 'linewidth', 2, 'displayname', 'Reference')
hold on
plot(bmesh.Nodes(:, 2), real(Fl1) ,'--r', 'linewidth', 2, 'displayname', 'l1-reg.')
hold off
xlabel('x (m)')
ylabel('Real part - Force (N)')
axis tight
legend

figure(12)
plot(res_l1.rho, res_l1.eta, 'b', 'linewidth', 2)
hold on
plot(res_l1.rho(pos_l1), res_l1.eta(pos_l1), 'rd', 'linewidth', 3, 'markersize', 6)
xlabel('Residual norm')
ylabel('Solution norm')
hold off
axis tight

figure(13)
semilogx(res_l1.lambda, res_l1.k, 'b', 'linewidth', 2)
hold on
semilogx(res_l1.lambda(pos_l1), res_l1.k(pos_l1), 'rd', 'linewidth', 3, 'markersize', 6)
xlabel('Regularization parameter value')
ylabel('Curvature')
ylim([min(res_l1.k), max(res_l1.k) + 0.5])
xlim([min(res_l1.lambda), max(res_l1.lambda)])
hold off

[~, ~, ~, WR] = reg_add_lambda(H, X, 1, res_l1.lambda(pos_l1));
fl1 = filter_factor_lp(H, sqrt(WR), res_l1.lambda(pos_l1));

%% lp regularization
npoints = 500;
p = 0.5;
[Flp, pos_lp, res_lp] = lcurve_reg(H, X, p, npoints);

figure(12)
plot(bmesh.Nodes(:, 2), F, 'b', 'linewidth', 2, 'displayname', 'Reference')
hold on
plot(bmesh.Nodes(:, 2), real(Flp) ,'--r', 'linewidth', 2, 'displayname', 'lp-reg.')
hold off
xlabel('x (m)')
ylabel('Real part - Force (N)')
axis tight
legend

figure(13)
plot(res_lp.rho, res_lp.eta, 'b', 'linewidth', 2)
hold on
plot(res_lp.rho(pos_lp), res_lp.eta(pos_lp), 'rd', 'linewidth', 3, 'markersize', 6)
xlabel('Residual norm')
ylabel('Solution norm')
hold off
axis tight

figure(14)
semilogx(res_lp.lambda, res_lp.k, 'b', 'linewidth', 2)
hold on
semilogx(res_lp.lambda(pos_lp), res_lp.k(pos_lp), 'rd', 'linewidth', 3, 'markersize', 6)
xlabel('Regularization parameter value')
ylabel('Curvature')
ylim([min(res_lp.k), max(res_lp.k) + 0.5])
xlim([min(res_lp.lambda), max(res_lp.lambda)])
hold off

[~, ~, ~, WR] = reg_add_lambda(H, X, p, res_lp.lambda(pos_lp));
flp = filter_factor_lp(H, sqrt(WR), res_lp.lambda(pos_lp));

%% MBF - Minimal Bayesian Formulation
qs = 2;
Nsamples = 1e4;

[Ymbf, Fmbf] = MBF_sampling(H, X, qs, Nsamples);

X1 = [bmesh.Nodes(:,2).', fliplr(bmesh.Nodes(:,2).')];
Y1 = [Ymbf(:,1).', fliplr(Ymbf(:,3).')];

figure(15)
h1 = fill(X1, Y1, [225, 225, 225]/255, 'Edgecolor', 'none');
hold on
h2 = plot(bmesh.Nodes(:, 2), F, 'b', 'linewidth', 2);
h3 = plot(bmesh.Nodes(:,2), real(Ymbf(:, 2)), '--r', 'linewidth', 2);
hold off
xlabel('x (m)')
ylabel('Real part - Force (N)')
axis tight
legend([h2 h3 h1], {'Reference', 'Median / MAP', '95% CI'})

%% OBR
q0 = 0.5;
[Fobr, tau_nobr, tau_fobr, q_obr] = OBR(H, X, q0);

figure(16)
plot(bmesh.Nodes(:, 2), F, 'b', 'linewidth', 2, 'displayname', 'Reference')
hold on
plot(bmesh.Nodes(:, 2), real(Fobr) ,'--r', 'linewidth', 2, 'displayname', 'MAP')
hold off
xlabel('x (m)')
ylabel('Real part - Force (N)')
axis tight
legend

%% CBF
rng(123)
Nsamples = 1e5;
Nburn = round(Nsamples/2);
q0 = 2;
[Fs, tau_ns, tau_fs, qs] = CBF_sampling(H, X, q0, Nsamples);

% Diagnostic
res_taun = raftery_lewis(tau_ns, 0.025, 0.01, 0.95, 0);
res_tauf = raftery_lewis(tau_fs, 0.025, 0.01, 0.95, 0);
res_q = raftery_lewis(qs, 0.025, 0.01, 0.95, 0);
res_F = raftery_lewis(real(Fs(Nexc, :)), 0.025, 0.01, 0.95, 0);

Fstats = real(Fs(:, Nburn:end));
taun_stats = tau_ns(Nburn:end);
tauf_stats = tau_fs(Nburn:end);
q_stats = qs(Nburn:end);

% Visualization
[pdf_taun, xn] = kernelDensity(taun_stats);
Ytn = prctile(taun_stats, [2.5, 50, 97.5]);
posn = (xn >= Ytn(1)) & (xn <= Ytn(3));

[pdf_tauf, xf] = kernelDensity(tauf_stats);
Ytf = prctile(tauf_stats, [2.5, 50, 97.5]);
posf = (xf >= Ytf(1)) & (xf <= Ytf(3));

[pdf_q, xq] = kernelDensity(q_stats);
Yq = prctile(q_stats, [2.5, 50, 97.5]);
posq = (xq >= Yq(1)) & (xq <= Yq(3));

[pdf_F, xF] = kernelDensity(Fstats(Nexc,:));
Yf = prctile(real(Fstats), [2.5, 50, 97.5], 2);
posF = (xF >= Yf(Nexc, 1)) & (xF <= Yf(Nexc, 3));

figure(17)
tiledlayout(2,2, 'TileSpacing', 'compact')

nexttile
area(xn(posn), pdf_taun(posn), 'facecolor', [225, 225, 225]/255, 'facealpha', 0.75, 'linestyle', 'none')
hold on
plot(xn, pdf_taun, 'b', 'linewidth', 2)
xlabel('\tau_n')
ylabel('Normalized density')
axis tight

nexttile
area(xf(posf), pdf_tauf(posf), 'facecolor', [225, 225, 225]/255, 'facealpha', 0.75, 'linestyle', 'none')
hold on
plot(xf, pdf_tauf, 'b', 'linewidth', 2)
xlabel('\tau_f')
ylabel('Normalized density')
axis tight

nexttile
area(xq(posq), pdf_q(posq), 'facecolor', [225, 225, 225]/255, 'facealpha', 0.75, 'linestyle', 'none')
hold on
plot(xq, pdf_q, 'b', 'linewidth', 2)
xlabel('q')
ylabel('Normalized density')
axis tight

nexttile
area(xF(posF), pdf_F(posF), 'facecolor', [225, 225, 225]/255, 'facealpha', 0.75, 'linestyle', 'none')
hold on
plot(xF, pdf_F, 'b', 'linewidth', 2)
xlabel('Force amplitude')
ylabel('Normalized density')
axis tight

% Markov chains
figure(18)
tiledlayout(2,2, 'TileSpacing', 'compact')

nexttile
plot(tau_ns, 'b')
ylabel('\tau_n')
axis tight

nexttile
plot(tau_fs, 'b')
ylabel('\tau_f')
axis tight

nexttile
plot(qs, 'b')
ylabel('q')
xlabel('Samples ID')
axis tight

nexttile
plot(real(Fs(Nexc, :)), 'b')
ylabel('Force amplitude')
xlabel('Samples ID')
axis tight

% Excitation field
X1 = [bmesh.Nodes(:,2).', fliplr(bmesh.Nodes(:,2).')];
Y1 = [Yf(:,1).', fliplr(Yf(:,3).')];

figure(19)
h1 = fill(X1, Y1, [225, 225, 225]/255, 'Edgecolor', 'none');
hold on
h2 = plot(bmesh.Nodes(:, 2), F, 'b', 'linewidth', 2);
h3 = plot(bmesh.Nodes(:,2), real(Yf(:, 2)), '--r', 'linewidth', 2);
hold off
xlabel('x (m)')
ylabel('Real part - Force (N)')
axis tight
legend([h2 h3 h1], {'Reference', 'Median / MAP', '95% CI'})

%% RVR
[Frvr, tau_nrvr, tau_frvr] = RVR(H, X);

figure(20)
plot(bmesh.Nodes(:, 2), F, 'b', 'linewidth', 2, 'displayname', 'Reference')
hold on
plot(bmesh.Nodes(:, 2), real(Frvr) ,'--r', 'linewidth', 2, 'displayname', 'MAP')
hold off
xlabel('x (m)')
ylabel('Real part - Force (N)')
axis tight
legend

%% RVR Sampling
rng(123)
Nsamples = 1e5;
Nburn = round(Nsamples/2);

[Fs_rvr, tauns_rvr, taufs_rvr] = RVR_sampling(H, X, Nsamples);

% Diagnostic
res_taun_rvm = raftery_lewis(tauns_rvr, 0.025, 0.01, 0.95, 0);
res_tauf_rvm = raftery_lewis(taufs_rvr(Nexc, :), 0.025, 0.01, 0.95, 0);
res_F_rvm = raftery_lewis(real(Fs_rvr(Nexc, :)), 0.025, 0.01, 0.95, 0);

Fstats_rvm = real(Fs_rvr(:, Nburn:end));
taun_stats_rvm = tauns_rvr(Nburn:end);
tauf_stats_rvm = taufs_rvr(:, Nburn:end);

% Visualization
[pdf_taun, xn] = kernelDensity(taun_stats_rvm);
Ytn = prctile(taun_stats_rvm, [2.5, 50, 97.5]);
posn = (xn >= Ytn(1)) & (xn <= Ytn(3));

[pdf_tauf, xf] = kernelDensity(tauf_stats_rvm(Nexc, :));
Ytf = prctile(tauf_stats_rvm, [2.5, 50, 97.5], 2);
posf = (xf >= Ytf(Nexc, 1)) & (xf <= Ytf(Nexc, 3));

[pdf_F, xF] = kernelDensity(Fstats_rvm(Nexc, :));
Yf = prctile(real(Fstats_rvm), [2.5, 50, 97.5], 2);
posF = (xF >= Yf(Nexc, 1)) & (xF <= Yf(Nexc, 3));

figure(21)
tiledlayout(2,2, 'TileSpacing', 'compact')

nexttile
area(xn(posn), pdf_taun(posn), 'facecolor', [225, 225, 225]/255, 'facealpha', 0.75, 'linestyle', 'none')
hold on
plot(xn, pdf_taun, 'b', 'linewidth', 2)
xlabel('\tau_n')
ylabel('Normalized density')
axis tight

nexttile
area(xf(posf), pdf_tauf(posf), 'facecolor', [225, 225, 225]/255, 'facealpha', 0.75, 'linestyle', 'none')
hold on
plot(xf, pdf_tauf, 'b', 'linewidth', 2)
xlabel('\tau_f')
ylabel('Normalized density')
axis tight

nexttile
area(xF(posF), pdf_F(posF), 'facecolor', [225, 225, 225]/255, 'facealpha', 0.75, 'linestyle', 'none')
hold on
plot(xF, pdf_F, 'b', 'linewidth', 2)
xlabel('Force amplitude')
ylabel('Normalized density')
axis tight

% Markov chains
figure(22)
tiledlayout(2,2, 'TileSpacing', 'compact')

nexttile
plot(tauns_rvr, 'b')
ylabel('\tau_n')
axis tight

nexttile
plot(taufs_rvr(Nexc, :), 'b')
ylabel('\tau_f')
axis tight

nexttile
plot(real(Fs_rvr(Nexc, :)), 'b')
ylabel('Force amplitude')
xlabel('Samples ID')
axis tight

% Excitation field
X1 = [bmesh.Nodes(:,2).', fliplr(bmesh.Nodes(:,2).')];
Y1 = [Yf(:,1).', fliplr(Yf(:,3).')];

figure(23)
h1 = fill(X1, Y1, [225, 225, 225]/255, 'Edgecolor', 'none');
hold on
h2 = plot(bmesh.Nodes(:, 2), F, 'b', 'linewidth', 2);
h3 = plot(bmesh.Nodes(:,2), real(Yf(:, 2)), '--r', 'linewidth', 2);
hold off
xlabel('x (m)')
ylabel('Real part - Force (N)')
axis tight
legend([h2 h3 h1], {'Reference', 'Median / MAP', '95% CI'})
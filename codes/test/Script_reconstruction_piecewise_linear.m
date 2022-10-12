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

Nexc = 12:17;                     % ID of the excitation node
Xexc = bmesh.Nodes(Nexc, 2);   % Coordinate of the excitation point
Xi = bmesh.Nodes(:,2);         % Reconstruction mesh
Ni = length(bmesh.Nodes(:,2)); % Number of reconstruction points

% Signal-to-Noise ratio of the Gaussian white noise 
SNR = 20; 
%% Generation of measured noisy data
load('data_piecewise_linear.mat');

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

%% OBR
q0 = 2;
[Fobr, tau_nobr, tau_fobr, q_obr] = OBR(H, X, q0);

figure(3)
plot(bmesh.Nodes(:, 2), F, 'b', 'linewidth', 2, 'displayname', 'Reference')
hold on
plot(bmesh.Nodes(:, 2), real(Fobr) ,'--r', 'linewidth', 2, 'displayname', 'MAP')
hold off
xlabel('x (m)')
ylabel('Real part - Force (N)')
axis tight
legend

%% RVR
[Frvr, tau_nrvr, tau_frvr] = RVR(H, X);

figure(4)
plot(bmesh.Nodes(:, 2), F, 'b', 'linewidth', 2, 'displayname', 'Reference')
hold on
plot(bmesh.Nodes(:, 2), real(Frvr) ,'--r', 'linewidth', 2, 'displayname', 'MAP')
hold off
xlabel('x (m)')
ylabel('Real part - Force (N)')
axis tight
legend

%% MBF - Minimal Bayesian Formulation
qs = 0.2;
Nsamples = 1e5;

Lx = get_l(Ni, 1);
[Ymbf, Fmbf] = MBFL_sampling(H, X, Lx, qs, Nsamples);

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
q0 = 1;
Lx = get_l(Ni, 1);
[FobrL, tau_nobrL, tau_fobrL, q_obrL] = OBR_L(H, X, q0, Lx);

figure(5)
plot(bmesh.Nodes(:, 2), F, 'b', 'linewidth', 2, 'displayname', 'Reference')
hold on
plot(bmesh.Nodes(:, 2), real(FobrL) ,'--r', 'linewidth', 2, 'displayname', 'MAP')
hold off
xlabel('x (m)')
ylabel('Real part - Force (N)')
axis tight
legend

%% CBF
rng(123)
Nsamples = 1e5;
Nburn = round(Nsamples/2);

q0 = 1;
Lx = get_l(Ni, 1);
[Fs, tau_ns, tau_fs, qs] = CBFL_sampling(H, X, Lx, q0, Nsamples);

% Diagnostic
res_taun = raftery_lewis(tau_ns, 0.025, 0.01, 0.95, 0);
res_tauf = raftery_lewis(tau_fs, 0.025, 0.01, 0.95, 0);
res_q = raftery_lewis(qs, 0.025, 0.01, 0.95, 0);

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

% Excitation field
Yf = prctile(real(Fstats), [2.5, 50, 97.5], 2);
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
Lx = get_l(Ni, 1);
[FrvrL, tau_nrvrL, tau_frvrL] = RVR_L(H, X, Lx);

figure(6)
plot(bmesh.Nodes(:, 2), F, 'b', 'linewidth', 2, 'displayname', 'Reference')
hold on
plot(bmesh.Nodes(:, 2), real(FrvrL) ,'--r', 'linewidth', 2, 'displayname', 'MAP')
hold off
xlabel('x (m)')
ylabel('Real part - Force (N)')
axis tight
legend

%% RVR Sampling
rng(123)
Nsamples = 1e5;
Nburn = round(Nsamples/2);

Lx = get_l(Ni, 1);
[Fs_rvrL, tauns_rvrL, taufs_rvrL] = RVRL_sampling(H, X, Lx, Nsamples);

% Diagnostic
res_taun_rvm = raftery_lewis(tauns_rvrL, 0.025, 0.01, 0.95, 0);

Fstats_rvm = real(Fs_rvrL(:, Nburn:end));
taun_stats_rvm = tauns_rvrL(Nburn:end);
tauf_stats_rvm = taufs_rvrL(:, Nburn:end);

% Visualization
[pdf_taun, xn] = kernelDensity(taun_stats_rvm);
Ytn = prctile(taun_stats_rvm, [2.5, 50, 97.5]);
posn = (xn >= Ytn(1)) & (xn <= Ytn(3));

figure(21)
area(xn(posn), pdf_taun(posn), 'facecolor', [225, 225, 225]/255, 'facealpha', 0.75, 'linestyle', 'none')
hold on
plot(xn, pdf_taun, 'b', 'linewidth', 2)
xlabel('\tau_n')
ylabel('Normalized density')
axis tight

% Markov chains
figure(22)
plot(tauns_rvrL, 'b')
ylabel('\tau_n')
axis tight

% Excitation field
Yf = prctile(real(Fstats_rvm), [2.5, 50, 97.5], 2);
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
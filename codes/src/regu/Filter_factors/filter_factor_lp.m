function f = filter_factor_lp(A, B, lambda)

[~,~,~, S, T] = gsvd(full(A), full(B), 0);
Ns = size(S, 1);

s = diag(S);
t = diag(T);

[s, pos] = sort(s, 'descend');
t = t(pos);

gam = s./t;

f = gam.^2./(gam.^2 + lambda);

figure
colororder({'k','b'})
yyaxis left
% plot(1:Ns, 1./s, 'ko', 'linewidth', 2, 'markersize', 6, 'displayname', '1/\sigma_i')
semilogy(1:Ns, 1./s, 'ko', 'linewidth', 2, 'markersize', 6, 'displayname', '1/\sigma_i')
hold on
% plot(1:Ns, f./s, 'rd', 'linewidth', 2, 'markersize', 10, 'displayname', 'f_i/\sigma_i')
semilogy(1:Ns, f./s, 'rd', 'linewidth', 2, 'markersize', 10, 'displayname', 'f_i/\sigma_i')
% ylim([-5, max(1./s)])
xlim([0.5, Ns + 0.5])
xlabel('Singular value index')
ylabel('Apparent inverse singular value')
legend('Location', 'northwest')

yyaxis right
plot(1:Ns, f, 'b', 'linewidth', 2, 'displayname', 'f_i')
ylim([-1e-3, 1+1e-3])
ylabel('Filter factor')
function f = filter_factor_tsvd(A, M)

s = svd(A);
Ns = length(s);

f = zeros(Ns, 1);
f(1:M) = 1;

figure
colororder({'k','b'})
yyaxis left
semilogy(1:Ns, 1./s, 'ko', 'linewidth', 2, 'markersize', 6, 'displayname', '1/\sigma_i')
hold on
semilogy(1:Ns, f./s, 'rd', 'linewidth', 2, 'markersize', 10, 'displayname', 'f_i/\sigma_i')
ylim([0, max(1./s)])
xlim([0.5, Ns + 0.5])
xlabel('Singular value index')
ylabel('Apparent inverse singular value')
legend('Location', 'northwest')

yyaxis right
plot(1:Ns, f, 'b', 'linewidth', 2, 'displayname', 'f_i')
ylim([-1e-3, 1+1e-3])
ylabel('Filter factor')
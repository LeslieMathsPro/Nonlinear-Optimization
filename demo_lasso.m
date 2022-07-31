clear;
seed = 97006855;
ss = RandStream('mt19937ar','Seed',seed);
RandStream.setGlobalStream(ss);

m = 512;
n = 1024;
A = randn(m, n);
u = sprandn(n, 1, 0.1);
b = A * u;

mu = 1;
x0 = u + 1e-1 * randn(n, 1);

% first find the f^\star
opts_g = struct();
opts_g.maxit = 50000;
opts_g.alpha0 = 1e-6;
opts_g.step_type = 'constant_step_size';
[x_g, out_g] = l1_subgrad(u, A, b, mu, opts_g);
f_star = out_g.f_hist_best(end);

% constant_step_size 
opts = struct();
opts.maxit = 3000;
opts.alpha0 = 0.0002;
opts.step_type = 'constant_step_size';
[x, out] = l1_subgrad(x0, A, b, mu, opts);
data1 = (out.f_hist_best - f_star) / f_star;

% constant_length_step
opts.gamma = 0.005;
opts.step_type = 'constant_step_length';
[x, out] = l1_subgrad(x0, A, b, mu, opts);
data2 = (out.f_hist_best - f_star) / f_star;

% square_summable_but_not_summable
opts.ssbns_a = 0.1;
opts.ssbns_b = 100;
opts.step_type = 'square_summable_but_not_summable';
[x, out] = l1_subgrad(x0, A, b, mu, opts);
data3 = (out.f_hist_best - f_star) / f_star;

% nonsummable diminishing
opts.nd_a = 0.002;
opts.step_type = 'nonsummable_diminishing';
[x, out] = l1_subgrad(x0, A, b, mu, opts);
data4 = (out.f_hist_best - f_star) / f_star;

% nonsummable diminishing step lengths
opts.ndsl_a = 10;
opts.step_type = 'nonsummable_diminishing_step_lengths';
[x, out] = l1_subgrad(x0, A, b, mu, opts);
data5 = (out.f_hist_best - f_star) / f_star;

%Polyak's step size
opts.f_star = f_star;
opts.step_type = 'polyak_step_size';
[x, out] = l1_subgrad(x0, A, b, mu, opts);
data6 = (out.f_hist_best - f_star) / f_star;

fig = figure;
semilogy(0:length(data1)-1, data1, '-', 'Color', [0.2 0.1 0.99], 'LineWidth', 2);
hold on
semilogy(0:length(data2)-1, data2, '--', 'Color', [0.99 0.1 0.2], 'LineWidth', 1.2);
hold on;
semilogy(0:length(data3)-1, data3, '-.', 'Color', [0.99 0.1 0.99], 'LineWidth', 1.5);
hold on;
semilogy(0:length(data4)-1, data4, ':', 'Color', [0.5 0.2 0.1], 'LineWidth', 1.8);
hold on;
semilogy(0:length(data5)-1, data5, '-', 'Color', [0.1 0.2 0.5], 'LineWidth', 1.9);
hold on;
semilogy(0:length(data6)-1, data6, '--', 'Color', [0.2 0.1 0.2], 'LineWidth', 1.4);
legend('$\alpha_k = 0.0002$', '$\alpha_k = 0.001/||g_k||_2$', '$\alpha_k = 0.1/(100+k)$', '$\alpha_k = 0.002/\sqrt{k}$', '$\alpha_k = 1/(k\times||g_k||_2)$', '$\alpha_k = \frac{f^{best,k} - f^\star}{||g_k||^2_2}$', 'interpreter', 'latex');
ylabel('$\hat{f}(x^k) - f^*)/f^*$', 'fontsize', 14, 'interpreter', 'latex');
xlabel('number of iterations');
print(fig, '-depsc','f.eps');
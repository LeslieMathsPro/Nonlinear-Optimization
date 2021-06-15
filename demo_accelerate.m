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
% opts = struct();
% opts.maxit = 100000;
% opts.alpha0 = 1e-6;
% opts.gamma = 0.005;
% opts.step_type = 'constant_step_length';
% [x, out] = l1_heavy_ball(x0, A, b, mu, opts);
% f_star = out.f_hist_best(end);
cvx_begin quiet
    cvx_precision
    variable x(n)
    minimize(0.5*sum_square(A*x - b) + mu*norm(x,1))
cvx_end
f_star = cvx_optval;

% constant_length_step
opts = struct();
opts.maxit = 3000;
opts.alpha0 = 0.0002;
opts.gamma = 0.005;
opts.step_type = 'constant_step_length';
[x, out] = l1_subgrad(x0, A, b, mu, opts);
data1 = (out.f_hist_best - f_star) / f_star;

% heavy ball acceleration
opts.step_type = 'constant_step_length';
[x, out] = l1_heavy_ball(x0, A, b, mu, opts);
data2 = (out.f_hist_best - f_star) / f_star;

% nesterov's acceleration
opts.step_type = 'constant_step_length';
[x, out] = l1_nesterov_acceleration(x0, A, b, mu, opts);
data3 = (out.f_hist_best - f_star) / f_star;

fig = figure;
semilogy(0:length(data1)-1, data1, '-', 'Color', [0.2 0.1 0.99], 'LineWidth', 1.5);
hold on
semilogy(0:length(data2)-1, data2, '-.', 'Color', [0.99 0.1 0.2], 'LineWidth', 1.5);
hold on
semilogy(0:length(data3)-1, data3, '--', 'Color', [0.1 0.99 0.2], 'LineWidth', 1.5);
legend('$\alpha_k = 0.001/||g_k||_2$', '$\alpha_k = 0.001/||g_k||_2$ with Heavy ball acceleration', '$\alpha_k = 0.001/||g_k||_2$ with Nesterov acceleration', 'interpreter', 'latex');
ylabel('$(\hat{f}(x^k) - f^*)/f^*$', 'fontsize', 14, 'interpreter', 'latex');
xlabel('number of iterations');
print(fig, '-depsc','f_accelerate.eps');
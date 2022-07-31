% L1-regularized least-squares example
randn('seed', 0);
rand('seed',0);

m = 1500;       % number of examples
n = 5000;       % number of features
p = 100/n;      % sparsity density

x0 = sprandn(n,1,p);
A = randn(m,n);
A = A*spdiags(1./sqrt(sum(A.^2))',0,n,n); % normalize columns
b = A*x0 + sqrt(0.001)*randn(m,1);

lambda_max = norm( A'*b, 'inf' );
lambda = 0.1*lambda_max;

[x history] = lasso(A, b, lambda, 1.0, 1.0, 1e-4, 1e-2);
[x0 history0] = lasso(A, b, lambda, 1.0, 1.0, 1e-5, 1e-3);
[x1 history1] = lasso(A, b, lambda, 1.0, 1.0, 1e-3, 1e-1);
K = length(history.objval);

h = figure;
plot(history0.time, history0.objval, 'MarkerSize', 10, 'LineWidth', 2, history.time, history.objval, 'MarkerSize', 10, 'LineWidth', 2, history1.time, history1.objval, 'MarkerSize', 10, 'LineWidth', 2);
ylabel('Objective function value'); xlabel('CPU time');

g = figure;
subplot(2,1,1);
semilogy(history0.time, max(1e-8, history0.r_norm), 'LineWidth', 2,history.time, max(1e-8, history.r_norm), 'LineWidth', 2, history1.time, max(1e-8, history1.r_norm), 'LineWidth', 2, ...
    history0.time, history0.eps_pri,'b--', 'LineWidth', 2, history.time, history.eps_pri,'r--', 'LineWidth', 2, history1.time, history1.eps_pri,'y--', 'LineWidth', 2);
ylabel('Primal residual'); xlabel('CPU time');

subplot(2,1,2);
semilogy(history0.time, max(1e-8, history0.s_norm), 'LineWidth', 2, history.time, max(1e-8, history.s_norm), 'LineWidth', 2,history1.time, max(1e-8, history1.s_norm), 'LineWidth', 2,...
    history0.time, history0.eps_dual,'b--', 'LineWidth', 2, history.time, history.eps_dual,'r--',  'LineWidth', 2, history1.time, history1.eps_dual,'y--',  'LineWidth', 2);
ylabel('Dual residual'); xlabel('CPU time');

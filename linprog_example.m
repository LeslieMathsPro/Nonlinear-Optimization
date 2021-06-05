randn('state', 0);
rand('state', 0);

n = 500;  % dimension of x
m = 400;  % number of equality constraints

c  = rand(n,1) + 0.5;    % create nonnegative price vector with mean 1
x0 = abs(randn(n,1));    % create random solution vector

A = abs(randn(m,n));     % create random, nonnegative matrix A
b = A*x0;
[x history] = linprog(c, A, b, 0.2, 1.0);
[x1 history1] = linprog(c, A, b, 1, 1.0);
[x2 history2] = linprog(c, A, b, 5, 1.0);
K = length(history.objval);
K1= length(history1.objval);
K2= length(history2.objval);

h = figure;
plot(history2.time, history2.objval,'LineWidth', 2, history1.time, history1.objval,'LineWidth', 2, history.time, history.objval,'LineWidth', 2);
ylabel('Objective function value'); xlabel('CPU time');
legend('\rho = 5', '\rho = 1', '\rho = 0.2' );

g = figure;
subplot(2,1,1);
semilogy(1:K, max(1e-8, history.r_norm), 'k', ...
    1:K, history.eps_pri, 'k--',  'LineWidth', 2);
ylabel('||r||_2');

subplot(2,1,2);
semilogy(1:K, max(1e-8, history.s_norm), 'k', ...
    1:K, history.eps_dual, 'k--', 'LineWidth', 2);
ylabel('||s||_2'); xlabel('iter (k)');
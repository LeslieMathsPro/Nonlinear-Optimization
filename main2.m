clc;
clear;

syms x1 x2;

% Drawing contour lines of the rosenbrock function near the minimum
[xgrid, ygrid] = meshgrid(-1.2:0.1:1.2, -0.5:0.1:1.5);
row = size(xgrid,1); col = size(xgrid,2);
for i = 1:row 
    for j = 1:col 
        x = [xgrid(i,j); ygrid(i,j)];
        func(i,j) = rosenbrock(x);
    end
end
figure;
contour(xgrid, ygrid,func,30); title('Contour Lines of Rosenbrock Function')

f = 10*(x2 - x1^2)^2 + (1 - x1)^2;
f = matlabFunction(f);
hold on
f_sym = sym(f);
F_td = matlabFunction(gradient(f_sym));
x0 = [-1.2 1]';
eps = 1e-5;
k = 0;
hist = [x0];

f_td = F_td(x0(1), x0(2));
	while norm(f_td) > eps
		d = - f_td;
		[newalpha] = linesearch_strongwolfe(x0, d, f, F_td);
		xnew = x0 + newalpha*d;
		hist = [hist xnew];
		refresh
		x0 = xnew;
		f_td = F_td(x0(1),x0(2));
		k = k + 1;
		
	end
	
for i = 1:length(hist) % no of iterations actually run
    funchist(i) = rosenbrock(hist(:,i));
end
figure;
semilogy(1:length(hist), funchist) % takes log of y-axis alone
xlabel('Iterations'); ylabel('Error of |f - f*|');
title('Convergence Evaluation')
	
% plotting the convergence on the contour plot ... 
figure;
contour(xgrid, ygrid,func,30); grid;title('Contour Lines of Rosenbrock Function')
hold on;
plot(hist(1,:), hist(2,:), 'ko-')
	
x = x0
result = f(x(1),x(2))

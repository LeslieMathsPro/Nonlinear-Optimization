% Main Script

clear all; close all;

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

% Convergence properties of BT steepest descent algorithm
x0 = [-1.2 1]'; % starting point of search
%[sol, hist, telapsed] = Linesearchcover(@rosenbrock, @rosenGradient, x0);
[sol, hist, telapsed] = Linesearchcover3(@rosenbrock, @rosenGradient, x0);

for i = 1:length(hist) % no of iterations actually run
    funchist(i) = rosenbrock(hist(:,i));
end
semilogy(1:length(hist), funchist) % takes log of y-axis alone
xlabel('Iterations'); ylabel('Error of |f - f*|');
title('Convergence Evaluation')

% plotting the convergence on the contour plot ... 
figure;
contour(xgrid, ygrid,func,30); grid;title('Contour Lines of Rosenbrock Function')
hold on;
plot(hist(1,:), hist(2,:), 'ko-')

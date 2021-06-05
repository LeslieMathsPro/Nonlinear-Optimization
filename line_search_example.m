F = @(x,y) 3*x.^2 + y.^4;
GF = @(x,y) [6*x ; 4*y.^3];

s = 1;
sigma = 0.1;
beta = 0.1;

x_k = 1.
y_k = -2.

function alpha = armijo(s, sigma, beta, x_k, y_k, F, GF, d, maxit = 10)
	m = 0;
	while m < maxit
		fprintf('Armijo rule beta^m is %5.2f, sigma is %5.2f, s is %4.2f\n',beta^m, sigma, s);
		f1 = F(x_k,y_k);
		f2 = F(x_k+beta^m*s*d(1),y_k+beta^m*s*d(2));
		diff = f1 - f2;
		lower_bound = -sigma*beta^m*s*GF(x_k,y_k)'*d;
		if diff >= lower_bound
			alpha = beta^m * s;
			break
		end
		m += 1;
	end
	alpha = beta^m*s;
end

function x_k1= steepest_descent(x_k, alpha_k, GF_k)
	x_k1 = x_k - alpha_k*GF_k;
end

d_k = -GF(x_k , y_k) / norm(GF(x_k,y_k));
armijo_stepsize = armijo(s, sigma, beta, x_k, y_k, F, GF, d_k);
fprintf('Armijo rule stepsize is %8.3f\n', armijo_stepsize);

out = steepest_descent([x_k,y_k], armijo_stepsize, -d_k);
fprintf('Steepest descent x_{k+1}: %8.3f y_{k+1}: %8.3f\n', out(1), out(2));

x_k1 = out(1);
y_k1 = out(2);
fprintf('Cost : %8.5f\n', F(x_k1,y_k1));

x = linspace(-2,2,100);
y = linspace(-2,2,100);
[X,Y] = meshgrid(x,y); %make the meshgrid
Z = 3*X.^2 + Y.^4;
surf(X,Y,Z);title('Steepest descent with Armijo rule f(x,y)=3x^2+y^4');
xlabel('x');
ylabel('y');
hold on
plot3([x_k,x_k1], [y_k,y_k1], [F(x_k,y_k),F(x_k1,y_k1)], 'r', 'Linewidth' , 3, 'MarkerSize', 2);
plot3(x_k, y_k, F(x_k, y_k), 'b', 'MarkerSize',10);
t = text(1,-2,22,'x_0 = 1.0, y_0 = -2.0, z_0 = 19.0','FontSize',15);
t = text(2,0,0.5,'x_{1} = 0.816, y_{1} = -1.017, z_{1} = 3.066','FontSize',15);
hold off
function func = rosenbrock(x)
	func = 10*(x(2) - x(1)^2)^2 + (1 - x(1))^2;
	%where x is a vector of variables, e.g. x = [x(1) x(2)];
	%function = 10(x_2 - x_1 ^2)^2 + (1 - x(1))^2;
end
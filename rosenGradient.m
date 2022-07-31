function grad = rosenGradient(x)
	grad = [2*x(1) - 40*x(1)*(x(2) - x(1)^2) - 2;...
			20*(x(2) - x(1)^2)];
end
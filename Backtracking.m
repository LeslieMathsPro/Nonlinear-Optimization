%Stepsize: 1. Backtracking 2. Armijo 3. Weak Wolfe 4. Wolfe 5. Strong Wolfe conditions
%Direction: 1. Steepest 2. Newton 3. Quasi-Newton

%in order to compare the effects of different stepsize method, we begin from 1. Steepest with different stepsize method
%1.a Backtracking + Steepest
%From the paper, typically one takes \gamma \in [0.5, 0.8], c \in [0.001, 1] with adjustments depending on the cost of function evaluation and degree of nonlinearity
function [sol, hist, time] = Backtracking(func, grad, x0)
	tic;
	gamma = 0.6;
	c = 0.01;
	tol = 10^(-5);
	t = 1;
	k = 0;
	x = x0;
	hist = [x0];
	while (norm(grad(x)) > tol) %convergence conditions 
		pk = -grad(x); %choose direction
		
		while (func(x+t*pk) > func(x) + t.*grad(x).*pk)
			t = gamma*t;
		end
		
		x = x + pk*t
		hist = [hist x];
		k = k + 1;
	end
	sol = x;
	time = toc;
end

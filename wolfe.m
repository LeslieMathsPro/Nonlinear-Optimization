function [newalpha] = wolfe(x, d, f0, f1)
	alpha = 1;
	sigma_1 = 0.25;
	sigma_2 = 0.75;
	
	a = 0;
	b = 1e6;
	while (1)
		x_new = x + d.*alpha;
		if ~(f0(x_new(1),x_new(2)) <= f0(x(1),x(2)) + sigma_1*alpha.*f1(x(1),x(2))'*d)
			b = alpha;
			alpha = (alpha + a)/2;
			continue;
		end
		
		if ~(f1(x_new(1),x_new(2))'*d >= sigma_2.*f1(x(1),x(2))'*d)
			a = alpha;
			alpha = min([2*alpha, (b+alpha)/2]);
			continue;
		end
	break
	end
	newalpha = alpha
end
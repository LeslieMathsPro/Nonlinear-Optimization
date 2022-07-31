function alphas = linesearch_strongwolfe(x0, d, f0, f1)
% function alphas = strongwolfe(f,d,x0,alpham)
% Line search algorithm satisfying strong Wolfe conditions
% Algorithms 3.5 on pages 60-61 in Nocedal and Wright

alpha0 = 0;
alphap = alpha0;
c1 = 0.25;
c2 = 0.75;
alpham = 1;
alphax = alpham*rand(1);
fx0 = f0(x0(1),x0(2));
gx0 = f1(x0(1),x0(2))'*d;
fxp = fx0;
gxp = gx0;
i=1;
% alphap is alpha_{i-1}
% alphax is alpha_i
while (1 ~= 2)
  xx = x0 + d.*alphax;
  fxx = f0(xx(1),xx(2));
  gxx = f1(xx(1),xx(2))'*d;
  if (fxx > fx0 + c1*alphax*gx0) || ((i > 1) && (fxx >= fxp)),
    alphas = zoom(f0,f1,x0,d,alphap,alphax)
    return;
  end
  if abs(gxx) <= -c2*gx0,
    alphas = alphax;
    return;
  end
  if gxx >= 0,
    alphas = zoom(f0,f1,x0,d,alphax,alphap);
    return;
  end
  alphap = alphax;
  fxp = fxx;
  gxp = gxx;
  alphax = alphax + (alpham-alphax)*rand(1);
  i = i+1;
end
function [x, out] = l1_subgrad(x0, A, b, mu, opts)

if ~isfield(opts, 'step_type'); opts.step_type = 'constant_step_size'; end
if ~isfield(opts, 'thres'); opts.thres = 1e-4; end 
if ~isfield(opts, 'alpha0'); opts.alpha0 = 0.01; end
if ~isfield(opts, 'maxit'); opts.maxit = 2000; end
if ~isfield(opts, 'ftol'); opts.ftol = 0; end

x = x0;
xprev = zeros(length(x),1);
out = struct();
out.f_hist = zeros(1, opts.maxit);
out.f_hist_best = zeros(1, opts.maxit);
out.g_hist = zeros(1, opts.maxit);
f_best = inf;

for k = 1:opts.maxit
    
    r = A * x - b;
    g = A' * r;
    out.g_hist(k) = norm(r, 2);
    
    f_now = 0.5 * norm(r, 2)^2 + mu * norm(x, 1);
    
    f_best = min(f_best, f_now);
    out.f_hist_best(k) = f_best;
    
    if k > 1 && abs(out.f_hist_best(k) - out.f_hist_best(k-1)) / abs(out.f_hist_best(1)) < opts.ftol
        break;
    end
    
    x(abs(x) < opts.thres) = 0;
    sub_g = g + mu * sign(x);
    
    alpha = set_step(k, opts, sub_g, f_now);
    xnew = x - alpha * sub_g + (x - xprev) * k / (k + 3);
    xprev = x;
    x = xnew;
end

out.itr = k;
out.f_hist = out.f_hist(1:k);
out.f_hist_best = out.f_hist_best(1:k);
out.g_hist = out.g_hist(1:k);
end
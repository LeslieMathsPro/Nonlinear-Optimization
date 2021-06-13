function a = set_step(k, opts, sub_g, f_now)
type = opts.step_type;
if strcmp(type, 'constant_step_size')
    a = opts.alpha0;
elseif strcmp(type, 'constant_step_length')
    a = opts.gamma / norm(sub_g,2);
elseif strcmp(type, 'square_summable_but_not_summable')
    a = opts.ssbns_a / (opts.ssbns_b + k);
elseif strcmp(type, 'nonsummable_diminishing')
    a = opts.nd_a / sqrt(k);
elseif strcmp(type, 'nonsummable_diminishing_step_lengths')
    a = opts.ndsl_a / (k * norm(sub_g,2));
elseif strcmp(type, 'polyak_step_size')
    a = (f_now - opts.f_star) / norm(sub_g,2)^2;
else
    error('unsupported type.');
end
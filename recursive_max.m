function f = recursive_max(v)
    if length(v) <= 1
        f = v;
    else
        f = v(1);
        if f < recursive_max(v(2:end));
            f = recursive_max(v(2:end));
        end
    end
end
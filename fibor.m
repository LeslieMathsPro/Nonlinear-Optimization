function f = fibor(n)
    if n <= 2
        f = ones(1,n);
    else
        u = fibor(n-1);
        f = [u,sum(u(end-1:end))];
    end
end
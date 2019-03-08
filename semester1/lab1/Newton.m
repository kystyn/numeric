% f - function
% df - functional derivative
% b - start point (!!!)
% eps - precision
function [x, steps, appr_roots] = Newton( f, df, a, b, eps )
    steps = 1;
    x = b;
    appr_roots(1) = x;
    while (f(x) * f(x - eps) > 0)
        steps = steps + 1;
        x = x - f(x) / df(x);
        appr_roots(steps) = x;
    end
end
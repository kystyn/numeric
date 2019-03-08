% f - function
% f(a) * f(b) < 0
% eps - precision
% root - root
% strps - count of steps
% appr_roots
function [root, steps, appr_roots] = BinomialDivision( f, df, a, b, eps )
    steps = 0;
    l = a;
    r = b;
    while (r - l > 2 * eps)
        steps = steps + 1;
        m = (l + r) / 2.0;
        appr_roots(steps) = m;
        if (f(l) * f(m) < 0)
            r = m;
        else
            if (f(l) * f(m) > 0)
                l = m;
            else
                root = m;
                return;
            end
        end
    end
    root = (l + r) / 2.0;
end

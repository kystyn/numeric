% Method - method function reference
% MethodName 
% f - equation function reference
% df - derivative equation function reference
% EqFormula
% a, b - borders
function GeomInterpretNewton( f, df, EqFormula, a, b, style )
    eps = 1e-6;
    X = a : 0.01 : b;
    
    %%% draw function
    for i = 1 : 1 : length(X)
        Y(i) = f(X(i));
    end
    
    hold on;
    grid on;
    
    plot(X, Y, style);
    
    %%% draw geometrical interpretation
    [root, steps, appr_roots] = Newton(f, df, a, b, eps);
    
    for i = 1 : 1 :min(3,  length(appr_roots))
        YA(1) = f(appr_roots(i)) + df(appr_roots(i)) * (b - appr_roots(i));
        YA(2) = f(appr_roots(i)) + df(appr_roots(i)) * (a - appr_roots(i));
        XA(1) = b;
        XA(2) = a;
        if (appr_roots(i) >= a && appr_roots(i) <= b)
          plot(XA, YA, 'b');
          plot(appr_roots(i), f(appr_roots(i)), 'go');
          text(appr_roots(i), f(b) / 2, sprintf('%i', i));
        end
    end
    plot(root, 0, 'black*', 'MarkerSize', 10);
    xlabel('x');
    ylabel('y');
    title(sprintf('Geometric interpretation.\nMethod: Newton.\nFunction: %s', EqFormula));
end
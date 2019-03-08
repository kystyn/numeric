% Method - method function reference
% MethodName 
% f - equation function reference
% df - derivative equation function reference
% EqFormula
% a, b - borders
function GeomInterpretBinDiv( f, df, EqFormula, a, b, style )
    eps = 1e-6;
    X = a : 0.01 : b;
    
    %%% draw function
    Y = f(X);
    
    hold on;
    grid on;
    
    plot(X, Y, style);
    
    %%% draw geometrical interpretation
    [root, steps, appr_roots] = BinomialDivision(f, df, a, b, eps);
    
    YA(1) = f(a);
    YA(2) = f(b);
    for i = 1 : 1 : min(5, length(appr_roots))
        XA(1) = appr_roots(i);
        XA(2) = appr_roots(i);
        plot(XA, YA, 'b');
        text(XA(1), f(a) / 2.0, sprintf('%i', i));
    end
    plot(root, 0, 'black*', 'MarkerSize', 10);
    xlabel('x');
    ylabel('y');
    title(sprintf('Geometric interpretation.\nMethod: bin. div.\nFunction: %s', EqFormula));
end
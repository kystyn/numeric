% Method - method function reference
% MethodName 
% Equation - equation function reference
% dEquation - derivative equation function reference
% EqFormula
% a, b - borders
function PrecisionGraphic( Method, Equation, dEquation, EqFormula, a, b, style )
    eps = 1e-1;
    i = 1;
    ln10 = log(10);
    while (eps >= 1e-10)
        [root, steps(i)] = Method(Equation, dEquation, a, b, eps);
        LogX(i) = log(eps) / ln10;
        X(i) = eps;
        eps = eps / 2.0;
        i = i + 1;
    end
    
    % draw plot
    semilogx(X, steps, style);
    xlabel('log_{10}(Epsilon)');
    ylabel('Number of steps');
    title(sprintf('Dependence of steps count on epsilon.\nFunction %s', EqFormula));       
end


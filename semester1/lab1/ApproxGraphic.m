% Method - numeric method
% f - equation function reference
% df - derivative equation function reference
% EqFormula
% a, b - start approximation
function ApproxGraphic( Method, f, df, EqFormula, a, b, style )    
    eps = 1e-8;
    [root, steps(1)] = Method(f, df, a, b, eps);
    offa = root - a;
    offb = b - root;
   
    i = 2;
    ln2 = log(2);
    
    logoff(1) = log(b - a) / ln2;
    off(1) = b - a;
    while (i < 10)
        a = a + offa / 2;
        b = b - offb / 2;
        
        offa = offa / 2;
        offb = offb / 2;
        
        logoff(i) = log(b - a) / ln2;
        off(i) = b - root;
        
        [root, steps(i)] = Method(f, df, a, b, eps);
        i = i + 1;
    end
   
    % draw plot
    semilogx(off, steps,  style);
    xlabel('log_{10}(Interval size)');
    ylabel('Number of steps');
    title(sprintf('Dependence of steps count on approximation.\nFunction %s', EqFormula));       
end


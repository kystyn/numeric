function main()
    eps = 1e-5;
    
    % Functions table 
    fprintf('Function\t\t\t\tMethod\tRoot\t\tEpsilon\tSteps\n');
    
    %%% Useful values
    f1_formula = '2x^3 - 9x^2 - 60x - 1 = 0';
    f1_b = 1 + 60 / 2;
    f1_a = 1 / f1_b;
    
    f2_formula = '2^x * x = 1';
    f2_a = 0;
    f2_b = 1;
    
    %%% Fill table 
    [x, st] = BinomialDivision(@f1, @df1, f1_a, f1_b, eps);
    fprintf('%s\tBinomial division\t%1.8f\t%g\t\t%i\n', f1_formula, x, eps, st);
    
    [x, st] = Newton(@f1, @df1, 8, f1_b, eps);
    fprintf('%s\t\tNewton\t%1.8f\t%g\t\t%i\n', f1_formula, x, eps, st);
    
    [x, fv, ex, out]  = fzero(@f1, 8);
    fprintf('%s\t\tfzero\t\t%1.8f\t%g\t\t%i\n', f1_formula, x, eps, out.iterations);
       
    [x, st] = BinomialDivision(@f2, @df2, f2_a, f2_b, eps);
    fprintf('%s\t\t\tBinomial division\t%1.8f\t%g\t\t%i\n', f2_formula, x, eps, st);
    
    [x, st] = Newton(@f2, @df2, f2_a, f2_b, eps);
    fprintf('%s\t\t\t\tNewton\t%1.8f\t%g\t\t%i\n', f2_formula, x, eps, st);
    
    [x, fv, ex, out] = fzero(@f2, f2_a);
    fprintf('%s\t\t\t\tfzero\t\t%1.8f\t%g\t\t%i\n', f2_formula, x, eps, out.iterations);
    
    %%%% Depndencies research
    % On epsilon
    figure;
    
    subplot(f2_b, 2, f2_b);
    PrecisionGraphic(@BinomialDivision, @f1, @df1, f1_formula, f1_a, f1_b, 'r:');
    hold on;
    PrecisionGraphic(@Newton, @f1, @df1, f1_formula, 8, f1_b, 'g-');
    legend('Binomial division', 'Newton');
    
    subplot(f2_b, 2, 2);
    PrecisionGraphic(@BinomialDivision, @f2, @df2, f2_formula, f2_a, f2_b, 'r:');
    hold on;
    PrecisionGraphic(@Newton, @f2, @df2, f2_formula, f2_a, f2_b, 'g-');
    legend('Binomial division', 'Newton');
    
    % On interval size
    figure;
    
    subplot(1, 2, 1);
    ApproxGraphic(@BinomialDivision, @f1, @df1, f1_formula, f1_a, f1_b, 'r:');
    hold on;
    ApproxGraphic(@Newton, @f1, @df1, f1_formula, f1_a, f1_b, 'g-');
    legend('Binomial division', 'Newton', 'Location', 'west');
    
    subplot(1, 2, 2);
    ApproxGraphic(@BinomialDivision, @f2, @df2, f2_formula, f2_a, f2_b, 'r:');
    hold on;
    ApproxGraphic(@Newton, @f2, @df2, f2_formula, f2_a, f2_b, 'g-');
    legend('Binomial division', 'Newton', 'Location', 'west');
    
    % Geometric interpretation
    figure;
    subplot(2, 2, 1);
    GeomInterpretBinDiv(@f1, @df1, f1_formula, 7, 9, 'r');
    subplot(2, 2, 2);
    GeomInterpretBinDiv(@f2, @df2, f2_formula, f2_a, f2_b, 'r');
    
    subplot(2, 2, 3);
    GeomInterpretNewton(@f1, @df1, f1_formula, 7, 9, 'r');
    subplot(2, 2, 4);
    GeomInterpretNewton(@f2, @df2, f2_formula, f2_a, f2_b, 'r');
end

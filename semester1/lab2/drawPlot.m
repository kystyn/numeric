function drawPlot( X, t, axisx )
    N = 15;
    F = fopen('solve.out', 'r');         
    format long;
    x = fscanf(F, '%f');
    fclose(F);
    
    s = size(x);
    s = s(1);
   
    x0 = transpose(genVector(N, 1));
    Y = zeros(1, s / (4 * N) );
    a = 1;
    for i = 1 : 4 * N : s
        b0 = x(i : i + N - 1);
        b1 = x(i + N: i + 2 * N - 1);
        x0 = x(i + 2 * N: i + 3 * N - 1);
        x1 = x(i + 3 * N: i + 4 * N - 1);
        Y(a) = norm(x1 - x0) / norm(x0);
        X(a) = norm(x1 - x0) / norm(b1);
        a = a + 1;
    end
    
    figure;
    loglog(X, Y, 'r.');
    %hold on;
    %plot([min(X) max(X)], [min(X) max(X)] , 'g-');
    title(t);
    xlabel(axisx);
    ylabel('{||\delta{x}||}/{||x||}');
end  

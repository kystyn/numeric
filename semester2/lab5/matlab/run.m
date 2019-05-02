function run()
    format long;
    a = 1;
    b = 3;
    cp = exp(1);
    
    N = 8;
    tolls = zeros(1, N);
    drawX = zeros(1, N);
    drawY = zeros(2, N);
    drawFrag = zeros(2, N);
    drawMinFrag = zeros(2, N);
    drawMaxFrag = zeros(2, N);

    f = @(x) (exp(x) .* (log(x) + 1));
    %f = @(x) (x * x);
    for t = 1 : 1 : 2
        for toll = 1 : 1 : N
          F = fopen('de.in', 'w');
          prepare(F, t, a, b, cp, 10 ^ -toll);
          tolls(toll) = 10 ^ -toll;
          fclose(F);

          %system('lab5.exe');
          system('./lab5');
          [X, Y, frag, minfrag, maxfrag, fr] = fileParser('de.out');
          drawY(t, toll) = deviation(Y, X, f, fr);
          drawMinFrag(t, toll) = minfrag;
          drawMaxFrag(t, toll) = maxfrag;
          drawFrag(t, toll) = frag;
        end
    end
    
    figure;
    drawSmth(drawFrag(1, 1 : N), drawY(1, 1 : N), 'r', 'Dependence of fact precision on fragmentation', ...
            'Fragmentation', 'Fact precision');
    grid on;
    hold on;
    drawSmth(drawFrag(2, 1 : N), drawY(2, 1 : N), 'g', 'Dependence of fact precision on fragmentation', ...
            'Fragmentation', 'Fact precision');
    legend('E', 'I');

    figure;
    drawSmth(tolls, drawMinFrag(1, 1 : N), 'r', '', ...
    'Fact precision', 'Step count');
    grid on;
    hold on;
    drawSmth(tolls, drawMinFrag(2, 1 : N), 'g', '', ...
    'Fact precision', 'Step count');
    drawSmth(tolls, drawMaxFrag(1, 1 : N), 'b', 'Dependence of min/max frag on precision', ...
    'Fact precision', 'Step count');
    drawSmth(tolls, drawMaxFrag(2, 1 : N), 'y', 'Dependence of min/max frag on precision', ...
    'Fact precision', 'Step count');
    legend('min E', 'min I', 'max E', 'max I');
    hold off;
    
    toll = 4;
    for t = 1 : 1 : 2
        i = 1;
        for error = -N : 1 : -1
          F = fopen('de.in', 'w');
          prepare(F, t, a, b, cp * (1 + 10 ^ error), 10 ^ -toll);
          fclose(F);

          %system('lab5.exe');
          system('./lab5');
          [X, Y, frag, minfrag, maxfrag, fr] = fileParser('de.out');
          drawX(i) = 10 ^ error;
          drawY(t, i) = deviation(Y, X, f, fr);
          i = i + 1;
        end
    end

    figure;
    drawSmth(drawX, drawY(1, 1 : N), 'r', 'Dependence of fact precision on Cauchy problem error', ...
    'Error', 'Fact precision');
    hold on;
    drawSmth(drawX, drawY(2, 1 : N), 'g', 'Dependence of fact precision on Cauchy problem error', ...
    'Error', 'Fact precision');
    grid on;
    legend('E', 'I');
end
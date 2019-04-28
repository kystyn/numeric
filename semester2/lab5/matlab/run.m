function run()
    format long;
    a = 1;
    b = 3;
    cp = exp(1);
    
    N = 8;
    tolls = zeros(1, N);
    drawX = zeros(1, N);
    drawY = zeros(1, N);
    drawMinFrag = zeros(1, N);
    drawMaxFrag = zeros(1, N);

    f = @(x) (exp(x) .* (log(x) + 1));
    %f = @(x) (x * x);
    for toll = 1 : 1 : N
      F = fopen('de.in', 'w');
      prepare(F, a, b, cp, 10 ^ -toll);
      tolls(toll) = 10 ^ -toll;
      fclose(F);

      system('lab5.exe');
      %system('./lab5');
      [X, Y, minfrag, maxfrag, N] = fileParser('de.out');
      drawY(toll) = deviation(Y, X, f, N);
      drawMinFrag(toll) = minfrag;
      drawMaxFrag(toll) = maxfrag;
    end
    
    figure;
      drawSmth(tolls, drawY, 'r', 'Dependence of fact precision on settable', ...
      'Settable precision', 'Fact precision');
      grid on;

      figure;
      drawSmth(tolls, drawMinFrag, 'r', '', ...
      'Fact precision', 'Step count');
      hold on;
      drawSmth(tolls, drawMaxFrag, 'b', 'Dependence of min/max frag on precision', ...
      'Fact precision', 'Step count');
      grid on;
      legend('min', 'max');
      hold off;
    
    toll = 6;
    i = 1;
    for error = -N : 1 : -1
      F = fopen('de.in', 'w');
      prepare(F, a, b, cp * (1 + 10 ^ error), 10 ^ -toll);
      fclose(F);

      system('lab5.exe');
      %system('./lab5');
      [X, Y, minfrag, maxfrag, N] = fileParser('de.out');
      drawX(i) = 10 ^ error;
      drawY(i) = deviation(Y, X, f, N);
      i = i + 1;
    end

    figure;
    drawSmth(drawX, drawY, 'r', 'Dependence of fact precision on Cauchy problem error', ...
    'Error', 'Fact precision');
    grid on;
end
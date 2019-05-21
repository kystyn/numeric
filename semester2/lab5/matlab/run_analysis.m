function run_analysis()
format long;
    a = 0.2;
    b = 1;
    cp = [5; -1/0.04];
    a0 = 1;
    a1 = 0;
    A = 5;
    b0 = 1;
    b1 = 0;
    B = 1;
    
    N = 8;
    tolls = zeros(1, N);
    drawX = zeros(3, N);
    drawY = zeros(3, N);
    drawFrag = zeros(3, N);
    drawErrLocal = zeros(3, N);
    drawVolume = zeros(3, N);

    f = @(x) (1 / x);
    %f = @(x) (x * x);
    for toll = 1 : 1 : N
      F = fopen('de.in', 'w');
      prepare(F, 0, a, b, cp, 10 ^ -toll);
      tolls(toll) = 10 ^ -toll;
      fclose(F);

      %system('lab5.exe cauchy');
      system('./lab5 cauchy');
      [X, Y, frag, minfrag, maxfrag, maxerrlocal, volume, fr] = fileParser('de.out');
      [drawY(1, toll), drawErrLocal(1, toll)] = deviation(Y, X, f, fr);
      drawFrag(1, toll) = frag;
      drawVolume(1, toll) = volume;
    end
    
    for t = 1 : 1 : 2
        for toll = 1 : 1 : N
          F = fopen('de.in', 'w');
          prepare_boundary(F, t, a, b, a0, a1, A, b0, b1, B, 10 ^ -toll);
          tolls(toll) = 10 ^ -toll;
          fclose(F);

          %system('lab5.exe boundary');
          system('./lab5 boundary');
          [X, Y, frag, minfrag, maxfrag, maxerrlocal, volume, fr] = fileParserBoundary('de.out');
          [drawY(t + 1, toll), drawErrLocal(t + 1, toll)] = deviation(Y, X, f, fr);
          drawFrag(t + 1, toll) = frag;
          drawVolume(t + 1, toll) = volume;
        end
    end
    
    figure;
    drawSmth(drawFrag(1, 1 : N), drawY(1, 1 : N), 'r', 'Dependence of fact precision on fragmentation', ...
            'Fragmentation', 'Fact precision');
    grid on;
    hold on;
    drawSmth(drawFrag(2, 1 : N), drawY(2, 1 : N), 'g', 'Dependence of fact precision on fragmentation', ...
            'Fragmentation', 'Fact precision');
    drawSmth(drawFrag(3, 1 : N), drawY(3, 1 : N), 'b', 'Dependence of fact precision on fragmentation', ...
            'Fragmentation', 'Fact precision');
    legend('Euler-cauchy', 'Finite diff', 'Reductor');

    figure;
    drawSmth(drawY(1, 1 : N), drawVolume(1, 1 : N), 'r', '', ...
    'Fact precision', 'Eval volume');
    grid on;
    hold on;
    drawSmth(drawY(2, 1 : N), drawVolume(2, 1 : N), 'g', '', ...
    'Fact precision', 'Eval volume');
    drawSmth(drawY(3, 1 : N), drawVolume(3, 1 : N), 'b', 'Dependence of evaluation volume on fact precision', ...
    'Fact precision', 'Eval volume');
    legend('Euler-cauchy', 'Finite diff', 'Reductor');
    hold off;
    
    figure;
    drawSmthlogXY(tolls, drawErrLocal(1, 1 : N), 'r', '', ...
    'Fact precision', 'Eval volume');
    grid on;
    hold on;
    drawSmthlogXY(tolls, drawErrLocal(2, 1 : N), 'g', '', ...
    'Fact precision', 'Eval volume');
    drawSmthlogXY(tolls, drawErrLocal(3, 1 : N), 'b', 'Dependence of max error localization on fact precision', ...
    'Fact precision', 'Eval volume');
    legend('Euler-cauchy', 'Finite diff', 'Reductor');
    
    toll = 4;
    i = 1;
    for error = -N : 1 : -1
      F = fopen('de.in', 'w');
      prepare(F, 0, a, b, cp .* (1 + 10 ^ error), 10 ^ -toll);
      fclose(F);

      %system('lab5.exe cauchy');
      system('./lab5 cauchy');
      [X, Y, frag, minfrag, maxfrag, maxerrlocal, volume, fr] = fileParser('de.out');
      drawX(1, i) = 10 ^ error;
      drawY(1, i) = deviation(Y, X, f, fr);
      i = i + 1;
    end
    
    for t = 1 : 1 : 2
        i = 1;
        for error = -N : 1 : -1
          F = fopen('de.in', 'w');
          prepare_boundary(F, t, a, b, a0, a1, A * (1 + 10 ^ error), b0, b1, B * (1 - 10 ^ error), 10 ^ -toll);
          fclose(F);

          %system('lab5.exe boundary');
          system('./lab5 boundary');
          [X, Y, frag, minfrag, maxfrag, maxerrlocal, volume, fr] = fileParserBoundary('de.out');
          drawX(t + 1, i) = 10 ^ error;
          drawY(t + 1, i) = deviation(Y, X, f, fr);
          i = i + 1;
        end
    end

    figure;
    drawSmth(drawX(1, 1 : N), drawY(1, 1 : N), 'r', 'Dependence of fact precision on Cauchy problem error', ...
    'Error', 'Fact precision');
    hold on;
    drawSmth(drawX(2, 1 : N), drawY(2, 1 : N), 'g', 'Dependence of fact precision on Cauchy problem error', ...
    'Error', 'Fact precision');
    grid on;
    drawSmth(drawX(3, 1 : N), drawY(3, 1 : N), 'b', 'Dependence of fact precision on Cauchy/boundary problem error', ...
    'Error', 'Fact precision');
    grid on;
    legend('Euler-cauchy', 'Finite diff', 'Reductor');
end
function run()
  a = 0.1;
  b = 1.1;
  
  grids = zeros(1, 3);

  grids(1) = 'R';
  grids(2) = 'U';
  grids(3) = 'C';

  weights = zeros(1, 4);

  weights(1) = 'B';
  weights(2) = 'M';
  weights(3) = 'E';
  weights(4) = 'U'; 

  # ALWAYS change grid size

  # uniform grid
  # uniform weight
  # change poly degree
  styles = ['r'; 'g'; 'b'; 'r:'; 'g:'; 'b:';
            'r--'; 'g--'; 'b--'; 'r-.'; 'g-.'; 'b-.'; 'r-'; 'b-'; 'g-'; 'y'];
  for funcT = 0 : 1 : 1
    figure;
    for polyDegree = 1 : 1 : 16
      F = fopen('lsm.in', 'w');
      prepare(F, a, b, 'U', funcT, 'U', polyDegree);
      fclose(F);
      system('./lab2');
      [X, Y] = fileParser('lsm.out');
      drawSmth(X, Y, styles(polyDegree, 1 : 3), strcat('Uniform grid, uniform weight, change poly deg ', char('0' + funcT),'-diff. function'));
      hold on;
      grid on;
    end
    legend('1', '2', '3', 
           '4', '5', '6',
           '7', '8', '9',
           '10', '11', '12', 
           '13', '14', '15', '16');
  end
  
  # uniform grid
  # change weight
  # fixed degree
  for funcT = 0 : 1 : 1
    figure;
    for w = 1 : 1 : 4
      F = fopen('lsm.in', 'w');
      prepare(F, a, b, 'U', funcT, weights(w), 16);
      fclose(F);
      system('./lab2');
      [X, Y] = fileParser('lsm.out');
      drawSmth(X, Y, styles(w, 1 : 3), strcat('Uniform grid, change weight, fixed grid, ', char('0' + funcT), '-diff. function'));
      hold on;
      grid on;
    end
    legend('normal begin', 'normal mid', 'normal end', 'uniform');
  end
  
  # change grid
  # uniform weight
  # fixed degree
  for funcT = 0 : 1 : 1
    figure;
    for g = 1 : 1 : 3
      F = fopen('lsm.in', 'w');
      prepare(F, a, b, grids(g), funcT, 'U', 16);
      fclose(F);
      system('./lab2');
      [X, Y] = fileParser('lsm.out');
      drawSmth(X, Y, styles(g, 1 : 3), strcat('Change grid, uniform weight, fixed grid, ', char('0' + funcT), '-diff. function'));
      hold on;
      grid on;
    end
    legend('random', 'uniform', 'chebyshev');
  end
end
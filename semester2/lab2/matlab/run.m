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
  for funcT = 0 : 1 : 1
    figure;
    for polyDegree = 1 : 1 : 12
      F = fopen('lsm.in', 'w');
      prepare(F, a, b, 'U', funcT, 'U', polyDegree);
      fclose(F);
      system('./lab1');
      drawSmth(fileParser('lsm.out'), 'r', 'Uniform grid, uniform weight, change poly deg ' + char(funcT) + '-diff. function');
      hold on;
    end
  end
  
  # uniform grid
  # change weight
  # fixed degree
  for funcT = 0 : 1 : 1
    figure;
    for w = 1 : 1 : 4
      F = fopen('lsm.in', 'w');
      prepare(F, a, b, 'U', funcT, weights(w), 10);
      fclose(F);
      system('./lab1');
      drawSmth(fileParser('lsm.out'), 'r', 'Uniform grid, change weight, fixed grid, ' + char(funcT) + '-diff. function');
      hold on;
    end
  end
  
  # change grid
  # uniform weight
  # fixed degree
  for funcT = 0 : 1 : 1
    figure;
    for g = 1 : 1 : 3
      F = fopen('lsm.in', 'w');
      prepare(F, a, b, grids(g), funcT, 'U', 10);
      fclose(F);
      system('./lab1');
      drawSmth(fileParser('lsm.out'), 'r', 'Change grid, uniform weight, fixed grid, ' + char(funcT) + '-diff. function');
      hold on;
    end
  end
end
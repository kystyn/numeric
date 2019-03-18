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
    for polyDegree = 1 : 1 : 12
      prepare(a, b, 'U', funcT, 'U', polyDegree);
      system('./lab1');
    end
  end
  
      
end
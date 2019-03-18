function N = prepare(a, b, gridT, funcT, weightT )
F = fopen('lsm.in', 'w');

for polyDim = 2 : 3 : 20
  for nodeCount = 2 : 1 : 100
    fprintf(F, '%f %f %i %c %c %i %i\n', a, b, ...
      nodeCount, gridT, weightT, funcT, polyDim);
    end
end

fclose(F);
end
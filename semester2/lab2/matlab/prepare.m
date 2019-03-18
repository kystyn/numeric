function prepare(F, a, b, gridT, funcT, weightT, polyDegree )
  for nodeCount = polyDegree + 1 : 1 : 100
    fprintf(F, '%f %f %i %c %c %i %i\n', a, b, ...
      nodeCount, gridT, weightT, funcT, polyDegree);
  end
end